/*
Kernels:
http://naif.jpl.nasa.gov/pub/naif/
*/
//////////////////////////////////////////
//HEADERS
//////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <SpiceUsr.h>
#include <novas.h>
#include <eph_manager.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_const_mksa.h>

//////////////////////////////////////////
//MACROS
//////////////////////////////////////////
#define D2R(x) (x*M_PI/180)
#define R2D(x) (x*180/M_PI)
#define POWI(x,n) gsl_pow_int(x,n)
#define SGN(x) (x<0?-1:+1)

//////////////////////////////////////////
//CSPICE CONSTANTS
//////////////////////////////////////////
#define EARTH_ID "EARTH"
#define ATTEMPTS 12 /*SEE NUMBER_OF_STEPS*/

//////////////////////////////////////////
//CONSTANTS
//////////////////////////////////////////
#define GCONST GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
#define GKMS (GCONST*1E-9)
#define MSUN 1.9885E30/*kg*/
#define YEAR (365.25*GSL_CONST_MKSA_DAY)
#define DAY GSL_CONST_MKSA_DAY
#define AU GSL_CONST_MKSA_ASTRONOMICAL_UNIT

//////////////////////////////////////////
//GLOBAL VARIABLES
//////////////////////////////////////////
static int NUMBER_OF_STEPS[]={2,4,6,8,12,16,24,32,48,64,96,128};
double REARTH;
double RPEARTH;
double FEARTH;
gsl_rng* RAND;
double GGLOBAL;
double UL,UM,UT,UV;

//////////////////////////////////////////
//ROUTINES
//////////////////////////////////////////
int initSpice(void)
{
  SpiceInt n;
  SpiceDouble radii[3];

  //KERNELS
  furnsh_c("kernels.txt");

  //EARTH RADII
  bodvrd_c(EARTH_ID,"RADII",3,&n,radii);
  REARTH=radii[0];
  RPEARTH=radii[0];
  FEARTH=(radii[0]-radii[2])/radii[0];

  //RANDOM NUMBERS
  RAND=gsl_rng_alloc(gsl_rng_default);

  return 0;
}

char* dec2sex(double dec)
{
  double d,m,s;
  int id,im,sgn;
  char *str=(char*)calloc(sizeof(char),100); 
  d=fabs(dec);
  sgn=dec/d;
  id=floor(d);
  m=(d-id)*60;
  im=floor(m);
  s=(m-im)*60;
  sprintf(str,"%+d:%02d:%.3f",sgn*id,im,s);
  return str;
}

double sex2dec(double d,double m,double s)
{
  double s2d;
  s2d=d+m/60.0+s/3600.0;
  return s2d;
}

char* vec2str(double vec[],char frm[]="%.8e")
{
  char format[100];
  char *str=(char*)calloc(sizeof(char),100); 
  sprintf(format,"%s %s %s",frm,frm,frm);
  sprintf(str,format,vec[0],vec[1],vec[2]);
  return str;
}

char* vec2strn(double vec[],int n,char frm[]="%.8e")
{
  int i;
  char format[100];
  char *str=(char*)calloc(sizeof(char),100*n);
  sprintf(format,"%ss%s","%",frm);
  for(i=0;i<n;i++) sprintf(str,format,str,vec[i]);
  return str;
}

double greatCircleDistance(double lam1,double lam2,
			   double phi1,double phi2)
{
  double d;

  //HARVESINE FORMULA
  double sf,sl;
  sf=sin((phi2-phi1)/2);
  sl=sin((lam2-lam1)/2);
  d=2*asin(sqrt(sf*sf+cos(phi1)*cos(phi2)*sl*sl));

  return d;
}

int bodyEphemerisApparent(ConstSpiceChar *body,
			  SpiceDouble t,
			  SpiceDouble lon,SpiceDouble lat,SpiceDouble alt,
			  SpiceDouble *range,
			  SpiceDouble *ltime,
			  SpiceDouble *raJ2000,
			  SpiceDouble *decJ2000,
			  SpiceDouble *ra,
			  SpiceDouble *dec
			  )
{
  SpiceDouble earthSSBJ2000[6];
  SpiceDouble bodyJ2000[6],bodySSBJ2000[6],ltbody;
  SpiceDouble bodyTOPOJ2000[3],bodyTOPOEpoch[3];
  SpiceDouble Dbody,RAbody,DECbody,RAbodyJ2000,DECbodyJ2000;
  SpiceDouble observerITRF93[3],observerJ2000[3],observerSSBJ2000[3];
  SpiceDouble M_J2000_Epoch[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  SpiceDouble M_ITRF93_J2000[3][3];
  SpiceDouble d,lt,ltmp,ltold,lttol=1E-2;
  int i,ie=0,ncn=10;
  double cspeed=clight_c();

  //ROTATION MATRIX AT THE TIME OF EPHEMERIS
  pxform_c("J2000","EARTHTRUEEPOCH",t,M_J2000_Epoch);
  pxform_c("ITRF93","J2000",t,M_ITRF93_J2000);

  //OBSERVER POSITION J2000 RELATIVE TO EARTH CENTER
  georec_c(D2R(lon),D2R(lat),alt/1000.0,REARTH,FEARTH,observerITRF93);
  mxv_c(M_ITRF93_J2000,observerITRF93,observerJ2000);

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //ASTROMETRIC POSITION
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  i=0;
  lt=0.0;ltold=1.0;
  spkezr_c(EARTH_ID,t,"J2000","NONE","SOLAR SYSTEM BARYCENTER",
	   earthSSBJ2000,&ltmp);
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);
  while((fabs(lt-ltold)/lt)>=lttol && i<ncn){
    ltold=lt;
    spkezr_c(body,t-lt,"J2000","NONE","SOLAR SYSTEM BARYCENTER",bodySSBJ2000,&ltmp);
    vsub_c(bodySSBJ2000,observerSSBJ2000,bodyTOPOJ2000);
    d=vnorm_c(bodyTOPOJ2000);
    lt=d/cspeed;
    i++;
  }
  recrad_c(bodyTOPOJ2000,&d,&RAbodyJ2000,&DECbodyJ2000);
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //CORRECTED POSITION
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //OBSERVER POSITION J2000
  spkezr_c(EARTH_ID,t,"J2000","NONE","EARTH BARYCENTER",
	   earthSSBJ2000,&lt);
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);

  //CORRECTED OBJECT POSITION
  spkezr_c(body,t,"J2000","LT+S","EARTH BARYCENTER",
	   bodySSBJ2000,&lt);
  vsub_c(bodySSBJ2000,observerSSBJ2000,bodyTOPOJ2000);

  //PRECESS POSITION
  mxv_c(M_J2000_Epoch,bodyTOPOJ2000,bodyTOPOEpoch);

  //RA & DEC PRECESSED
  recrad_c(bodyTOPOEpoch,&d,&RAbody,&DECbody);

  //RETURNED
  *range=d;
  *ltime=lt;
  *ra=RAbody*180/M_PI/15;
  *dec=DECbody*180/M_PI;
  *raJ2000=RAbodyJ2000*180/M_PI/15;
  *decJ2000=DECbodyJ2000*180/M_PI;
  return 0;
}

/*
  Calculates the julian date.  If et=0 it gives the Ephemeris Julian
  Date.  If et=1 it gives the Jul ian Date in the International Atomic
  time reference.
 */
SpiceDouble t2jd(SpiceDouble t,int et=0)
{
  SpiceDouble deltat;
  deltet_c(t,"ET",&deltat);
  double tjd=unitim_c(t,"ET","JED");
  return tjd-et*deltat/86400.0;
}

/*
  Matrix to convert from planetocentric celestial coordinates to
  horizontal coordinates and viceversa.

  See discussion at:
  https://naif.jpl.nasa.gov/pipermail/spice_discussion/2010-July/000307.html

  h2m: converts from geocentric to topocentric
  h2i: converts from topocentric to geocentric
 */
void hormat(SpiceDouble lat,SpiceDouble lon,SpiceDouble t,SpiceDouble h2m[3][3],SpiceDouble h2i[3][3])
{
  SpiceDouble geopos[3],normal[3],normalJ2000[3],normalEpoch[3];
  SpiceDouble ux[]={1,0,0},uy[]={0,1,0},uz[]={0,0,1},uzJ2000[3],uzEpoch[3];
  SpiceDouble M_J2000_Epoch[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  SpiceDouble M_ITRF93_J2000[3][3];

  //TRANSFORM MATRICES
  pxform_c("J2000","EARTHTRUEEPOCH",t,M_J2000_Epoch);
  pxform_c("ITRF93","J2000",t,M_ITRF93_J2000);
  
  //NORMAL ITRF93
  georec_c(D2R(lon),D2R(lat),0.0,REARTH,FEARTH,geopos);
  surfnm_c(REARTH,REARTH,RPEARTH,geopos,normal);

  //NORMAL EPOCH
  mxv_c(M_ITRF93_J2000,normal,normalJ2000);
  mxv_c(M_J2000_Epoch,normalJ2000,normalEpoch);

  //Z EPOCH
  mxv_c(M_ITRF93_J2000,uz,uzJ2000);
  mxv_c(M_J2000_Epoch,uzJ2000,uzEpoch);

  //TRANSFORM MATRICES
  ucrss_c(normalEpoch,uzEpoch,uy);
  ucrss_c(uy,normalEpoch,ux);
  h2m[0][0]=ux[0];h2m[0][1]=ux[1];h2m[0][2]=ux[2];
  h2m[1][0]=uy[0];h2m[1][1]=uy[1];h2m[1][2]=uy[2];
  h2m[2][0]=normalEpoch[0];h2m[2][1]=normalEpoch[1];h2m[2][2]=normalEpoch[2];
  invert_c(h2m,h2i);
}

/*
  Matrix to convert from planetocentric coordinates to horizontal
  coordinates and viceversa.

  See discussion at:
  https://naif.jpl.nasa.gov/pipermail/spice_discussion/2010-July/000307.html

  h2m: converts from geocentric to topocentric
  h2i: converts from topocentric to geocentric
 */
void horgeo(SpiceDouble lat,SpiceDouble lon,SpiceDouble h2m[3][3],SpiceDouble h2i[3][3])
{
  SpiceDouble geopos[3],normal[3];
  SpiceDouble ux[]={1,0,0},uy[]={0,1,0},uz[]={0,0,1};

  //NORMAL ITRF93
  georec_c(D2R(lon),D2R(lat),0.0,REARTH,FEARTH,geopos);
  surfnm_c(REARTH,REARTH,RPEARTH,geopos,normal);

  //TRANSFORM MATRICES
  ucrss_c(normal,uz,uy);
  ucrss_c(uy,normal,ux);
  h2m[0][0]=ux[0];h2m[0][1]=ux[1];h2m[0][2]=ux[2];
  h2m[1][0]=uy[0];h2m[1][1]=uy[1];h2m[1][2]=uy[2];
  h2m[2][0]=normal[0];h2m[2][1]=normal[1];h2m[2][2]=normal[2];
  invert_c(h2m,h2i);
}

/*
  Transform from rectangular position to spherical position using the
  astronomical convention of Azimuth.
 */
int rec2hor(double pos[],double *Az,double *h)
{
  double phi,tmp;
  pos[1]*=-1;
  reclat_c(pos,&tmp,Az,h);
  if(*Az<0) *Az+=2*M_PI;
  *h=R2D(*h);
  *Az=R2D(*Az);
  return 0;
}

int copyVec(double tgt[],double src[],int n)
{
  memcpy(tgt,src,n*sizeof(double));
  return 0;
}

int sumVec(double c[],double ca,double a[],double cb,double b[],int n)
{
  int i;
  for(i=n;i-->0;) c[i]=ca*a[i]+cb*b[i];
  return 0;
}

double maxAbsVec(double a[],int n)
{
  int i;
  double max=-1E100;
  for(i=n;i-->0;) if(fabs(a[i])>max) max=fabs(a[i]);
  return max;
}

static int Rational_Extrapolation_to_Zero(double *fzero,double tableau[],
					  double x[],double f,int n) 
{
  double t, up, across, denominator, dum;
  int col;

  if (n==0) {  *fzero = f; tableau[0] = f; return 0; }
  if ( x[n] == 0.0 ) { *fzero = f; return -2; }
   
  across = 0.0;                                                        
  up = tableau[0];                                                    
  tableau[0] = f;                                               

  for (col = 1; col <= n; col++) {
    if(tableau[col-1]==0 && across==0){t=0;break;}
    denominator = tableau[col-1] - across;
    if (denominator == 0.0) return -1;
    dum = 1.0 - (tableau[col-1] - up) / denominator;
    denominator = (x[n - col] / x[n]) * dum - 1.0;
    if (denominator == 0.0) return -1;
    t = tableau[col-1] + ( tableau[col-1] - up ) / denominator;
    across = up;
    up = tableau[col];
    tableau[col] = t;
  }
  *fzero = t;
  return 0;
}

static int Polynomial_Extrapolation_to_Zero(double *fzero,double tableau[],
					    double x[], double f, int n )
{
  double back_two_columns;    //  T[row,col-2];
  double old_aux;             //  T[row-1,col];
  double new_value;           //  T[row,col];
  double vertical_diff;       //  T[row,col]-T[row-1,col]
  double backslant_diff;      //  T[row,col]-T[row,col-1]
  double forwardslant_diff;   //  T[row,col]-T[row-1,col-1];
  double denominator;        
  int i;

  if (n == 0) { tableau[0] = f; return 0; }
  if ( x[n] == 0.0 ) { *fzero = f; return -2; }

  back_two_columns = 0.0;
  old_aux = tableau[0];
  tableau[0] = f;
  for (i = 0; i < n; i++) {
    if(tableau[i]==0 && old_aux==0){tableau[n]=0.0;break;}
    vertical_diff = tableau[i] - old_aux;
    backslant_diff = tableau[i] - back_two_columns;
    forwardslant_diff = backslant_diff - vertical_diff;
    denominator = (x[n-i-1]/x[n]) * forwardslant_diff - backslant_diff;
    if (denominator == 0.0) return -1;
    back_two_columns = old_aux;
    old_aux = tableau[i+1];
    tableau[i+1] = tableau[i] + vertical_diff * backslant_diff / denominator;
  }
  *fzero = tableau[n];
  return 0;
}

static int Graggs_Method(int (*f)(double,double*,double*,void*),
			 double y0[],
			 double t0,double t,
			 int NUMBER_OF_STEPS,
			 void *params,
			 double yres[]) {
  
  double* pars=(double*)params;
  int order=(int)pars[0],i;
  double y1[order],dydt[order],y2[order],yaux[order];
  double h=(t-t0)/(double)NUMBER_OF_STEPS;
  double h2=h+h;

  copyVec(yaux,y0,order);
  (*f)(t0,yaux,dydt,params);
  sumVec(y1,1,yaux,h,dydt,order);

  while(--NUMBER_OF_STEPS) {
    t0+=h;
    (*f)(t0,y1,dydt,params);
    sumVec(y2,1,yaux,h2,dydt,order);
    copyVec(yaux,y1,order);
    copyVec(y1,y2,order);
  } 

  (*f)(t,y1,dydt,params);

  sumVec(yres,0.5,yaux,0.5,y1,order);
  sumVec(yres,1,yres,0.5*h,dydt,order);
  return 0;
}

int Gragg_Bulirsch_Stoer(int (*f)(double,double*,double*,void*), 
			 double y0[], double y1[],
			 double t, double h, double *h_new, 
			 double epsilon, double yscale, 
			 int rational_extrapolate,
			 void *params)
{
  double* pars=(double*)params;
  int order=(int)pars[0];
  double step_size2[ATTEMPTS];
  double tableau[order][ATTEMPTS+1];
  double dum;
  double est[order],dest[order],destmax;
  double old_est[order];
  
  int (*Extrapolate)(double*,double*,double*,double,int);
  int i,j;
  int err;

  if(yscale==0.0) return -3;
  if(rational_extrapolate) Extrapolate=Rational_Extrapolation_to_Zero;
  else Extrapolate=Polynomial_Extrapolation_to_Zero;
 
  Graggs_Method(f,y0,t,t+h,NUMBER_OF_STEPS[0],params,est);
  step_size2[0]=(dum=h/(double)NUMBER_OF_STEPS[0],dum*dum);
  
  copyVec(y1,est,order);
  
  for(i=order;i-->0;){
    err=Extrapolate(&y1[i],tableau[i],step_size2,est[i],0);
    if(err<0) return err-1;
  }

  for(i = 1; i < ATTEMPTS; i++) {
    copyVec(old_est,y1,order);
    Graggs_Method(f,y0,t,t+h,NUMBER_OF_STEPS[i],params,est);
    step_size2[i]=(dum=h/(double)NUMBER_OF_STEPS[i],dum*dum);

    for(j=order;j-->0;){
      err=Extrapolate(&y1[j],tableau[j],step_size2,est[j],i);
      if(err<0) return err-1;
    }
    
    sumVec(dest,1.0/yscale,y1,-1.0/yscale,old_est,order);
    destmax=maxAbsVec(dest,order);

    if(destmax<epsilon){
      if(i>1) *h_new=8.0*h/(double)NUMBER_OF_STEPS[i-1];
      else *h_new=h;
      return 0;
    }
  }
  return -1;
}

double energy2B(double X[])
{
  double v=vnorm_c(X+3);
  double r=vnorm_c(X);
  double E=0.5*v*v-1/r;
  return E;
}
