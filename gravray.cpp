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
#define MAX(x,y) (x>y?x:y)

//////////////////////////////////////////
//CSPICE CONSTANTS
//////////////////////////////////////////
#define EARTH_ID "EARTH"
#define ATTEMPTS 12 /*SEE NUMBER_OF_STEPS*/
#define SSB "SOLAR SYSTEM BARYCENTER"

//FOR ABSOLUTE EPHEMERIS
#define ECJ2000 "ECLIPJ2000"
//FOR LOCAL EPHEMERIS
#define J2000 "J2000"

/*
  Range of Ephemeris Times where data to calculate precise ITRF93
  frame directions is availanle
 */

// 01/01/2000 00:00:00.000 UTC
//#define ETINI -4.313581609e+04

// 01/21/1962 00:00:00.000 UTC
#define ETINI -1.197460759e+09
// 07/17/2037 00:00:00.000 UTC
#define ETEND 1.184673668e+09

//////////////////////////////////////////
//CONSTANTS
//////////////////////////////////////////
#define GCONST GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
#define GKMS (GCONST*1E-9)
#define MSUN 1.9885E30/*kg*/
#define YEAR (365.25*GSL_CONST_MKSA_DAY)
#define DAY GSL_CONST_MKSA_DAY
#define AU GSL_CONST_MKSA_ASTRONOMICAL_UNIT
#define VESC_EARTH 11.217 //km/s
//////////////////////////////////////////
//BEHAVIOR
//////////////////////////////////////////
#define TOLERANCE 1E-10
#define EXTMET 1
#define VERBOSE 0
#define HTOL 1E-6
#define MAXSTALL 100

//////////////////////////////////////////
//OBJECTS
//////////////////////////////////////////
#define NUMOBJS 10

//THESE ARE THE LABELS FOR KERNEL de430
static char* LABELS[]={
  "SUN",
  "MERCURY",
  "VENUS",
  "EARTH",
  "MOON",
  "MARS",
  "JUPITER",
  "SATURN",
  "URANUS",
  "NEPTUNE"
};


static char* OBJS[]={
  "10",/*SUN*/
  "1",/*MERCURY*/
  "2",/*VENUS*/
  "399",/*EARTH*/
  "301",/*MOON*/
  "4",/*MARS*/
  "5",/*JUPITER*/
  "6",/*SATURN*/
  "7",/*URANUS*/
  "8"/*NEPTUNE*/
};

//SEE WIKIPEDIA
static double MASSES[]={
  1.9891E30/*SUN*/,
  3.3022E23/*MERCURY*/,
  4.8685E24/*VENUS*/,
  5.9736E24/*EARTH*/,
  7.349E22/*MOON*/,
  6.4185E23/*MARS*/,
  1.8986E27/*JUPITER*/,
  5.6846E26/*SATURN*/,
  8.6810E25/*URANUS*/,
  1.0243E26/*NEPTUNE*/
};

/*
  Source: documentation DE421
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421_announcement.pdf
 */
static double GMASSES[]={
  132712440040.944000/*SUN*/,
  22032.090000/*MERCURY*/,
  324858.592000/*VENUS*/,
  398600.436233/*EARTH*/,
  4902.800076/*MOON*/,
  42828.375214/*MARS*/,
  126712764.800000/*JUPITER*/,
  37940585.200000/*SATURN*/,
  5794548.600000/*URANUS*/,
  6836535.000000/*NEPTUNE*/
};

/*
  Gotten from kernes DE430
 */
static double RADII[]={
  696000/*SUN*/,
  2439.70/*MERCURY*/,
  6051.80/*VENUS*/,
  6378.14/*EARTH*/,
  1737.40/*MOON*/,
  3396.19/*MARS*/,
  71492.00/*JUPITER*/,
  37940585.200000/*SATURN*/,
  5794548.600000/*URANUS*/,
  6836535.000000/*NEPTUNE*/
};


struct ObserverStruct{

  SpiceDouble t;
  SpiceDouble lat,lon,alt;

  //Conversion matrix from ITRF93 to ECLIPJ2000 at time t
  SpiceDouble MEJ[3][3];

  //Conversion matrix from ECLIPJ2000 to EPOCHEQUINOX at time t
  SpiceDouble MEE[3][3];

  //hm convert from itrf93 to local and hi is the inverse
  SpiceDouble hm[3][3];
  SpiceDouble hi[3][3];

  //Position with respect to ITRF93 (Earth) J2000
  SpiceDouble posearth[6];

  //Position with respect to ITRF93 (Earth) ECLIPJ2000
  SpiceDouble posj2000[6];

  //Position with respect to SSB in ECLIPJ2000
  SpiceDouble posabs[6];

  //Position with respect to SSB in ECLIPEPOCH
  SpiceDouble posepoch[6];

  //Rotaion velocity of a still observer with respect to ITRF93
  SpiceDouble v[3];

  //Direction of velocity in the direction (A,a) with respect to ECLIPJ2000
  SpiceDouble uv[3];

  //Earth position at observer epoch
  SpiceDouble earth[6];
};

//////////////////////////////////////////
//GLOBAL VARIABLES
//////////////////////////////////////////
enum COMPONENTS {CX,CY,CZ,CVX,CVY,CVZ};
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
  SpiceInt i,n;
  SpiceDouble radii[3];

  //KERNELS
  furnsh_c("kernels.txt");

  //EARTH RADII
  bodvrd_c(EARTH_ID,"RADII",3,&n,radii);
  REARTH=radii[0];
  RPEARTH=radii[0];
  FEARTH=(radii[0]-radii[2])/radii[0];

  //PLANETARY RADII
  for(i=0;i<10;i++){
    bodvrd_c(LABELS[i],"RADII",3,&n,radii);
    RADII[i]=radii[0];
  }

  //RANDOM NUMBERS
  RAND=gsl_rng_alloc(gsl_rng_default);

  return 0;
}

int initRadii(void)
{

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

char* vec2str(double vec[],char frm[]="%.8e ")
{
  char format[100];
  char *str=(char*)calloc(sizeof(char),100); 
  sprintf(format,"%s %s %s",frm,frm,frm);
  sprintf(str,format,vec[0],vec[1],vec[2]);
  return str;
}

char* vec2strn(double vec[],int n,char frm[]="%.8e ")
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
  spkezr_c(EARTH_ID,t,"J2000","NONE",SSB,
	   earthSSBJ2000,&ltmp);
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);
  while((fabs(lt-ltold)/lt)>=lttol && i<ncn){
    ltold=lt;
    spkezr_c(body,t-lt,"J2000","NONE",SSB,bodySSBJ2000,&ltmp);
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

/*
Adapted from: http://www.mymathlib.com/diffeq/bulirsch_stoer.html
 */

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

#include <objects.hpp>

int EoM(double t,double y[],double dydt[],void *params) 
{ 
  //COMPUTE THE CONTRIBUTION OF EVERY OBJECT
  int i;
  double r,object[6],R[3],Rmag,tmp,GM,fac;

  fac=UT*UT/(UL/1E3*UL/1E3*UL/1E3);
  dydt[CX]=y[CVX];
  dydt[CY]=y[CVY];
  dydt[CZ]=y[CVZ];
  dydt[CVX]=0.0;
  dydt[CVY]=0.0;
  dydt[CVZ]=0.0;

  for(i=NUMOBJS;i-->0;){
    if(!ACTIVE[i]) continue;
    spkezr_c(OBJS[i],t*UT,ECJ2000,"NONE",SSB,object,&tmp);
    vscl_c(1E3/UL,object,object);
    sumVec(R,1.0,y,-1.0,object,3);
    Rmag=vnorm_c(R);
    if(Rmag*UL/1e3<=RADII[i]){
      fprintf(stderr,"\t\tObject has collided with %s at t = %e days (Rmag = %e, RADII = %e)\n",
	     LABELS[i],t*UT/DAY,Rmag,RADII[i]);
      throw(1);
    }else{
      GM=GMASSES[i]*fac;
      /*
      if(i==3){
	fprintf(stdout,"t=%e, i = %d, R = %e\n",t,i,Rmag);
	getc(stdin);
      }
      */
      sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);
    }
  }
  return 0;
}

int initObserver(SpiceDouble t,struct ObserverStruct* observer)
{
  SpiceDouble rho,vcirc,vrot[3];
  SpiceDouble lt;
  
  observer->t=t;
  
  //CHECK DATES
  SpiceDouble tref;
  if(t>=ETINI && t<=ETEND){
    tref=t;
  }else{
    SpiceChar UTC[100];
    SpiceDouble dt;
    deltet_c(t,"et",&dt);
    et2utc_c(t+dt,"ISOC",2,100,UTC);
    if(t<ETINI){
      UTC[0]='1';UTC[1]='9';UTC[2]='6';UTC[3]='3';
      str2et_c(UTC,&tref);
      fprintf(stdout,"ETINI = %.10e, t = %.10e, tref = %.10e\n",ETINI,t,tref);
    }
    else{
      UTC[0]='2';UTC[1]='0';UTC[2]='3';UTC[3]='2';
      str2et_c(UTC,&tref);
      fprintf(stdout,"ETEND = %.10e, t = %.10e, tref = %.10e\n",ETEND,t,tref);
    }
  }
  //DEBUGGING
  //printf("t = %e, tref = %e\n",t,tref);

  //CONVERSION FROM EARTH SYSTEM TO ECLIPTIC SYSTEM AT TIME T
  pxform_c("ITRF93",ECJ2000,tref,observer->MEJ);
  /*
  printf("%e,%e,%e\n%e,%e,%e\n%e,%e,%e\n",
	 observer->MEJ[0][1],observer->MEJ[0][2],observer->MEJ[0][3],
	 observer->MEJ[1][1],observer->MEJ[1][2],observer->MEJ[1][3],
	 observer->MEJ[2][1],observer->MEJ[2][2],observer->MEJ[2][3]);
  exit(0);
  */

  //CONVERSION FROM EARTH SYSTEM TO ECLIPTIC SYSTEM AT TIME T
  pxform_c("ECLIPJ2000","EARTHTRUEEPOCH",tref,observer->MEE);

  //LOCATE OBSERVER 
  georec_c(D2R(observer->lon),D2R(observer->lat),observer->alt/1000.0,
	   REARTH,FEARTH,observer->posearth);

  //DEBUGGING
  //printf("Observer = %s\n",vec2str(observer->posearth));

  //TOPOCENTRIC CONVERSION MATRICES
  horgeo(observer->lat,observer->lon,observer->hm,observer->hi);

  //VELOCITY OF OBSERVER DUE TO EARTH ROTATION
  rho=sqrt(observer->posearth[0]*observer->posearth[0]+
	   observer->posearth[1]*observer->posearth[1]);
  vcirc=2*M_PI*rho/GSL_CONST_MKSA_DAY;
  vpack_c(0.0,-vcirc,0.0,vrot);
  mxv_c(observer->hi,vrot,observer->v);

  //POSITION OF THE EARTH
  spkezr_c(EARTH_ID,t,ECJ2000,"NONE","SOLAR SYSTEM BARYCENTER",
	   observer->earth,&lt);

  //DEBUGGING
  //printf("Earth ECJ2000 = %s\n",vec2str(observer->earth,"%.17e "));

  //POSITION WITH RESPECT TO SSB IN ECLIPJ2000
  mxv_c(observer->MEJ,observer->posearth,observer->posj2000);
  vadd_c(observer->earth,observer->posj2000,observer->posabs);

  //DEBUGGING
  //printf("Observer ECJ2000 = %s\n",vec2str(observer->posabs,"%.17e "));

  //POSITION WITH RESPECT TO SSB IN ECLIPEPOCH
  //This is not working
  mxv_c(observer->MEE,observer->posabs,observer->posepoch);
  mxv_c(observer->MEE,(observer->posabs)+3,(observer->posepoch)+3);

}

int observerVelocity(struct ObserverStruct *observer,
		     SpiceDouble elev,SpiceDouble Az,SpiceDouble v)
{
  ////////////////////////////////////////////////////////////// 
  //DETERMINE THE OBSERVER VELOCITY WITH RESPECT TO SSB
  ////////////////////////////////////////////////////////////// 

  //VELOCITY OF OBSERVER IN SPACE W.R.T. TO LOCAL REFERENCE
  SpiceDouble cA=cos(D2R(Az)),sA=sin(D2R(Az)),ch=cos(D2R(elev)),sh=sin(D2R(elev));
  SpiceDouble vloc[3];
  if(Az==0 && elev==0 && v==0){
    vpack_c(0,0,0,vloc);
  }else{
    vpack_c(v*ch*cA,-v*ch*sA,v*sh,vloc);
  }

  //IMPACT VELOCITY IS THE INVERSE
  vscl_c(-1,vloc,vloc);

  //VELOCITY OF OBSERVER IN SPACE W.R.T. TO ITRF93
  SpiceDouble vmot[3];
  mxv_c(observer->hi,vloc,vmot);

  //TOTAL VELOCITY WITH RESPECT ITRF93
  vadd_c(observer->v,vmot,observer->posearth+3);

  //VELOCITY W.R.T. EARTH CENTER IN ECLIPJ2000 RF
  mxv_c(observer->MEJ,observer->posearth+3,observer->posj2000+3);

  //VELOCITY W.R.T. SOLAR SYSTEM BARYCENTER IN J2000 RF
  vadd_c(observer->earth+3,observer->posj2000+3,observer->posabs+3);

  /*NEW*/
  //************************************************************
  //COMPUTING DIRECTION OF INCOMING VELOCITY IN ECLIPJ2000
  /*
    It does not take into account Earth rotation effect on velocity
   */
  //************************************************************
  SpiceDouble uv[3],nuv;
  if(Az==0 && elev==0 && v==0){
    vpack_c(0,0,0,vloc);
  }else{
    vpack_c(ch*cA,-ch*sA,sh,uv);
  }
  //IMPACT DIRECTION IS THE INVERSE
  vscl_c(-1,uv,uv);
  //DIRECTION IN SPACE W.R.T. TO ITRF93
  mxv_c(observer->hi,uv,vmot);
  //DIRECTION IN SPACE W.R.T. ECLIPJ2000
  mxv_c(observer->MEJ,vmot,observer->uv);

  return 0;
}

int rayPropagation(struct ObserverStruct *observer,
		   SpiceDouble deltat,
		   SpiceDouble elements[6])
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIAL CONDITIONS FOR PROPAGATION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SpiceDouble x=observer->posabs[0];
  SpiceDouble y=observer->posabs[1];
  SpiceDouble z=observer->posabs[2];
  SpiceDouble vx=observer->posabs[3];
  SpiceDouble vy=observer->posabs[4];
  SpiceDouble vz=observer->posabs[5];

  deltat*=365.25*GSL_CONST_MKSA_DAY;
  SpiceDouble tini=observer->t;
  double direction=deltat/abs(deltat);
  double params[]={6};

  //UNITS
  UL=GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  UM=MSUN;
  GGLOBAL=1.0;
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  UV=UL/UT;

  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],E[8],a;
  vpack_c(x*1E3/UL,y*1E3/UL,z*1E3/UL,X0);
  vpack_c(vx*1E3/UV,vy*1E3/UV,vz*1E3/UV,X0+3);

  //DYNAMICAL TIMESCALE
  a=vnorm_c(X0);
  double tdyn=2*M_PI*sqrt(a*a*a/(GGLOBAL*MSUN/UM));

  //TIME LIMITS
  deltat/=UT;

  //GUESS TIME-STEP AS 1/1000 OF THE CHARACTERISTIC DYNAMICAL TIME
  double h=direction*tdyn/1000.0,h_used,h_next,h_adjust,delt;

  double t_start=tini/UT;
  double t_step=deltat;
  double tend=t_start+deltat;
  double t_stop=tend;

  double t=t_start;

  //INTEGRATION
  int status;

  t_stop = t_start + t_step;

  h_used = h;
  int nstall=0;
  do {
    //ADJUST H UNTIL OBTAINING A PROPER TIMESTEP
    while(1){
      status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_used,&h_next,1.0,TOLERANCE,EXTMET,params);
      if(status) h_used/=4.0;
      else break;
    }
    if(fabs(h_used/t_step)<HTOL) nstall++;
    else nstall=0;
    if(nstall>MAXSTALL){
      fprintf(stderr,"\t\tIntegration has stalled at t = %e days with h/DT = %e\n",
	      t*UT/DAY,h_used/t_step);
      throw(1);
    }
    t+=h_used;
    copyVec(X0,X,6);
    if(direction*(t+h_next-t_stop)>0) h_used=t+h_next-t_stop;
    else h_used=h_next;
  }while(direction*(t-(t_stop-direction*1.e-10))<0);

  //PREVIOUS PROCEDURE WILL LEAVE YOU STILL APART FROM FINAL TIME, SO ADJUST
  if(direction*(t-t_stop)>0){
    h_adjust=(t_stop-t);
    status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_adjust,&h_next,1.0,TOLERANCE,EXTMET,params);
    copyVec(X0,X,6);
    t=t_stop;
  }

  //CONVERTING TO CLASSICAL ELEMENTS IN KM AND KM/S
  vscl_c(UL/1E3,X0,Xu);vscl_c(UV/1E3,X0+3,Xu+3);
  oscelt_c(Xu,t*UT,GKMS*MSUN,E);

  vsclg_c(180/M_PI,E+2,4,E+2);
  vsclg_c(1E3/UL,E,1,E);

  //STORE THE ELEMENTS
  copyVec(elements,E,6);

  return 0;
}

int argsError(char* pname,char* msg="Bad options.")
{
    char cmd[1000];
    sprintf(cmd,"cat .help/%s.help",pname);
    fprintf(stderr,msg);
    system(cmd);
    exit(1);
}

char *str_replace(char *orig, char *rep, char *with) 
{
  /*
    Source:
    http://stackoverflow.com/questions/779875/what-is-the-function-to-replace-string-in-c
    You must free the result if result is non-NULL.
  */
  char *result; // the return string
  char *ins;    // the next insert point
  char *tmp;    // varies
  int len_rep;  // length of rep
  int len_with; // length of with
  int len_front; // distance between rep and end of last rep
  int count;    // number of replacements

  if (!orig)
    return NULL;
  if (!rep)
    rep = "";
  len_rep = strlen(rep);
  if (!with)
    with = "";
  len_with = strlen(with);

  ins = orig;
  for (count = 0; tmp = strstr(ins, rep); ++count) {
    ins = tmp + len_rep;
  }

  // first time through the loop, all the variable are set correctly
  // from here on,
  //    tmp points to the end of the result string
  //    ins points to the next occurrence of rep in orig
  //    orig points to the remainder of orig after "end of rep"
  tmp = result = (char*) malloc(strlen(orig) + (len_with - len_rep) * count + 1);

  if (!result)
    return NULL;

  while (count--) {
    ins = strstr(orig, rep);
    len_front = ins - orig;
    tmp = strncpy(tmp, orig, len_front) + len_front;
    tmp = strcpy(tmp, with) + len_with;
    orig += len_front + len_rep; // move to next "end of rep"
  }
  strcpy(tmp, orig);
  return result;
}
