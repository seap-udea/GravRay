#include <gravray.cpp>
using namespace std;

#define TOLERANCE 1E-10
#define EXTMET 1

#include <objects.cpp>

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
    spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
    vscl_c(1E3/UL,object,object);
    sumVec(R,1.0,y,-1.0,object,3);
    Rmag=vnorm_c(R);
    GM=GMASSES[i]*fac;
    sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);
  }
  return 0;
}

int main(int argc,char* argv[])
{
  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();

  ////////////////////////////////////////////////////
  //GET ARGUMENTS
  ////////////////////////////////////////////////////
  /*
    Funcion: 

      Propagate the orbit of an object in time starting at a given
      position and time and including in the integration the combined
      gravitational effect of the major planets in the Solar System.

    Arguments are: 

      ephemeris time (in seconds), x, y, z, vx, vy, vz (in km and
      km/s), duration (in sideral years = 365.25), number of
      intermediate points (sampling points),

    Example:

      ./wherewillitbe.exe 4.80211267185610712e+08  -1.48536453607321531e+08 -8.62147140760204988e+05 -3.96733796100291074e+05 1.97063437949220379e+01 -2.74993804268180604e+01 -1.19226355280104350e+01 2.5 100

  */
  SpiceDouble tini=atof(argv[1]);
  SpiceDouble x=atof(argv[2]);
  SpiceDouble y=atof(argv[3]);
  SpiceDouble z=atof(argv[4]);
  SpiceDouble vx=atof(argv[5]);
  SpiceDouble vy=atof(argv[6]);
  SpiceDouble vz=atof(argv[7]);
  SpiceDouble duration=atof(argv[8])*365.25*GSL_CONST_MKSA_DAY;
  SpiceInt npoints=atoi(argv[9]);
  fprintf(stderr,"Integrating during %.17e seconds\n",duration);

  ////////////////////////////////////////////////////
  //PROCEED WITH THE INTEGRATION
  ////////////////////////////////////////////////////
  double params[]={6};

  //UNITS
  UL=GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  UM=MSUN;
  GGLOBAL=1.0;
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  UV=UL/UT;
  fprintf(stderr,"Units: UL = %.17e, UT = %.17e, UM = %.17e, UV = %.17e\n",UL,UT,UM,UV);

  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],E[8],a;
  vpack_c(x*1E3/UL,y*1E3/UL,z*1E3/UL,X0);
  vpack_c(vx*1E3/UV,vy*1E3/UV,vz*1E3/UV,X0+3);
  a=vnorm_c(X0);
  fprintf(stderr,"Initial conditions: %s\n",vec2strn(X0,6));

  //DYNAMICAL TIMESCALE
  double tdyn=2*M_PI*sqrt(a*a*a/(GGLOBAL*MSUN/UM));
  fprintf(stderr,"Dynamical time =  %e\n",tdyn);

  //TIME LIMITS
  duration/=UT;

  double h=tdyn/1000.0,h_used,h_next,h_adjust,deltat;
  double t_start=tini/UT;
  double t_step=duration/npoints;
  double tend=t_start+duration;
  double t_stop=tend;
  double t=t_start;

  fprintf(stderr,"tini = %lf\n",tini/UT);
  fprintf(stderr,"t_start = %lf\n",t_start);
  fprintf(stderr,"h = %e\n",h);
  fprintf(stderr,"t_step = %lf\n",t_step);
  fprintf(stderr,"t_stop = %lf\n",t_stop);
  fprintf(stderr,"tend = %lf\n",tend);
  //getc(stdin);

  //INTEGRATION
  int i,status;
  FILE *f=fopen("ray.dat","w");
  fprintf(f,"%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s\n",
	  "#1:t",
	  "2:x","3:y","4:z",
	  "5:vx","6:vy","7:vz",
	  "8:q","9:e","10:i",
	  "11:W","12:w",
	  "13:M","14:t0","15:mu");
  
  h_used=h;

  for(i=0;i<npoints;i++) {
    deltat=(t-tini/UT)*UT/YEAR;
    fprintf(stderr,"Step %d: t-t_start = %e yrs (last h = %e days)\n",i,deltat,h_used*UT/DAY);

    //CONVERTING TO CLASSICAL ELEMENTS
    vscl_c(UL/1E3,X0,Xu);vscl_c(UV/1E3,X0+3,Xu+3);
    oscelt_c(Xu,t*UT,GKMS*MSUN,E);
    E[0]*=1E3/AU;
    vscl_c(180/M_PI,E+2,E+2);

    //STORING RESULTS
    fprintf(f,"%-+26.17e%s%s\n",deltat,vec2strn(X0,6,"%-+26.17e"),vec2strn(E,8,"%-+26.17e"));

    if((t_start+t_step)>tend) t_step=tend-t_start;
    t_stop = t_start + t_step;

    fprintf(stderr,"\tt_start = %lf, t_step = %lf, t_stop = %lf\n",t_start,t_step,t_stop);
    h_used = h;
    do {
      while(1){
	status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_used,&h_next,1.0,TOLERANCE,EXTMET,params);
	if(status) h_used/=4.0;
	else break;
      }
      t+=h_used;
      copyVec(X0,X,6);
      if(t+h_next>t_stop) h_used=t+h_next-t_stop;
      else h_used=h_next;
    }while(t<t_stop-1.e-10);
    if(t>t_stop){
      h_adjust=t_stop-t;
      fprintf(stderr,"\tAdjusting fom t=%lf to t_stop=%lf using h = %e\n",t,t_stop,h_adjust);
      status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_adjust,&h_next,1.0,TOLERANCE,EXTMET,params);
      copyVec(X0,X,6);
      t=t_stop;
    }
    fprintf(stderr,"\tt = %lf, t_stop = %lf\n",t,t_stop);
    t_start = t;
    if(t_start>tend) break;
    //getc(stdin);
  }
  fclose(f);

  return 0;
}
