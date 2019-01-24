#include <gravray.cpp>
using namespace std;

int main(int argc,char* argv[])
{
  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();

  ////////////////////////////////////////////////////
  //INPUTS
  ////////////////////////////////////////////////////
  SpiceChar file[1000],outfile[1000];
  SpiceDouble t,lat,lon,alt;
  int i,qvel;
  SpiceDouble vproj,qapex,dapex,earth[6],uvearth[3],ltmp;
  SpiceDouble vloc[3],vmot[3],uvimpact[3],uvimpactEJ2000[3];
  SpiceDouble cA,sA,ch,sh;
  int nvel,iapex;
  double vels[100];

  if(argc==8){
    t=atof(argv[1]);
    lat=atof(argv[2]);
    lon=atof(argv[3]);
    alt=atof(argv[4]);
    strcpy(file,argv[5]);
    qvel=atoi(argv[6]);
    strcpy(outfile,argv[7]);
  }else
    argsError(argv[0]);
  
  ////////////////////////////////////////////////////
  //INITIALIZE OBSERVER
  ////////////////////////////////////////////////////
  struct ObserverStruct observer;
  observer.lat=lat;
  observer.lon=lon;
  observer.alt=alt;
  initObserver(t,&observer);
  fprintf(stdout,"%s\n",vec2str(observer.v));
  exit(0);
  
  ////////////////////////////////////////////////////
  //READ INITIAL CONDITIONS
  ////////////////////////////////////////////////////
  int ncond=1;
  int ncoll=0;
  SpiceDouble h,Az,v,speed;
  SpiceDouble elements[6];
  SpiceDouble deltat=-0.5;
  char line[1000];

  FILE* fi=fopen(file,"r");
  fgets(line,1000,fi);

  FILE* fe=fopen(outfile,"w");
  fprintf(fe,"%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s\n",
	  "#1:h","2:Az","3:vimp",
	  "4:xobs","5:yobs","6:zobs","7:vxobs","8:vyobs","9:vzobs",
	  "10:q(AU)","11:e","12:i","13:Omega","14:omega","15:M","16:qapex");

  ////////////////////////////////////////////////////
  //VELOCITY OF THE EARTH
  ////////////////////////////////////////////////////
  spkezr_c("EARTH",t,ECJ2000,"NONE",SSB,earth,&ltmp);
  unorm_c(earth+3,uvearth,&ltmp);

  ////////////////////////////////////////////////////
  //SIMULATION
  ////////////////////////////////////////////////////
  while(1){

    fscanf(fi,"%lf %lf",&h,&Az);
    if(feof(fi)) break;

    v=0;i=0;
    while(1){
      fscanf(fi,"%lf",&v);
      if(v==-1) break;
      vels[i++]=v;
    }
    nvel=i;
    if(feof(fi)) break;

    //==================================================
    //DETERMINE DIRECTION WITH RESPECT TO APEX
    //==================================================
    cA=cos(D2R(Az));sA=sin(D2R(Az));
    ch=cos(D2R(h));sh=sin(D2R(h));

    //DIRECTION OF IMPACT RESPECT TO OBSERVER REFERENCE FRAME
    if(Az==0 && h==0){vpack_c(1,0,0,vloc);}
    else{vpack_c(ch*cA,-ch*sA,sh,vloc);}

    //DIRECTION OF IMPACT W.R.T. TO ITR93
    mxv_c(observer.hi,vloc,uvimpact);
      
    //DIRECTION OF IMPACT W.R.T. TO ECLIPJ2000
    mxv_c(observer.MEJ,uvimpact,uvimpactEJ2000);

    //COSINE QAPEX
    vproj=vdot_c(uvimpactEJ2000,uvearth);

    //QAPEX
    qapex=R2D(acos(vproj));

    if(nvel>1 && qvel){
      dapex=180.0/(nvel-1);
      iapex=floor(qapex/dapex)+1;
      v=vels[iapex];
    }else{
      v=vels[0];
    }

    fprintf(stdout,"Condition %d at lat=%.5lf,lon=%.5lf: Integrating for h = %.2lf, A = %.2lf, v = %.2lf (qvel? %d), qapex = %.2f\n",
	    ncond,lat,lon,h,Az,v,qvel,qapex);

    ////////////////////////////////////////////////////
    //DETERMINE POSITION FOR THIS INITIAL CONDITION
    ////////////////////////////////////////////////////
    observerVelocity(&observer,h,Az,v);
    fprintf(stdout,"\tObserver position (ECLIPJ2000): %s\n",vec2strn(observer.posabs,6,"%.8e "));
    speed=vnorm_c(observer.posabs+3);
    fprintf(stdout,"\tObserver speed: %e\n",speed);

    ////////////////////////////////////////////////////
    //PROPAGATE
    ////////////////////////////////////////////////////
    try {
      rayPropagation(&observer,deltat,elements);
    } 
    catch (int e) {
      fprintf(stdout,"\t\tA collision occurred. Skipping point\n");
      ncoll++;
      continue;
    }

    ////////////////////////////////////////////////////
    //SAVE
    ////////////////////////////////////////////////////
    fprintf(stdout,"\tFinal elements for ray %d: %s\n",
	    ncond,
	    vec2strn(elements,6,"%lf "));
    fprintf(fe,"%-+26.17e%-+26.17e%-+26.17e%s%s%-+26.17e\n",
	    h,Az,v,
	    vec2strn(observer.posabs,6,"%-26.17e"),
	    vec2strn(elements,6,"%-+26.17e"),
	    qapex);
    ncond++;
  }
  fclose(fi);
  fclose(fe);

  fprintf(stdout,"\nNumber of initial conditions: %d\n",ncond-1);
  fprintf(stdout,"Number of collisions: %d\n",ncoll);

  return 0;
}
