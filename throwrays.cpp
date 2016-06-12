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

  if(argc==7){
    t=atof(argv[1]);
    lat=atof(argv[2]);
    lon=atof(argv[3]);
    alt=atof(argv[4]);
    strcpy(file,argv[5]);
    strcpy(outfile,argv[6]);
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

  ////////////////////////////////////////////////////
  //READ INITIAL CONDITIONS
  ////////////////////////////////////////////////////
  int ncond=1;
  int ncoll=0;
  SpiceDouble h,Az,v,speed;
  SpiceDouble elements[6];
  SpiceDouble deltat=-2.0;
  char line[1000];

  FILE* fi=fopen(file,"r");
  fgets(line,1000,fi);

  FILE* fe=fopen(outfile,"w");
  fprintf(fe,"%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s\n",
	  "#1:h","2:Az","3:vimp",
	  "4:xobs","5:yobs","6:zobs","7:vxobs","8:vyobs","9:vzobs",
	  "10:q(AU)","11:e","12:i","13:Omega","14:omega","15:M");

  while(1){
    fscanf(fi,"%lf %lf %lf",&h,&Az,&v);
    if(feof(fi)) break;
    
    fprintf(stdout,"Condition %d:Integrating for h = %.2lf, A = %.2lf, v = %.2lf\n",
	    ncond,h,Az,v);

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
    fprintf(fe,"%-+26.17e%-+26.17e%-+26.17e%s%s\n",
	    h,Az,v,
	    vec2strn(observer.posabs,6,"%-26.17e"),
	    vec2strn(elements,6,"%-+26.17e"));
    ncond++;
  }
  fclose(fi);
  fclose(fe);

  fprintf(stdout,"\nNumber of initial conditions: %d\n",ncond-1);
  fprintf(stdout,"Number of collisions: %d\n",ncoll);

  return 0;
}
