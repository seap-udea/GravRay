#include <gravray.cpp>
using namespace std;

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

      Launch many objects from a given location on Earth.

    Arguments are: 

      ephemeris time (in seconds), latitude (degrees), longitude
      (degrees), elevation (meters), file initial conditions

    Example:

      ./launchmany.exe 4.141704340e+08 54.456093 63.492323 8.234765e+04 initial.dat.temp elements.dat

      These are the initial conditions for the Chelyabinsk
      meteoroid. If everything is correct the result will be:

      Integrating for h = 17.66, A = 104.98, v = 20.46 Observer
      position: -1.232674e+08 8.134791e+07 -3.465527e+03 -3.359174e+01
      -1.543290e+01 -6.384145e+00 

      Final elements: 0.702040 0.612004 5.848491 326.515954 102.986331
      81.443436

  */
  ////////////////////////////////////////////////////
  //INPUTS
  ////////////////////////////////////////////////////
  SpiceChar file[1000],outfile[1000];
  SpiceDouble t=atof(argv[1]);
  SpiceDouble lat=atof(argv[2]);
  SpiceDouble lon=atof(argv[3]);
  SpiceDouble alt=atof(argv[4]);
  strcpy(file,argv[5]);
  strcpy(outfile,argv[6]);

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
  SpiceDouble h,Az,v;
  SpiceDouble elements[6];
  SpiceDouble deltat=-2.0;
  char line[1000];

  FILE* fi=fopen(file,"r");
  fgets(line,1000,fi);

  FILE* fe=fopen(outfile,"w");
  fprintf(fe,"%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s%-26s\n",
	  "#1:xobs","2:yobs","3:zobs","4:vxobs","5:vyobs","6:vzobs",
	  "7:q(AU)","8:e","9:i","10:Omega","11:omega","12:M");

  while(1){
    fscanf(fi,"%lf %lf %lf",&h,&Az,&v);
    if(feof(fi)) break;
    
    fprintf(stdout,"Condition %d:Integrating for h = %.2lf, A = %.2lf, v = %.2lf\n",
	    ncond,h,Az,v);

    ////////////////////////////////////////////////////
    //DETERMINE POSITION FOR THIS INITIAL CONDITION
    ////////////////////////////////////////////////////
    observerVelocity(&observer,h,Az,v);
    fprintf(stdout,"\tObserver position: %s\n",vec2strn(observer.posabs,6,"%.17e "));
    ncond++;

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
    //Components of elements vector: q,e,i,Omega,omega,M
    fprintf(stdout,"\tFinal elements: %s\n",vec2strn(elements,6,"%lf "));
    fprintf(fe,"%s%s\n",
	    vec2strn(observer.posabs,6,"%-26.17e"),
	    vec2strn(elements,6,"%-+26.17e"));
  }
  fclose(fi);
  fclose(fe);

  fprintf(stdout,"\nNumber of initial conditions: %d\n",ncond);
  fprintf(stdout,"Number of collisions: %d\n",ncoll);

  return 0;
}
