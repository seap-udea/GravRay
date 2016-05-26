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

      ./launchmany.exe 4.141704340e+08 54.456093 63.492323 8.234765e+04 initial.dat

  */
  ////////////////////////////////////////////////////
  //INPUTS
  ////////////////////////////////////////////////////
  SpiceChar file[1000];
  SpiceDouble t=atof(argv[1]);
  SpiceDouble lat=atof(argv[2]);
  SpiceDouble lon=atof(argv[3]);
  SpiceDouble alt=atof(argv[4]);
  strcpy(file,argv[5]);

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
  int ncond=0;
  SpiceDouble h,Az,v;
  SpiceDouble elements[6];
  SpiceDouble deltat=-2.0;
  char line[1000];

  FILE* fi=fopen(file,"r");
  fgets(line,1000,fi);
  while(1){
    fscanf(fi,"%lf %lf %lf",&h,&Az,&v);
    if(feof(fi)) break;
    
    fprintf(stdout,"Integrating for h = %.2lf, A = %.2lf, v = %.2lf\n",
	    h,Az,v);

    ////////////////////////////////////////////////////
    //DETERMINE POSITION FOR THIS INITIAL CONDITION
    ////////////////////////////////////////////////////
    observerVelocity(&observer,h,Az,v);
    fprintf(stdout,"Observer position: %s\n",vec2strn(observer.posabs,6,"%e "));

    rayPropagation(&observer,deltat,elements);
    fprintf(stdout,"Final elements: %s\n",vec2strn(elements,6,"%lf "));

    ncond++;
  }
  fclose(fi);

  return 0;
}
