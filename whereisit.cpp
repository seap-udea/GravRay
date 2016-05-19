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

      Calculate the state vector (position and velocity) with respect
      to the Solar System Barycenter of a Solar System object.

    Arguments are: 

      date, object

    Date format:

       MM/DD/CCYY HH:MM:SS.dcm UTC-L

    Example:

       ./whereisit MARS "07/19/2015 00:00:00.000 UTC-5"
  */
  SpiceChar date[100],obj[100];
  strcpy(obj,argv[1]);
  strcpy(date,argv[2]);
  
  ////////////////////////////////////////////////////
  //GET EPHEMERIS TIME
  ////////////////////////////////////////////////////
  SpiceDouble t,tjd,ltmp;
  str2et_c(date,&t);
  printf("TT = %e\n",t);
  tjd=t2jd(t);
  printf("Julian Date = %.6lf\n",tjd);

  ////////////////////////////////////////////////////
  //GET POSITION AT t
  ////////////////////////////////////////////////////
  SpiceDouble objectSSBJ2000[6];
  spkezr_c(obj,t,"J2000","NONE","SOLAR SYSTEM BARYCENTER",objectSSBJ2000,&ltmp);
  printf("State vector of %s: %s\n",obj,vec2strn(objectSSBJ2000,6,"%.17e"));
}
