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
       ./whereisit MARS_BARYCENTER "07/19/2015 00:00:00.000 UTC-5"
       ./whereisit 4 "07/19/2015 00:00:00.000 UTC-5"

  */
  SpiceChar date[100],obj[100];
  strcpy(obj,argv[1]);
  strcpy(date,argv[2]);
  
  ////////////////////////////////////////////////////
  //GET EPHEMERIS TIME
  ////////////////////////////////////////////////////
  SpiceDouble t,tjd,ltmp,dt;

  //TIME IS INPUT IN UTC
  str2et_c(date,&t);

  //COMPUTE DELTAT FOR THE DATE
  deltet_c(t,"et",&dt);
  fprintf(stderr,"DT = %.2lf\n",dt);
  fprintf(stderr,"TT = %.9e\n",t);

  //CORRECT TIME FOR DELTAT.  NOW T IS TDB
  t-=dt;
  fprintf(stderr,"TDB = ");
  fprintf(stdout,"%.9e ",t);

  //JULIAN DATE TO VERIFY
  tjd=t2jd(t);
  fprintf(stderr,"\nJulian Date = %.6lf\n",tjd);

  ////////////////////////////////////////////////////
  //GET POSITION AT t
  ////////////////////////////////////////////////////
  SpiceDouble objectSSBJ2000[6];
  spkezr_c(obj,t,ABSJ2000,"NONE",SSB,objectSSBJ2000,&ltmp);
  fprintf(stderr,"State vector of %s: ",obj);
  fprintf(stdout,"%s",vec2strn(objectSSBJ2000,6,"%.17e "));
  fprintf(stderr,"\nLT = %.17e\n",ltmp);
}
