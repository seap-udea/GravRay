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

      object, (ET), date

    Date format:

       MM/DD/CCYY HH:MM:SS.dcm UTC-L

       If ET is used the date should be provided as ephemeris time.

    Example:

       ./whereisit.exe MARS "07/19/2015 00:00:00.000 UTC-5"
       ./whereisit.exe MARS_BARYCENTER "07/19/2015 00:00:00.000 UTC-5"
       ./whereisit.exe 4 "07/19/2015 00:00:00.000 UTC-5"
       ./whereisit.exe MARS_BARYCENTER ET 4.905360682e+08

  */
  SpiceChar date[100],obj[100];
  SpiceDouble t,tjd,ltmp,dt;

  strcpy(obj,argv[1]);
  strcpy(date,argv[2]);
  if(strcmp(date,"ET")==0){
    ////////////////////////////////////////////////////
    //GET EPHEMERIS TIME
    ////////////////////////////////////////////////////
    t=atof(argv[3]);
  }else{
    ////////////////////////////////////////////////////
    //GET DATE
    ////////////////////////////////////////////////////
    str2et_c(date,&t);
    deltet_c(t,"et",&dt);
    fprintf(stdout,"DT = %.2lf\n",dt);
    fprintf(stdout,"TT = %.9e\n",t);
    t-=dt;
  }
  //CORRECT TIME FOR DELTAT.  NOW T IS TDB
  fprintf(stdout,"TDB = %.9e\n",t);

  //JULIAN DATE TO VERIFY
  tjd=t2jd(t);
  fprintf(stdout,"Julian Date = %.6lf\n",tjd);

  ////////////////////////////////////////////////////
  //GET POSITION AT t
  ////////////////////////////////////////////////////
  SpiceDouble objectSSBJ2000[6];
  spkezr_c(obj,t,ECJ2000,"NONE",SSB,objectSSBJ2000,&ltmp);
  fprintf(stdout,"State vector of %s: %s\n",obj,vec2strn(objectSSBJ2000,6,"%.10e "));
  fprintf(stdout,"LT = %.17e\n",ltmp);
  
  ////////////////////////////////////////////////////
  //PLAIN OUTPUT
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
  fprintf(stderr,"TDB,JD,DT,STATE(6),LT\n");
  fprintf(stderr,"%.9e\n%.6lf\n%.2lf\n",t,tjd,dt);
  fprintf(stderr,"%s\n",vec2strn(objectSSBJ2000,6,"%+.17e "));
  fprintf(stderr,"%.6e\n",ltmp);
}
