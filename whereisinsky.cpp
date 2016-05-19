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

      Calculate the position in sky of a Solar System Object w.r.t. an
      observer in the surface of the Earth at latitude (degrees),
      longitud (degrees) and altitude (meters) and at a given date.

    Arguments are: 

      object, latitude (degrees), longitude (degrees), elevation
      (meters), date

    Date format:

       MM/DD/CCYY HH:MM:SS.dcm UTC-L

    Example:

       ./whereisinsky.exe MARS 6.2 -75.34 1450.0 "07/19/2015 00:00:00.000 UTC-5"
       ./whereisinsky.exe MARS_BARYCENTER 6.2 -75.34 1450.0 "07/19/2015 00:00:00.000 UTC-5"
       ./whereisinsky.exe 4 6.2 -75.34 1450.0 "07/19/2015 00:00:00.000 UTC-5"
  */
  SpiceChar date[100],obj[100];
  strcpy(obj,argv[1]);
  SpiceDouble lat=atof(argv[2]);
  SpiceDouble lon=atof(argv[3]);
  SpiceDouble alt=atof(argv[4]);
  strcpy(date,argv[5]);
  
  ////////////////////////////////////////////////////
  //GET EPHEMERIS TIME
  ////////////////////////////////////////////////////
  SpiceDouble t,tjd,ltmp;
  str2et_c(date,&t);
  tjd=t2jd(t);
  printf("Julian Date = %.6lf\n",tjd);

  ////////////////////////////////////////////////////
  //GET OBJECT ABSOLUTE EPHEMERIS
  ////////////////////////////////////////////////////
  double distance,ltime,raJ2000,decJ2000,ra,dec;
  bodyEphemerisApparent(obj,t,lon,lat,alt,
			&distance,&ltime,
			&raJ2000,&decJ2000,
			&ra,&dec
			);

  ////////////////////////////////////////////////////
  //GET OBJECT LOCAL EPHEMERIS
  ////////////////////////////////////////////////////
  SpiceDouble objObserverEquatorialEpoch[3];
  SpiceDouble objObserverLocalEpoch[3],Az,h;
  SpiceDouble hm[3][3],hi[3][3];
  hormat(lat,lon,t,hm,hi);
  latrec_c(1.0,D2R(ra*15.0),D2R(dec),objObserverEquatorialEpoch);
  mxv_c(hm,objObserverEquatorialEpoch,objObserverLocalEpoch);
  rec2hor(objObserverLocalEpoch,&Az,&h);
  printf("\tEPOCH: Az = %lf, h = %lf\n",Az,h);

  ////////////////////////////////////////////////////
  //SHOW OBJECT POSITION
  ////////////////////////////////////////////////////
  printf("Position of %s in sky w.r.t. Observer:\n",obj);
  printf("\tDistance = %.17e km\n",distance);
  printf("\tLight time = %.17e seconds\n",ltime);
  printf("\tJ2000: RA = %s, DEC = %s\n",dec2sex(raJ2000),dec2sex(decJ2000));
  printf("\tEPOCH: RA = %s, DEC = %s\n",dec2sex(ra),dec2sex(dec));
  printf("\tEPOCH: Az = %s, h = %s\n",dec2sex(Az),dec2sex(h));
  printf("\tEPOCH: Az = %lf, h = %lf\n",Az,h);

  return 0;
}
