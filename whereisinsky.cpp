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
  SpiceChar date[100],obj[100];
  SpiceDouble lat,lon,alt;
  
  if(argc>=6){
    strcpy(obj,argv[1]);
    lat=atof(argv[2]);
    lon=atof(argv[3]);
    alt=atof(argv[4]);
    strcpy(date,argv[5]);
  }else
    argsError(argv[0]);
  
  ////////////////////////////////////////////////////
  //GET EPHEMERIS TIME
  ////////////////////////////////////////////////////
  SpiceDouble t,tjd,ltmp;
  str2et_c(date,&t);
  tjd=t2jd(t);
  fprintf(stdout,"Julian Date = %.6lf\n",tjd);

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
  fprintf(stdout,"\tEPOCH: Az = %lf, h = %lf\n",Az,h);

  ////////////////////////////////////////////////////
  //SHOW OBJECT POSITION
  ////////////////////////////////////////////////////
  fprintf(stdout,"Position of %s in sky w.r.t. Observer:\n",obj);
  fprintf(stdout,"\tDistance = %.17e km\n",distance);
  fprintf(stdout,"\tLight time = %.17e seconds\n",ltime);
  fprintf(stdout,"\tJ2000: RA = %s, DEC = %s\n",dec2sex(raJ2000),dec2sex(decJ2000));
  fprintf(stdout,"\tEPOCH: RA = %s, DEC = %s\n",dec2sex(ra),dec2sex(dec));
  fprintf(stdout,"\tEPOCH: Az = %s, h = %s\n",dec2sex(Az),dec2sex(h));
  fprintf(stdout,"\tEPOCH: Az = %lf, h = %lf\n",Az,h);

  ////////////////////////////////////////////////////
  //PLAIN OUTPUT
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
  fprintf(stderr,"ET,JD,Az(t),h(t),d,LT,RA(J2000),DEC(J2000),RA(t),DEC(t)\n");
  fprintf(stderr,"%.9e\n%.6lf\n",t,tjd);
  fprintf(stderr,"%.10e\n%.10e\n",Az,h);
  fprintf(stderr,"%.15e\n%.2lf\n",distance,ltime);
  fprintf(stderr,"%.10e\n%.10e\n",raJ2000,decJ2000);
  fprintf(stderr,"%.10e\n%.10e\n",ra,dec);

  return 0;
}
