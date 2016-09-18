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
  SpiceChar date[100];
  SpiceDouble lat,lon,alt,h,Az,v;
  if(argc==8){
    lat=atof(argv[1]);
    lon=atof(argv[2]);
    alt=atof(argv[3]);
    h=atof(argv[4]);
    Az=atof(argv[5]);
    v=atof(argv[6]);
    strcpy(date,argv[7]);
  }else argsError(argv[0]);

  ////////////////////////////////////////////////////
  //GET EPHEMERIS TIME
  ////////////////////////////////////////////////////
  SpiceDouble t,tjd,ltmp,dt;
  str2et_c(date,&t);
  //CONVERT UTC TO TDB
  deltet_c(t,"et",&dt);
  t-=dt;
  tjd=t2jd(t);
  fprintf(stdout,"Julian Date = %.6lf\n",tjd);

  ////////////////////////////////////////////////////
  //CHECK TIME
  ////////////////////////////////////////////////////
  //CHECK DATES
  SpiceDouble tref;
  if(t>=ETINI && t<=ETEND){
    tref=t;
  }else{
    SpiceChar UTC[100];
    SpiceDouble dt;
    deltet_c(t,"et",&dt);
    et2utc_c(t+dt,"ISOC",2,100,UTC);
    if(t<ETINI){
      UTC[0]='1';UTC[1]='9';UTC[2]='6';UTC[3]='3';
      str2et_c(UTC,&tref);
      fprintf(stdout,"ETEND = %.10e, t = %.10e, tref = %.10e\n",ETINI,t,tref);
    }
    else{
      UTC[0]='2';UTC[1]='0';UTC[2]='3';UTC[3]='2';
      str2et_c(UTC,&tref);
      fprintf(stdout,"ETEND = %.10e, t = %.10e, tref = %.10e\n",ETEND,t,tref);
    }
  }

  ////////////////////////////////////////////////////
  //GET OBSERVER POSITION IN TIME
  ////////////////////////////////////////////////////
  SpiceDouble earthSSBJ2000[6];
  SpiceDouble M_ITRF93_J2000[3][3];
  SpiceDouble observerITRF93[6],observerJ2000[6],observerSSBJ2000[6];
  pxform_c("ITRF93",ECJ2000,tref,M_ITRF93_J2000);
  /*
  printf("%e,%e,%e\n%e,%e,%e\n%e,%e,%e\n",
	 M_ITRF93_J2000[0][0],M_ITRF93_J2000[0][0],M_ITRF93_J2000[0][0],
	 M_ITRF93_J2000[1][0],M_ITRF93_J2000[1][0],M_ITRF93_J2000[1][0],
	 M_ITRF93_J2000[2][0],M_ITRF93_J2000[2][0],M_ITRF93_J2000[2][0]);
  */

  //OBSERVER POSITION W.R.T. EARTH CENTER IN ITRF93
  georec_c(D2R(lon),D2R(lat),alt/1000.0,REARTH,FEARTH,observerITRF93);
  fprintf(stdout,"Position observer w.r.t. ITRF93: %s\n",vec2str(observerITRF93,"%.17e"));

  ////////////////////////////////////////////////////
  //GET TOPOCENTRIC TRANSFORM MATRIX
  ////////////////////////////////////////////////////
  SpiceDouble hm[3][3],hi[3][3];
  horgeo(lat,lon,hm,hi);

  ////////////////////////////////////////////////////
  //VELOCITY OF OBSERVER
  ////////////////////////////////////////////////////

  //ROTATIONAL VELOCITY
  SpiceDouble rho,vcirc,vrot[3],vrotitrf93[3];
  rho=sqrt(observerITRF93[0]*observerITRF93[0]+observerITRF93[1]*observerITRF93[1]);
  vcirc=2*M_PI*rho/GSL_CONST_MKSA_DAY;
  //THE MINUS SIGN INDICATED 
  vpack_c(0.0,-vcirc,0.0,vrot);
  mxv_c(hi,vrot,vrotitrf93);
  fprintf(stdout,"\tVelocity of rotation w.r.t. ITRF93: %s\n",vec2str(vrotitrf93,"%.17e"));
  
  //VELOCITY OF OBSERVER IN SPACE W.R.T. TO LOCAL REFERENCE
  SpiceDouble vloc[3],vmotitrf93[3];
  SpiceDouble cA=cos(D2R(Az)),sA=sin(D2R(Az)),ch=cos(D2R(h)),sh=sin(D2R(h));
  vpack_c(v*ch*cA,-v*ch*sA,v*sh,vloc);
  fprintf(stdout,"\tVelocity observer w.r.t. LOCAL (outwards): %s\n",vec2str(vloc,"%.17e"));
  /*
    //IN CASE YOU NEED THE NEGATIVE VELOCITY
    vscl_c(-1,vloc,vloc);
    fprintf(stdout,"\tVelocity observer w.r.t. LOCAL (inwards): %s\n",vec2str(vloc,"%.17e"));
  */

  //VELOCITY OF OBSERVER IN SPACE W.R.T. TO ITRF93
  mxv_c(hi,vloc,vmotitrf93);
  fprintf(stdout,"\tVelocity of motion w.r.t. ITRF93: %s\n",vec2str(vmotitrf93,"%.17e"));

  //TOTAL VELOCITY WITH RESPECT ITRF93
  vadd_c(vrotitrf93,vmotitrf93,observerITRF93+3);
  fprintf(stdout,"\tVelocity total w.r.t. ITRF93: %s\n",vec2str(observerITRF93+3,"%.17e"));

  //OBSERVER POSITION AND VELOCITY W.R.T. EARTH CENTER IN J2000 RF
  mxv_c(M_ITRF93_J2000,observerITRF93,observerJ2000);
  fprintf(stdout,"Position observer w.r.t. J2000: %s\n",vec2str(observerJ2000,"%.17e"));
  mxv_c(M_ITRF93_J2000,observerITRF93+3,observerJ2000+3);
  fprintf(stdout,"\tVelocity observer w.r.t. J2000: %s\n",vec2str(observerJ2000+3,"%.17e"));

  //EARTH POSITION W.R.T. SOLAR SYSTEM BARYCENTER IN J2000 RF
  spkezr_c(EARTH_ID,t,ECJ2000,"NONE","SOLAR SYSTEM BARYCENTER",earthSSBJ2000,&ltmp);
  fprintf(stdout,"Position earth w.r.t. SSB J2000: %s\n",vec2str(earthSSBJ2000,"%.17e"));
  fprintf(stdout,"\tVelocity earth w.r.t. SSB J2000: %s\n",vec2str(earthSSBJ2000+3,"%.17e"));

  //OBSERVER POSITION W.R.T. SOLAR SYSTEM BARYCENTER IN J2000 RF
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);
  fprintf(stdout,"Position observer w.r.t. SSB J2000: %s\n",vec2str(observerSSBJ2000,"%.17e"));
  vadd_c(earthSSBJ2000+3,observerJ2000+3,observerSSBJ2000+3);
  fprintf(stdout,"\tVelocity observer w.r.t. SSB J2000: %s\n",vec2str(observerSSBJ2000+3,"%.17e"));

  ////////////////////////////////////////////////////
  //PLAIN INFORMATION
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
  fprintf(stderr,"TDB,JD,DT,ITRF93(6),ECJ2000(6)\n");
  fprintf(stderr,"%.9e\n%.6lf\n%.2lf\n",t,tjd,dt);
  fprintf(stderr,"%s\n",vec2strn(observerITRF93,6,"%+.17e "));
  fprintf(stderr,"%s\n",vec2strn(observerSSBJ2000,6,"%+.17e "));
  return 0;
}
