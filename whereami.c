#include <gravray.c>

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
      Calculate the position and velocity of an observer at a given
      geographic position and having a given velocity in space with
      respect to the local reference frame.

    Arguments are: 

      latitude (degrees), longitude (degrees), elevation (meters),
      elevation (degrees), azimuth (degrees), velocity (km/s), date

    Date format:
       MM/DD/CCYY HH:MM:SS.dcm UTC-L

    Example:

       ./whereami.exe 6.2 -75.34 1450.0 45.0 0.0 1.0 "07/19/2015 00:00:00.000 UTC-5"

       ./whereami.exe 0.0 0.0 1560.0 45.0 0.0 1.0 "07/19/2015 10:00:00.000 UTC"
  */
  SpiceChar date[100];
  SpiceDouble lat=atof(argv[1]);
  SpiceDouble lon=atof(argv[2]);
  SpiceDouble alt=atof(argv[3]);
  SpiceDouble h=atof(argv[4]);
  SpiceDouble Az=atof(argv[5]);
  SpiceDouble v=atof(argv[6]);
  strcpy(date,argv[7]);

  ////////////////////////////////////////////////////
  //GET EPHEMERIS TIME
  ////////////////////////////////////////////////////
  SpiceDouble t,tjd,ltmp;
  str2et_c(date,&t);
  tjd=t2jd(t,0);
  printf("jd = %.17e\n",tjd);

  ////////////////////////////////////////////////////
  //GET OBSERVER POSITION IN TIME
  ////////////////////////////////////////////////////
  SpiceDouble earthSSBJ2000[6];
  SpiceDouble M_J2000_Epoch[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  SpiceDouble M_ITRF93_J2000[3][3];
  SpiceDouble observerITRF93[6],observerJ2000[6],observerSSBJ2000[6],observerSSBEpoch[6];

  pxform_c("J2000","EARTHTRUEEPOCH",t,M_J2000_Epoch);
  pxform_c("ITRF93","J2000",t,M_ITRF93_J2000);

  //OBSERVER POSITION W.R.T. EARTH CENTER IN ITRF93
  georec_c(D2R(lon),D2R(lat),alt/1000.0,REARTH,FEARTH,observerITRF93);
  printf("Position observer w.r.t. ITRF93: %s\n",vec2str(observerITRF93,"%.17e"));

  ////////////////////////////////////////////////////
  //GET TOPOCENTRIC TRANSFORM MATRIX
  ////////////////////////////////////////////////////
  SpiceDouble hm[3][3],hi[3][3];
  hormat(lat,lon,hm,hi);

  ////////////////////////////////////////////////////
  //VELOCITY OF OBSERVER
  ////////////////////////////////////////////////////

  //ROTATIONAL VELOCITY
  SpiceDouble rho,vcirc,vrot[3],vrotitrf93[3];
  rho=sqrt(observerITRF93[0]*observerITRF93[0]+observerITRF93[1]*observerITRF93[1]);
  vcirc=2*M_PI*rho/GSL_CONST_MKSA_DAY*0;
  vpack_c(0.0,-vcirc,0.0,vrot);
  mxv_c(hi,vrot,vrotitrf93);
  printf("\tVelocity of rotation w.r.t. ITRF93: %s\n",vec2str(vrotitrf93,"%.17e"));
  
  //VELOCITY OF OBSERVER IN SPACE
  SpiceDouble vloc[3],vmotitrf93[3];
  SpiceDouble cA=cos(D2R(Az)),sA=sin(D2R(Az)),ch=cos(D2R(h)),sh=sin(D2R(h));
  vpack_c(v*ch*cA,-v*ch*sA,v*sh,vloc);
  printf("\tVelocity observer w.r.t. LOCAL: %s\n",vec2str(vloc,"%.17e"));
  mxv_c(hi,vloc,vmotitrf93);
  printf("\tVelocity of motion w.r.t. ITRF93: %s\n",vec2str(vmotitrf93,"%.17e"));

  //TOTAL VELOCITY WITH RESPECT ITRF93
  vadd_c(vrotitrf93,vmotitrf93,observerITRF93+3);
  printf("\tVelocity total w.r.t. ITRF93: %s\n",vec2str(observerITRF93+3,"%.17e"));


  //OBSERVER POSITION AND VELOCITY W.R.T. EARTH CENTER IN J2000 RF
  mxv_c(M_ITRF93_J2000,observerITRF93,observerJ2000);
  printf("Position observer w.r.t. J2000: %s\n",vec2str(observerJ2000,"%.17e"));

  mxv_c(M_ITRF93_J2000,observerITRF93+3,observerJ2000+3);
  printf("\tVelocity observer w.r.t. J2000: %s\n",vec2str(observerJ2000+3,"%.17e"));

  //*
  SpiceDouble vsky[3];
  ter2cel(tjd,0.0,0.0,1,0,1,0.0,0.0,vloc,vsky);
  printf("\tVelocity observer w.r.t. CELESTRIAL: %s\n",vec2str(vsky,"%.17e"));
  //*/
 
  exit(0);
  

  //EARTH POSITION W.R.T. SOLAR SYSTEM BARYCENTER IN J2000 RF
  spkezr_c(EARTH_ID,t,"J2000","NONE","SOLAR SYSTEM BARYCENTER",earthSSBJ2000,&ltmp);
  printf("Position earth w.r.t. SSB J2000: %s\n",vec2str(earthSSBJ2000,"%.17e"));
  printf("\tVelocity earth w.r.t. SSB J2000: %s\n",vec2str(earthSSBJ2000+3,"%.17e"));

  //OBSERVER POSITION W.R.T. SOLAR SYSTEM BARYCENTER IN J2000 RF
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);
  printf("Position observer w.r.t. SSB J2000: %s\n",vec2str(observerSSBJ2000,"%.17e"));
  vadd_c(earthSSBJ2000+3,observerJ2000+3,observerSSBJ2000+3);
  printf("\tVelocity observer w.r.t. SSB J2000: %s\n",vec2str(observerSSBJ2000+3,"%.17e"));

  //OBSERVER POSITION W.R.T. SOLAR SYSTEM BARYCENTER IN EPOCH RF
  mxv_c(M_J2000_Epoch,observerSSBJ2000,observerSSBEpoch);
  printf("Position observer w.r.t. SSB Epoch: %s\n",vec2str(observerSSBEpoch,"%.17e"));
  mxv_c(M_J2000_Epoch,observerSSBJ2000+3,observerSSBEpoch+3);
  printf("\tVelocity observer w.r.t. SSB Epoch: %s\n",vec2str(observerSSBEpoch+3,"%.17e"));

  return 0;
}
