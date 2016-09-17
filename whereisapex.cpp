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
  SpiceDouble t,tjd,ltmp,dt;
  SpiceDouble lat,lon;
  int qmatrix=1;

  ////////////////////////////////////////////////////
  //INPUTS
  ////////////////////////////////////////////////////
  int iarg=1;
  if(argc>=4){
    strcpy(date,argv[iarg++]);
    if(strcmp(date,"ET")==0){
      //==============================
      //GET EPHEMERIS TIME
      //==============================
      if(argc==5)
	t=atof(argv[iarg++]);
      else argsError(argv[0],"No ET provided\n");
    }else{
      //==============================
      //GET DATA STRINGS
      //==============================
      str2et_c(date,&t);
      deltet_c(t,"et",&dt);
      fprintf(stdout,"DT = %.2lf\n",dt);
      fprintf(stdout,"TT = %.9e\n",t);
      t-=dt;
    }
    lat=atof(argv[iarg++]);
    lon=atof(argv[iarg++]);
  }else
    argsError(argv[0]);

  fprintf(stdout,"TDB = %.9e\n",t);

  ////////////////////////////////////////////////////
  //EARTH POSITION
  ////////////////////////////////////////////////////
  SpiceDouble earth[6];
  SpiceDouble urearth[3],uvearth[3],normal[3],uantisun[3];
  struct ObserverStruct observer;
  double colat,alt=0;
  SpiceDouble urobserver[3],uvobserver[3],dv[3];
  SpiceDouble rproj,vproj,aproj;
  int nlat,nlon;
  int i,j;

  spkezr_c("EARTH",t,ECJ2000,"NONE",SSB,earth,&ltmp);

  //UNITARY POSITION AND VELOCITY VECTORS
  unorm_c(earth,urearth,&ltmp);
  unorm_c(earth+3,uvearth,&ltmp);

  //VECTOR NORMAL TO ECLIPTIC
  ucrss_c(urearth,uvearth,normal);

  //CROSS PRODUCT OF THE ECLIPTIC AND NORMAL TO ECLIPTIC VECTORS
  ucrss_c(uvearth,normal,uantisun);
  vproj=vdot_c(urearth,uantisun);
  printf("Antisun: %s\n",vec2str(uantisun));
  printf("rearth: %s\n",vec2str(urearth));
  printf("Difference: %e\n",R2D(acos(vproj)));

  ////////////////////////////////////////////////////
  //CREATING A GRID OF OBSERVERS
  ////////////////////////////////////////////////////

  observer.lat=lat;
  observer.lon=lon;
  observer.alt=alt;
  initObserver(t,&observer);

  //............................................................
  //COMPUTING ECLIPJ2000 POSITION AND VELOCITY OF OBSERVER
  //............................................................
  observerVelocity(&observer,0,0,0);//Zero means a rest observer
  
  //UNIT VECTORS RESPECT EARTH CENTER OF THE OBSERVER
  unorm_c(observer.posj2000,urobserver,&ltmp);
  unorm_c(observer.posj2000+3,uvobserver,&ltmp);
  
  //............................................................
  //COSINE ECLIPTIC CO-LATITUDE
  //............................................................
  rproj=vdot_c(urobserver,normal);
  /*
    This is also the cosine of the angle of the point with respect to
    apex based x-axis.
   */

  //............................................................
  //COSINE OF APEX CO-LATITUDE
  //............................................................
  vproj=vdot_c(urobserver,uvearth);
  /*
    This is also the cosine of the angle of the point with respect to
    apex based z-axis.
   */

  //............................................................
  //COSINE OF ANTISUN COLATITUDE
  //............................................................
  aproj=vdot_c(urobserver,uantisun);
  /*
    This is also the cosine of the angle of the point with respect to
    apex based y-axis.
   */
  fprintf(stdout,"Ecliptic projection: %e\n",rproj);
  fprintf(stdout,"Apex projection: %e\n",vproj);

  //............................................................
  //SPHERICAL COORDINATES OF APEX
  //............................................................
  double xapex,yapex,rhoapex;
  double latecl,latapex,lonapex;
  latecl=90.0-R2D(acos(rproj));
  latapex=90.0-R2D(acos(vproj));
  lonapex=R2D(atan2(aproj,rproj));

  ////////////////////////////////////////////////////
  //PLAIN OUTPUT
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
  fprintf(stderr,"TDB,COSECL,COSAPEX,COSASUN,LATECL,LATAPEX,LONAPEX\n");
  fprintf(stderr,"%.9e\n",t);
  fprintf(stderr,"%.9e\n",rproj);
  fprintf(stderr,"%.9e\n",vproj);
  fprintf(stderr,"%.9e\n",aproj);
  fprintf(stderr,"%.9e\n",latecl);
  fprintf(stderr,"%.9e\n",latapex);
  fprintf(stderr,"%.9e\n",lonapex);

  /*
    Apex Reference System
    ---------------------
    
    Axis: 

      zapex: Direction of Earth's velocity
      xapex: North of the ecliptic
      yapex: Antisun direction

    Special points:

      (latapex,lonapex)=(90,0) : Apex
      (latapex,lonapex)=(0,0) : Northest point of Earth above ecliptic
      (latapex,lonapex)=(0,+90) : Antisun
      (latapex,lonapex)=(0,-90) : Solar direction
  */
}
