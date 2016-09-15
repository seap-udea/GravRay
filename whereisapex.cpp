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
  SpiceDouble urearth[3],uvearth[3],normal[3];
  spkezr_c("EARTH",t,ECJ2000,"NONE",SSB,earth,&ltmp);

  //UNITARY POSITION AND VELOCITY VECTORS
  unorm_c(earth,urearth,&ltmp);
  unorm_c(earth+3,uvearth,&ltmp);

  //VECTOR NORMAL TO ECLIPTIC
  ucrss_c(urearth,uvearth,normal);

  ////////////////////////////////////////////////////
  //CREATING A GRID OF OBSERVERS
  ////////////////////////////////////////////////////
  struct ObserverStruct observer;
  double colat,alt=0;
  SpiceDouble urobserver[3],uvobserver[3],dv[3];
  SpiceDouble rproj,vproj;
  int nlat,nlon;
  int i,j;

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

  //............................................................
  //COSINE OF APEX CO-LATITUDE
  //............................................................
  vproj=vdot_c(urobserver,uvearth);
  
  fprintf(stdout,"Ecliptic projection: %e\n",rproj);
  fprintf(stdout,"Apex projection: %e\n",vproj);

  ////////////////////////////////////////////////////
  //PLAIN OUTPUT
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
  fprintf(stderr,"TDB,COSECL,COSAPEX\n");
  fprintf(stderr,"%.9e\n",t);
  fprintf(stderr,"%.9e\n",rproj);
  fprintf(stderr,"%.9e\n",vproj);
}
