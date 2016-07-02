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
  int qmatrix=1;

  ////////////////////////////////////////////////////
  //INPUTS
  ////////////////////////////////////////////////////
  if(argc>=2){
    strcpy(date,argv[1]);
    if(strcmp(date,"ET")==0){
      //==============================
      //GET EPHEMERIS TIME
      //==============================
      if(argc==3)
	t=atof(argv[2]);
      else argsError(argv[0],"No ET provided");
    }else if(strcmp(date,"NOMAT")==0){
      qmatrix=0;
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
  double colat,lat,lon,dlat,dlon,alt;
  SpiceDouble urobserver[3],uvobserver[3],dv[3];
  SpiceDouble rproj,vproj;
  int nlat,nlon;
  int i,j;

  //GRID PROPERTIES
  alt=0;
  nlat=nlon=50; /* It is always recommended to use the same division */
  dlat=dlon=90./nlat;

  //MATRIX FILE
  FILE* fm;
  if(qmatrix){
    fm=fopen("scratch/observers-matrices.dat","w");
    fprintf(fm,"%-06s%-06s%-12s%-12s%-12s%-12s\n","#1:i","2:j","3:lat","4:lon","5:rproj","6:vproj");
  }

  if(qmatrix)
    fprintf(stdout,"Generating matrices.\n");
  else
    fprintf(stdout,"Computing places.\n");
  i=0;
  
  double vprojantapex=1e100;
  double latantapex=0,lonantapex=0;
  double vprojapex=0;
  double latapex=0,lonapex=0;
  for(colat=0.+dlat;colat<=180.-dlat;colat+=dlat){
    lat=colat-90.;
    j=0;
    /*
      Last and beginning longitude must be repeated for mapping the
      matrices
    */
    for(lon=0.;lon<=360.;lon+=dlon){
      //............................................................
      //INITIALIZING OBSERVER 
      //............................................................
      observer.lat=lat;
      observer.lon=lon;
      observer.alt=alt;
      initObserver(t,&observer);

      //............................................................
      //COMPUTING ECLIPJ2000 POSITION AND VELOCITY OF OBSERVER
      //............................................................
      observerVelocity(&observer,0,0,0);/*Zero means a rest observer*/

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

      //ANTIAPEX
      if(vproj<=vprojantapex){
	vprojantapex=vproj;
	latantapex=lat;
	lonantapex=lon;
      }

      //APEX
      if(vproj>=vprojapex){
	vprojapex=vproj;
	latapex=lat;
	lonapex=lon;
      }

      if(qmatrix)
	fprintf(fm,"%-06d%-06d%-+12.5f%-+12.5f%-+12.3e%-+12.3e\n",
		i,j,lat,lon,rproj,vproj);
      j++;
    }
    i++;
  }
  if(qmatrix) fclose(fm);

  fprintf(stdout,"Matrix size: %d x %d\n",i,j);
  fprintf(stdout,"Grid precision: dlat = %.2lf, dlon %.2lf\n",dlat,dlon);
  fprintf(stdout,"Apex location: lat = %.5e, lon = %.5e\n",latapex,lonapex);
  fprintf(stdout,"Antiapex location: lat = %.5e, lon = %.5e\n",latantapex,lonantapex);
  fprintf(stdout,"Done.\n");

  ////////////////////////////////////////////////////
  //PLAIN OUTPUT
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
  fprintf(stderr,"TDB,NLAT,NLON,dlat,dlon,lat.apex,lon.apex,lat.antiapex,lon.antiapex\n");
  fprintf(stderr,"%.9e\n",t);
  fprintf(stderr,"%d\n%d\n",i,j);
  fprintf(stderr,"%.2lf\n%.2lf\n",dlat,dlon);
  fprintf(stderr,"%.5lf\n%.5lf\n",latapex,lonapex);
  fprintf(stderr,"%.5lf\n%.5lf\n",latantapex,lonantapex);
}
