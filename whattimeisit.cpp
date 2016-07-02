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
  SpiceChar date[100],outsys[100];
  if(argc==3){
    strcpy(date,argv[1]);
    strcpy(outsys,argv[2]);
    fprintf(stdout,"Input date: %s\n",date);
  }else argsError(argv[0]);

  ////////////////////////////////////////////////////
  //GET EPHEMERIS TIME
  ////////////////////////////////////////////////////
  SpiceDouble et,tdb,t,tjd,tjdb,ltmp,dt;

  //EPHEMERIS TIME
  str2et_c(date,&et);
  fprintf(stdout,"ET = %.9e\n",et);
  tjd=t2jd(et);
  fprintf(stdout,"Julian Date at ET = %.6lf\n",tjd);

  //TDB
  deltet_c(et,"et",&dt);
  tdb=et-dt;
  fprintf(stdout,"TDB = %.9e\n",tdb);
  tjdb=t2jd(tdb);
  fprintf(stdout,"Julian Date at TDB = %.6lf\n",tjdb);

  //CUSTOM SYSTEM
  t=unitim_c(et,"ET",outsys);
  fprintf(stdout,"Custom system (%s) = %.17e\n",outsys,t);

  ////////////////////////////////////////////////////
  //PLAIN OUTPUT
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
  fprintf(stderr,"ET,JD,DT,TDB,JDB,%s\n",outsys);
  fprintf(stderr,"%.9e\n",et);
  fprintf(stderr,"%.6lf\n",tjd);
  fprintf(stderr,"%.2lf\n",dt);
  fprintf(stderr,"%.9e\n",tdb);
  fprintf(stderr,"%.6lf\n",tjdb);
  fprintf(stderr,"%.17e\n",t);
  return 0;
}
