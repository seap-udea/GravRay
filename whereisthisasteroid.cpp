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
      to the Solar System Barycenter of an asteroid

    Arguments are: 

      asteroid, dat

    Date format:

       MM/DD/CCYY HH:MM:SS.dcm UTC-L

    Example:

       ./whereisthisasteroid EROS "07/19/2015 00:00:00.000 UTC-5"

    Available asteroids see list at the end

  */
  SpiceChar date[100],obj[100];
  SpiceDouble t,tjd,ltmp,dt;

  strcpy(obj,argv[1]);
  strcpy(date,argv[2]);
  if(strcmp(date,"ET")==0){
    t=atof(argv[3]);
  }else{
    ////////////////////////////////////////////////////
    //GET EPHEMERIS TIME
    ////////////////////////////////////////////////////
    str2et_c(date,&t);
    deltet_c(t,"et",&dt);
    fprintf(stderr,"DT = %.2lf\n",dt);
    fprintf(stderr,"TT = %.9e\n",t);
    t-=dt;
  }
  //CORRECT TIME FOR DELTAT.  NOW T IS TDB
  fprintf(stderr,"TDB = ");
  fprintf(stdout,"%.9e ",t);

  //JULIAN DATE TO VERIFY
  tjd=t2jd(t);
  fprintf(stderr,"\nJulian Date = %.6lf\n",tjd);

  ////////////////////////////////////////////////////
  //GET POSITION AT t
  ////////////////////////////////////////////////////
  SpiceDouble objectSSBJ2000[6],objectSSBEJ2000[6],objectSUNJ2000[6],sunSSBJ2000[6];
  spkezr_c(obj,t,"ECLIPJ2000_DE405","NONE","SUN",objectSUNJ2000,&ltmp);
  spkezr_c("SUN",t,"ECLIPJ2000_DE405","NONE",SSB,sunSSBJ2000,&ltmp);
  vadd_c(objectSUNJ2000,sunSSBJ2000,objectSSBJ2000);
  vadd_c(objectSUNJ2000+3,sunSSBJ2000+3,objectSSBJ2000+3);
  
  SpiceDouble M_405_ABS[3][3];
  pxform_c("ECLIPJ2000_DE405","ECLIPJ2000",t,M_405_ABS);
  mxv_c(M_405_ABS,objectSSBJ2000,objectSSBEJ2000);
  mxv_c(M_405_ABS,objectSSBJ2000+3,objectSSBEJ2000+3);

  fprintf(stderr,"State vector of %s: ",obj);
  fprintf(stdout,"%s",vec2strn(objectSSBEJ2000,6,"%.17e "));
  fprintf(stderr,"\nLT = %.17e\n",ltmp);
}

/*

List of asteroids:

2000001 CERES* w.r.t. 10 SUN      1799 DEC 30 12:00:00.000        2199 DEC 13 12:00:00.000
2000002* w.r.t. 10 SUN                        Same coverage as previous object
2000003* w.r.t. 10 SUN                        Same coverage as previous object
2000004 VESTA* w.r.t. 10 SUN                  Same coverage as previous object
2000005* w.r.t. 10 SUN                        Same coverage as previous object
2000006* w.r.t. 10 SUN                        Same coverage as previous object
2000007* w.r.t. 10 SUN                        Same coverage as previous object
2000008* w.r.t. 10 SUN                        Same coverage as previous object
2000009* w.r.t. 10 SUN                        Same coverage as previous object
2000010* w.r.t. 10 SUN                        Same coverage as previous object
2000011* w.r.t. 10 SUN                        Same coverage as previous object
2000012* w.r.t. 10 SUN                        Same coverage as previous object
2000013* w.r.t. 10 SUN                        Same coverage as previous object
2000014* w.r.t. 10 SUN                        Same coverage as previous object
2000015* w.r.t. 10 SUN                        Same coverage as previous object
2000016* w.r.t. 10 SUN                        Same coverage as previous object
2000017* w.r.t. 10 SUN                        Same coverage as previous object
2000018* w.r.t. 10 SUN                        Same coverage as previous object
2000019* w.r.t. 10 SUN                        Same coverage as previous object
2000020* w.r.t. 10 SUN                        Same coverage as previous object
2000021 LUTETIA* w.r.t. 10 SUN                Same coverage as previous object
2000022* w.r.t. 10 SUN                        Same coverage as previous object
2000023* w.r.t. 10 SUN                        Same coverage as previous object
2000024* w.r.t. 10 SUN                        Same coverage as previous object
2000025* w.r.t. 10 SUN                        Same coverage as previous object
2000026* w.r.t. 10 SUN                        Same coverage as previous object
2000027* w.r.t. 10 SUN                        Same coverage as previous object
2000028* w.r.t. 10 SUN                        Same coverage as previous object
2000029* w.r.t. 10 SUN                        Same coverage as previous object
2000030* w.r.t. 10 SUN                        Same coverage as previous object
2000031* w.r.t. 10 SUN                        Same coverage as previous object
2000032* w.r.t. 10 SUN                        Same coverage as previous object
2000034* w.r.t. 10 SUN                        Same coverage as previous object
2000035* w.r.t. 10 SUN                        Same coverage as previous object
2000036* w.r.t. 10 SUN                        Same coverage as previous object
2000037* w.r.t. 10 SUN                        Same coverage as previous object
2000038* w.r.t. 10 SUN                        Same coverage as previous object
2000039* w.r.t. 10 SUN                        Same coverage as previous object
2000040* w.r.t. 10 SUN                        Same coverage as previous object
2000041* w.r.t. 10 SUN                        Same coverage as previous object
2000042* w.r.t. 10 SUN                        Same coverage as previous object
2000043* w.r.t. 10 SUN                        Same coverage as previous object
2000044* w.r.t. 10 SUN                        Same coverage as previous object
2000045* w.r.t. 10 SUN                        Same coverage as previous object
2000046* w.r.t. 10 SUN                        Same coverage as previous object
2000047* w.r.t. 10 SUN                        Same coverage as previous object
2000048* w.r.t. 10 SUN                        Same coverage as previous object
2000049* w.r.t. 10 SUN                        Same coverage as previous object
2000050* w.r.t. 10 SUN                        Same coverage as previous object
2000051* w.r.t. 10 SUN                        Same coverage as previous object
2000052* w.r.t. 10 SUN                        Same coverage as previous object
2000053* w.r.t. 10 SUN                        Same coverage as previous object
2000054* w.r.t. 10 SUN                        Same coverage as previous object
2000056* w.r.t. 10 SUN                        Same coverage as previous object
2000057* w.r.t. 10 SUN                        Same coverage as previous object
2000058* w.r.t. 10 SUN                        Same coverage as previous object
2000059* w.r.t. 10 SUN                        Same coverage as previous object
2000062* w.r.t. 10 SUN                        Same coverage as previous object
2000063* w.r.t. 10 SUN                        Same coverage as previous object
2000065* w.r.t. 10 SUN                        Same coverage as previous object
2000068* w.r.t. 10 SUN                        Same coverage as previous object
2000069* w.r.t. 10 SUN                        Same coverage as previous object
2000070* w.r.t. 10 SUN                        Same coverage as previous object
2000071* w.r.t. 10 SUN                        Same coverage as previous object
2000072* w.r.t. 10 SUN                        Same coverage as previous object
2000074* w.r.t. 10 SUN                        Same coverage as previous object
2000075* w.r.t. 10 SUN                        Same coverage as previous object
2000076* w.r.t. 10 SUN                        Same coverage as previous object
2000077* w.r.t. 10 SUN                        Same coverage as previous object
2000078* w.r.t. 10 SUN                        Same coverage as previous object
2000080* w.r.t. 10 SUN                        Same coverage as previous object
2000081* w.r.t. 10 SUN                        Same coverage as previous object
2000083* w.r.t. 10 SUN                        Same coverage as previous object
2000084* w.r.t. 10 SUN                        Same coverage as previous object
2000085* w.r.t. 10 SUN                        Same coverage as previous object
2000086* w.r.t. 10 SUN                        Same coverage as previous object
2000087* w.r.t. 10 SUN                        Same coverage as previous object
2000088* w.r.t. 10 SUN                        Same coverage as previous object
2000089* w.r.t. 10 SUN                        Same coverage as previous object
2000090* w.r.t. 10 SUN                        Same coverage as previous object
2000091* w.r.t. 10 SUN                        Same coverage as previous object
2000092* w.r.t. 10 SUN                        Same coverage as previous object
2000093* w.r.t. 10 SUN                        Same coverage as previous object
2000094* w.r.t. 10 SUN                        Same coverage as previous object
2000095* w.r.t. 10 SUN                        Same coverage as previous object
2000096* w.r.t. 10 SUN                        Same coverage as previous object
2000097* w.r.t. 10 SUN                        Same coverage as previous object
2000098* w.r.t. 10 SUN                        Same coverage as previous object
2000099* w.r.t. 10 SUN                        Same coverage as previous object
2000102* w.r.t. 10 SUN                        Same coverage as previous object
2000103* w.r.t. 10 SUN                        Same coverage as previous object
2000104* w.r.t. 10 SUN                        Same coverage as previous object
2000105* w.r.t. 10 SUN                        Same coverage as previous object
2000106* w.r.t. 10 SUN                        Same coverage as previous object
2000107* w.r.t. 10 SUN                        Same coverage as previous object
2000109* w.r.t. 10 SUN                        Same coverage as previous object
2000110* w.r.t. 10 SUN                        Same coverage as previous object
2000111* w.r.t. 10 SUN                        Same coverage as previous object
2000112* w.r.t. 10 SUN                        Same coverage as previous object
2000114* w.r.t. 10 SUN                        Same coverage as previous object
2000115* w.r.t. 10 SUN                        Same coverage as previous object
2000117* w.r.t. 10 SUN                        Same coverage as previous object
2000120* w.r.t. 10 SUN                        Same coverage as previous object
2000121* w.r.t. 10 SUN                        Same coverage as previous object
2000124* w.r.t. 10 SUN                        Same coverage as previous object
2000127* w.r.t. 10 SUN                        Same coverage as previous object
2000128* w.r.t. 10 SUN                        Same coverage as previous object
2000129* w.r.t. 10 SUN                        Same coverage as previous object
2000130* w.r.t. 10 SUN                        Same coverage as previous object
2000134* w.r.t. 10 SUN                        Same coverage as previous object
2000135* w.r.t. 10 SUN                        Same coverage as previous object
2000137* w.r.t. 10 SUN                        Same coverage as previous object
2000139* w.r.t. 10 SUN                        Same coverage as previous object
2000140* w.r.t. 10 SUN                        Same coverage as previous object
2000141* w.r.t. 10 SUN                        Same coverage as previous object
2000143* w.r.t. 10 SUN                        Same coverage as previous object
2000144* w.r.t. 10 SUN                        Same coverage as previous object
2000145* w.r.t. 10 SUN                        Same coverage as previous object
2000146* w.r.t. 10 SUN                        Same coverage as previous object
2000147* w.r.t. 10 SUN                        Same coverage as previous object
2000148* w.r.t. 10 SUN                        Same coverage as previous object
2000150* w.r.t. 10 SUN                        Same coverage as previous object
2000154* w.r.t. 10 SUN                        Same coverage as previous object
2000156* w.r.t. 10 SUN                        Same coverage as previous object
2000159* w.r.t. 10 SUN                        Same coverage as previous object
2000160* w.r.t. 10 SUN                        Same coverage as previous object
2000162* w.r.t. 10 SUN                        Same coverage as previous object
2000163* w.r.t. 10 SUN                        Same coverage as previous object
2000164* w.r.t. 10 SUN                        Same coverage as previous object
2000165* w.r.t. 10 SUN                        Same coverage as previous object
2000168* w.r.t. 10 SUN                        Same coverage as previous object
2000171* w.r.t. 10 SUN                        Same coverage as previous object
2000173* w.r.t. 10 SUN                        Same coverage as previous object
2000175* w.r.t. 10 SUN                        Same coverage as previous object
2000176* w.r.t. 10 SUN                        Same coverage as previous object
2000181* w.r.t. 10 SUN                        Same coverage as previous object
2000185* w.r.t. 10 SUN                        Same coverage as previous object
2000187* w.r.t. 10 SUN                        Same coverage as previous object
2000191* w.r.t. 10 SUN                        Same coverage as previous object
2000192* w.r.t. 10 SUN                        Same coverage as previous object
2000194* w.r.t. 10 SUN                        Same coverage as previous object
2000195* w.r.t. 10 SUN                        Same coverage as previous object
2000196* w.r.t. 10 SUN                        Same coverage as previous object
2000200* w.r.t. 10 SUN                        Same coverage as previous object
2000201* w.r.t. 10 SUN                        Same coverage as previous object
2000203* w.r.t. 10 SUN                        Same coverage as previous object
2000205* w.r.t. 10 SUN                        Same coverage as previous object
2000206* w.r.t. 10 SUN                        Same coverage as previous object
2000209* w.r.t. 10 SUN                        Same coverage as previous object
2000210* w.r.t. 10 SUN                        Same coverage as previous object
2000211* w.r.t. 10 SUN                        Same coverage as previous object
2000212* w.r.t. 10 SUN                        Same coverage as previous object
2000213* w.r.t. 10 SUN                        Same coverage as previous object
2000216 KLEOPATRA* w.r.t. 10 SUN              Same coverage as previous object
2000221* w.r.t. 10 SUN                        Same coverage as previous object
2000224* w.r.t. 10 SUN                        Same coverage as previous object
2000225* w.r.t. 10 SUN                        Same coverage as previous object
2000230* w.r.t. 10 SUN                        Same coverage as previous object
2000233* w.r.t. 10 SUN                        Same coverage as previous object
2000236* w.r.t. 10 SUN                        Same coverage as previous object
2000238* w.r.t. 10 SUN                        Same coverage as previous object
2000240* w.r.t. 10 SUN                        Same coverage as previous object
2000241* w.r.t. 10 SUN                        Same coverage as previous object
2000247* w.r.t. 10 SUN                        Same coverage as previous object
2000250* w.r.t. 10 SUN                        Same coverage as previous object
2000259* w.r.t. 10 SUN                        Same coverage as previous object
2000266* w.r.t. 10 SUN                        Same coverage as previous object
2000268* w.r.t. 10 SUN                        Same coverage as previous object
2000275* w.r.t. 10 SUN                        Same coverage as previous object
2000276* w.r.t. 10 SUN                        Same coverage as previous object
2000283* w.r.t. 10 SUN                        Same coverage as previous object
2000287* w.r.t. 10 SUN                        Same coverage as previous object
2000303* w.r.t. 10 SUN                        Same coverage as previous object
2000304* w.r.t. 10 SUN                        Same coverage as previous object
2000308* w.r.t. 10 SUN                        Same coverage as previous object
2000313* w.r.t. 10 SUN                        Same coverage as previous object
2000322* w.r.t. 10 SUN                        Same coverage as previous object
2000324* w.r.t. 10 SUN                        Same coverage as previous object
2000326* w.r.t. 10 SUN                        Same coverage as previous object
2000328* w.r.t. 10 SUN                        Same coverage as previous object
2000329* w.r.t. 10 SUN                        Same coverage as previous object
2000334* w.r.t. 10 SUN                        Same coverage as previous object
2000335* w.r.t. 10 SUN                        Same coverage as previous object
2000336* w.r.t. 10 SUN                        Same coverage as previous object
2000337* w.r.t. 10 SUN                        Same coverage as previous object
2000338* w.r.t. 10 SUN                        Same coverage as previous object
2000344* w.r.t. 10 SUN                        Same coverage as previous object
2000345* w.r.t. 10 SUN                        Same coverage as previous object
2000346* w.r.t. 10 SUN                        Same coverage as previous object
2000347* w.r.t. 10 SUN                        Same coverage as previous object
2000349* w.r.t. 10 SUN                        Same coverage as previous object
2000350* w.r.t. 10 SUN                        Same coverage as previous object
2000354* w.r.t. 10 SUN                        Same coverage as previous object
2000356* w.r.t. 10 SUN                        Same coverage as previous object
2000357* w.r.t. 10 SUN                        Same coverage as previous object
2000358* w.r.t. 10 SUN                        Same coverage as previous object
2000360* w.r.t. 10 SUN                        Same coverage as previous object
2000362* w.r.t. 10 SUN                        Same coverage as previous object
2000363* w.r.t. 10 SUN                        Same coverage as previous object
2000365* w.r.t. 10 SUN                        Same coverage as previous object
2000366* w.r.t. 10 SUN                        Same coverage as previous object
2000369* w.r.t. 10 SUN                        Same coverage as previous object
2000372* w.r.t. 10 SUN                        Same coverage as previous object
2000373* w.r.t. 10 SUN                        Same coverage as previous object
2000375* w.r.t. 10 SUN                        Same coverage as previous object
2000377* w.r.t. 10 SUN                        Same coverage as previous object
2000381* w.r.t. 10 SUN                        Same coverage as previous object
2000385* w.r.t. 10 SUN                        Same coverage as previous object
2000386* w.r.t. 10 SUN                        Same coverage as previous object
2000387* w.r.t. 10 SUN                        Same coverage as previous object
2000388* w.r.t. 10 SUN                        Same coverage as previous object
2000389* w.r.t. 10 SUN                        Same coverage as previous object
2000393* w.r.t. 10 SUN                        Same coverage as previous object
2000404* w.r.t. 10 SUN                        Same coverage as previous object
2000405* w.r.t. 10 SUN                        Same coverage as previous object
2000407* w.r.t. 10 SUN                        Same coverage as previous object
2000409* w.r.t. 10 SUN                        Same coverage as previous object
2000410* w.r.t. 10 SUN                        Same coverage as previous object
2000412* w.r.t. 10 SUN                        Same coverage as previous object
2000416* w.r.t. 10 SUN                        Same coverage as previous object
2000419* w.r.t. 10 SUN                        Same coverage as previous object
2000420* w.r.t. 10 SUN                        Same coverage as previous object
2000423* w.r.t. 10 SUN                        Same coverage as previous object
2000424* w.r.t. 10 SUN                        Same coverage as previous object
2000426* w.r.t. 10 SUN                        Same coverage as previous object
2000431* w.r.t. 10 SUN                        Same coverage as previous object
2000433 EROS* w.r.t. 10 SUN                   Same coverage as previous object
2000442* w.r.t. 10 SUN                        Same coverage as previous object
2000444* w.r.t. 10 SUN                        Same coverage as previous object
2000449* w.r.t. 10 SUN                        Same coverage as previous object
2000451* w.r.t. 10 SUN                        Same coverage as previous object
2000454* w.r.t. 10 SUN                        Same coverage as previous object
2000455* w.r.t. 10 SUN                        Same coverage as previous object
2000466* w.r.t. 10 SUN                        Same coverage as previous object
2000469* w.r.t. 10 SUN                        Same coverage as previous object
2000471* w.r.t. 10 SUN                        Same coverage as previous object
2000476* w.r.t. 10 SUN                        Same coverage as previous object
2000481* w.r.t. 10 SUN                        Same coverage as previous object
2000488* w.r.t. 10 SUN                        Same coverage as previous object
2000489* w.r.t. 10 SUN                        Same coverage as previous object
2000490* w.r.t. 10 SUN                        Same coverage as previous object
2000491* w.r.t. 10 SUN                        Same coverage as previous object
2000498* w.r.t. 10 SUN                        Same coverage as previous object
2000505* w.r.t. 10 SUN                        Same coverage as previous object
2000506* w.r.t. 10 SUN                        Same coverage as previous object
2000508* w.r.t. 10 SUN                        Same coverage as previous object
2000511* w.r.t. 10 SUN                        Same coverage as previous object
2000514* w.r.t. 10 SUN                        Same coverage as previous object
2000521* w.r.t. 10 SUN                        Same coverage as previous object
2000532* w.r.t. 10 SUN                        Same coverage as previous object
2000535* w.r.t. 10 SUN                        Same coverage as previous object
2000536* w.r.t. 10 SUN                        Same coverage as previous object
2000545* w.r.t. 10 SUN                        Same coverage as previous object
2000554* w.r.t. 10 SUN                        Same coverage as previous object
2000566* w.r.t. 10 SUN                        Same coverage as previous object
2000568* w.r.t. 10 SUN                        Same coverage as previous object
2000595* w.r.t. 10 SUN                        Same coverage as previous object
2000596* w.r.t. 10 SUN                        Same coverage as previous object
2000602* w.r.t. 10 SUN                        Same coverage as previous object
2000618* w.r.t. 10 SUN                        Same coverage as previous object
2000626* w.r.t. 10 SUN                        Same coverage as previous object
2000635* w.r.t. 10 SUN                        Same coverage as previous object
2000654* w.r.t. 10 SUN                        Same coverage as previous object
2000663* w.r.t. 10 SUN                        Same coverage as previous object
2000674* w.r.t. 10 SUN                        Same coverage as previous object
2000683* w.r.t. 10 SUN                        Same coverage as previous object
2000690* w.r.t. 10 SUN                        Same coverage as previous object
2000691* w.r.t. 10 SUN                        Same coverage as previous object
2000694* w.r.t. 10 SUN                        Same coverage as previous object
2000702* w.r.t. 10 SUN                        Same coverage as previous object
2000704* w.r.t. 10 SUN                        Same coverage as previous object
2000705* w.r.t. 10 SUN                        Same coverage as previous object
2000709* w.r.t. 10 SUN                        Same coverage as previous object
2000712* w.r.t. 10 SUN                        Same coverage as previous object
2000713* w.r.t. 10 SUN                        Same coverage as previous object
2000739* w.r.t. 10 SUN                        Same coverage as previous object
2000740* w.r.t. 10 SUN                        Same coverage as previous object
2000747* w.r.t. 10 SUN                        Same coverage as previous object
2000751* w.r.t. 10 SUN                        Same coverage as previous object
2000762* w.r.t. 10 SUN                        Same coverage as previous object
2000769* w.r.t. 10 SUN                        Same coverage as previous object
2000772* w.r.t. 10 SUN                        Same coverage as previous object
2000773* w.r.t. 10 SUN                        Same coverage as previous object
2000776* w.r.t. 10 SUN                        Same coverage as previous object
2000780* w.r.t. 10 SUN                        Same coverage as previous object
2000788* w.r.t. 10 SUN                        Same coverage as previous object
2000790* w.r.t. 10 SUN                        Same coverage as previous object
2000791* w.r.t. 10 SUN                        Same coverage as previous object
2000804* w.r.t. 10 SUN                        Same coverage as previous object
2000814* w.r.t. 10 SUN                        Same coverage as previous object
2000849* w.r.t. 10 SUN                        Same coverage as previous object
2000895* w.r.t. 10 SUN                        Same coverage as previous object
2000909* w.r.t. 10 SUN                        Same coverage as previous object
2000914* w.r.t. 10 SUN                        Same coverage as previous object
2000980* w.r.t. 10 SUN                        Same coverage as previous object
2001015* w.r.t. 10 SUN                        Same coverage as previous object
2001021* w.r.t. 10 SUN                        Same coverage as previous object
2001036* w.r.t. 10 SUN                        Same coverage as previous object
2001093* w.r.t. 10 SUN                        Same coverage as previous object
2001467* w.r.t. 10 SUN                        Same coverage as previous object

*/
