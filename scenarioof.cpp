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
  char rayfile[100];
  if(argc==2)
    strcpy(rayfile,argv[1]);
  else
    argsError(argv[0]);
  
  char line[1000];
  FILE* fr=fopen(rayfile,"r");
  fgets(line,1000,fr);
  SpiceDouble et,tmp;

  ////////////////////////////////////////////////////
  //PREPARE FILE
  ////////////////////////////////////////////////////
  int i,k;
  SpiceDouble object[6];
  char header[100];
  FILE* fs=fopen("scenario.dat","w");
  fprintf(fs,"%-20s","#1:t");
  k=2;
  for(i=0;i<NUMOBJS;i++){
      if(!ACTIVE[i]) continue;
      sprintf(header,"#%d:i",k++);
      fprintf(fs,"%-7s",header);
      sprintf(header,"#%d:x%d",k++,i);
      fprintf(fs,"%-26s",header);
      sprintf(header,"#%d:y%d",k++,i);
      fprintf(fs,"%-26s",header);
      sprintf(header,"#%d:z%d",k++,i);
      fprintf(fs,"%-26s",header);
      sprintf(header,"#%d:vx%d",k++,i);
      fprintf(fs,"%-26s",header);
      sprintf(header,"#%d:vy%d",k++,i);
      fprintf(fs,"%-26s",header);
      sprintf(header,"#%d:vz%d",k++,i);
      fprintf(fs,"%-26s",header);
  }
  fprintf(fs,"\n\n");  
  while(1){
    ////////////////////////////////////////////////////
    //GET TIME
    ////////////////////////////////////////////////////
    fscanf(fr,"%lf",&et);
    if(feof(fr)) break;
    for(i=2;i<=17;i++) fscanf(fr,"%lf",&tmp);
    printf("T = %.9e\n",et);
    
    ////////////////////////////////////////////////////
    //COMPUTE EPHEMERIS
    ////////////////////////////////////////////////////
    fprintf(fs,"%-+20.9e",et);
    for(i=0;i<NUMOBJS;i++){
      if(!ACTIVE[i]) continue;
      spkezr_c(OBJS[i],et,ECJ2000,"NONE",SSB,object,&tmp);
      fprintf(fs,"%-7d%s",i,vec2strn(object,6,"%-+26.17e"));
    }
    fprintf(fs,"\n\n");
  }
  fclose(fs);

  ////////////////////////////////////////////////////
  //PLAIN OUTPUT
  ////////////////////////////////////////////////////
  fprintf(stdout,"--PLAIN--\n");
}
