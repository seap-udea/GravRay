  //*
  printf("ux : %s, uy: %s, uz: %s\n",
	 vec2str(ux,"%.17f"),
	 vec2str(uy,"%.17f"),
	 vec2str(normal,"%.17f")
	 );
  //*/

int EoM(double t,double y[],double dydt[],void *params) 
{ 
  //COMPUTE THE CONTRIBUTION OF EVERY OBJECT
  int i;
  double r,object[6],R[3],Rmag,tmp,GM,fac;

  fac=UT*UT/(UL/1E3*UL/1E3*UL/1E3);
  dydt[CX]=y[CVX];
  dydt[CY]=y[CVY];
  dydt[CZ]=y[CVZ];
  dydt[CVX]=0.0;
  dydt[CVY]=0.0;
  dydt[CVZ]=0.0;

  for(i=NUMOBJS;i-->0;){
    if(!ACTIVE[i]) continue;
    spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
    vscl_c(1E3/UL,object,object);
    sumVec(R,1.0,y,-1.0,object,3);
    Rmag=vnorm_c(R);
    GM=GMASSES[i]*fac;
    sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);
  }
  return 0;
}

int (*EoM)(double,double*,double*,void*);
  //EoM=EoM_unfold;
  EoM=EoM_loop;

int (*EoM)(double,double*,double*,void*);

int EoM_unfold(double t,double y[],double dydt[],void *params) 
{ 
  //COMPUTE THE CONTRIBUTION OF EVERY OBJECT
  static int q=1;
  int i;
  double r,object[6],R[3],Rmag,tmp,GM,fac;
  if(q){
    fprintf(stdout,"Unfold\n");
    q=0;
  }

  fac=UT*UT/(UL/1E3*UL/1E3*UL/1E3);
  dydt[CX]=y[CVX];
  dydt[CY]=y[CVY];
  dydt[CZ]=y[CVZ];
  dydt[CVX]=0.0;
  dydt[CVY]=0.0;
  dydt[CVZ]=0.0;

  i=-1;

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  i++;
  spkezr_c(OBJS[i],t*UT,ABSJ2000,"NONE",SSB,object,&tmp);
  vscl_c(1E3/UL,object,object);
  sumVec(R,1.0,y,-1.0,object,3);
  Rmag=vnorm_c(R);
  GM=ACTIVE[i]*GMASSES[i]*fac;
  sumVec(dydt+3,1.0,dydt+3,-GM/(Rmag*Rmag*Rmag),R,3);

  return 0;
}

