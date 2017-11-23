#include "variables.h"
#include "cosmoparam.h"
#include "colores.h"
#include "leesloanreal.h"

int    ngrup;
struct particle_data *P;
struct grup_data *Gr;
struct cosmoparam cp;
struct galsloanrealg *gal;

void readoutsloanreal()
{

  FILE *pf;
  char filename[200];
  int i,j,k;
  int ngalssloan;
  float disred;
  double dismin = 1.E26, dismax = -1.E26;
  double alfamin = 1.E26, alfamax = -1.E26;
  double deltamin = 1.E26, deltamax = -1.E26;
  #ifdef GAL_LUM
  float mtest;
  ngrup = 0;
  #endif

  cp.omegam = 0.3                       ;  /* OMEGA MATERIA                              */
  cp.omegal = 0.7                       ;  /* OMEGA LAMBDA                               */
  cp.omegak = 1.0-cp.omegam-cp.omegal   ;  /* OMEGA CURVATURA                            */
  cp.h0     = 100.                      ;  /* ESTO DEJA TODO EN UNIDADES DE H^-1         */

  RED("Read Out Sloan...\n");

  sprintf(filename,"%s%s",snap.root,snap.name);

  pf = fopen(filename,"r");

  if(pf == NULL)
  {
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  ngalssloan = cp.npart;
  fread(&cp.npart,sizeof(int),1,pf);

  if(cp.npart!=ngalssloan){ 
    fprintf(stdout,"se rompe %d %d\n",cp.npart,ngalssloan);
    exit(0);
    }

  /* ALOCATACION Y LECTURA */
  gal = (struct galsloanrealg *) malloc(ngalssloan*sizeof(struct galsloanrealg));
  P   = (struct particle_data *) malloc(ngalssloan*sizeof(struct particle_data));
  #ifdef GAL_LUM
  Gr = (struct grup_data *) malloc(cp.npart*sizeof(struct grup_data));
  #endif

  for(i=0;i<3;i++)
  {
    pmin[i] = 1.E26; 
    pmax[i] = -1.E26;
  }

  j = k = 0;
  for(i=0;i<ngalssloan;i++)
  {
    fread(&gal[i].gal.id,sizeof(int),1,pf);
    fread(&gal[i].gal.id_parent,sizeof(int),1,pf);
    fread(&gal[i].gal.id_NYU_VAGC,sizeof(int),1,pf);
    fread(&gal[i].gal.alfa,sizeof(float),1,pf);
    fread(&gal[i].gal.delta,sizeof(float),1,pf);
    fread(&gal[i].gal.red,sizeof(float),1,pf);
    fread(&gal[i].gal.mr,sizeof(float),1,pf);
    fread(&gal[i].gal.mr_lim,sizeof(float),1,pf);
    fread(&gal[i].gal.completenness,sizeof(float),1,pf);
    fread(&gal[i].gal.KE_petro[0],sizeof(float),1,pf);
    fread(&gal[i].gal.KE_petro[1],sizeof(float),1,pf);
    fread(&gal[i].gal.KE_model[0],sizeof(float),1,pf);
    fread(&gal[i].gal.KE_model[1],sizeof(float),1,pf);
    fread(&gal[i].gal.id_red_source_type,sizeof(int),1,pf);
    fread(&gal[i].gal.red_FOG,sizeof(float),1,pf);
    fread(&gal[i].gal.red_kaiser,sizeof(float),1,pf);
    fread(&gal[i].gal.red_real,sizeof(float),1,pf);
    fread(&gal[i].grupo[0],sizeof(int),1,pf);  
    fread(&gal[i].grupo[1],sizeof(int),1,pf);  
    fread(&gal[i].npgrupo[0],sizeof(int),1,pf);  
    fread(&gal[i].npgrupo[1],sizeof(int),1,pf);  

    gal[i].gal.alfa  *= PI180;
    gal[i].gal.delta *= PI180;

    disred = red2dis(gal[i].gal.red_real); /*EN MPC*/

    P[i].Pos[0] = disred*cos(gal[i].gal.delta)*cos(gal[i].gal.alfa) ;
    P[i].Pos[1] = disred*cos(gal[i].gal.delta)*sin(gal[i].gal.alfa) ;
    P[i].Pos[2] = disred*sin(gal[i].gal.delta)                  ;
    P[i].Dis    = disred;
    P[i].sub    = gal[i].grupo[0];

    if(P[i].Pos[0] > pmax[0]) pmax[0] = P[i].Pos[0];
    if(P[i].Pos[0] < pmin[0]) pmin[0] = P[i].Pos[0];
    if(P[i].Pos[1] > pmax[1]) pmax[1] = P[i].Pos[1];
    if(P[i].Pos[1] < pmin[1]) pmin[1] = P[i].Pos[1];
    if(P[i].Pos[2] > pmax[2]) pmax[2] = P[i].Pos[2];
    if(P[i].Pos[2] < pmin[2]) pmin[2] = P[i].Pos[2];

    if(P[i].Dis > dismax) dismax = P[i].Dis;
    if(P[i].Dis < dismin) dismin = P[i].Dis;

    if(gal[i].gal.alfa > alfamax) alfamax = gal[i].gal.alfa;
    if(gal[i].gal.alfa < alfamin) alfamin = gal[i].gal.alfa;
    if(gal[i].gal.delta > deltamax) deltamax = gal[i].gal.delta;
    if(gal[i].gal.delta < deltamin) deltamin = gal[i].gal.delta;

    if(gal[i].grupo[0] > j) j = gal[i].grupo[0];
    if(gal[i].grupo[1] > k) k = gal[i].grupo[1];

    #ifdef GAL_LUM

    if(P[i].sub==0) continue;

    mtest = gal[i].gal.mr-25.0-5.0*log10(disred*(1.0+gal[i].gal.red_real));

    if(mtest<mcut)
    {
      Gr[ngrup].save = P[i].sub;
      Gr[ngrup].id = i;
      memcpy(Gr[ngrup].Pos,P[i].Pos,3*sizeof(float));
      Gr[ngrup].mr = mtest;
      ngrup++;
    }

    #endif
  }

  cp.vol = (alfamax-alfamin)*(cos(deltamin)-cos(deltamax))*(pow(dismax,3.)-pow(dismin,3.))/3.0;
  cp.vol = fabs(cp.vol);

  #ifdef GAL_LUM
  fprintf(stdout,"%d Num de glx mas lum que %f\n",ngrup,mcut);
  fflush(stdout);
  Gr = (struct grup_data *) realloc (Gr,ngrup*sizeof(struct grup_data));
  #endif

  fprintf(stdout,"DisMin %.8e DisMax %.8e\n",dismin,dismax);
  fprintf(stdout,"AlfaMin %.8e AlfaMax %.8e\n",alfamin,alfamax);
  fprintf(stdout,"DeltaMin %.8e DeltaMax %.8e\n",deltamin,deltamax);
  fprintf(stdout,"Volumen aprox %.8e\n",cp.vol);
  fprintf(stdout,"Num Total %d\n",cp.npart);
  fprintf(stdout,"Num Total de grupos en la primera identificacion %d\n",j);
  fprintf(stdout,"Num Total de grupos en la segunda identificacion %d\n",k);
  fflush(stdout);

  #ifdef GAL_LUM
  fprintf(stdout,"%d Num de glx mas lum que %f\n",ngrup,mcut);
  fflush(stdout);
  #endif

  RED("End Read Sloan\n");

  fclose(pf);

}

#ifndef GAL_LUM

void read_grup_fof(double *fof)
{
  char filename[200];
  int  i;
  FILE *pfin;
 
  GREEN("Read Grups...\n");
  sprintf(filename,"%s%.2f_%.2f_centros.bin",snap.root,fof[1],fof[0]);
  pfin = fopen(filename,"rb"); 

  fread(&ngrup,sizeof(int),1,pfin);

  fprintf(stdout,"Grupos %d\n",ngrup);
  fflush(stdout);

  Gr = (struct grup_data *) malloc(ngrup*sizeof(struct grup_data));

  for(i=0;i<ngrup;i++)
  {
    fread(&Gr[i].save,sizeof(int),1,pfin);
    fread(&Gr[i].id,sizeof(int),1,pfin);
    fread(&Gr[i].Pos[0],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[1],sizeof(float),1,pfin);
    fread(&Gr[i].Pos[2],sizeof(float),1,pfin);
    fread(&Gr[i].NumPart,sizeof(int),1,pfin);
  }

  fclose(pfin);

  GREEN("End Read Grups\n");

}

#endif

double funlum(double x, void *p)
{
  struct paramfl *par=(struct paramfl *)p ;
  double ma=(par->ma)                     ;
  double alfa=(par->alfa)                 ;
  double t                                ;

  t=pow(10.0,0.4*(ma-x)) ;
  t=pow(t,1.0+alfa)*exp(-t) ;

  return(t) ;
}

double intfl(double x1, double x2)
{
  double resultado, error;
  size_t neval;
  struct paramfl pfl;
  gsl_function F; 

  pfl.ma=flma;
  pfl.alfa=flalfa;
  pfl.fia=flfia;

  F.function = &funlum;
  F.params = &pfl;

  if(x1<MAGMENOSINF)x1=MAGMENOSINF;
  gsl_integration_qng(&F,x1,x2,1.0e-7,1.0e-7,&resultado,&error,&neval);
  return(resultado);
}

double f(double z, void *p) /* FUNCION A INTEGRAR PARA LA DISTANCIA EN FUNCION DE Z*/
{
  struct paramcos *par=(struct paramcos *)p ;
  double om=(par->omegam);
  double ol=(par->omegal);
  double ok=(par->omegak);
  double q;
  q=pow(1.+z,3.)*om+pow(1.+z,2.)*ok+ol;
  return(1.0/sqrt(q));
}

double red2dis(double z)
{
  double resultado, error;
  size_t neval;
  struct paramcos pcos;
  gsl_function F;

  pcos.omegam=cp.omegam;
  pcos.omegal=cp.omegal;
  pcos.omegak=cp.omegak;

  F.function=&f;
  F.params=&pcos;

  gsl_integration_qng(&F,0.0,z,1.0e-7,1.0e-7,&resultado,&error,&neval);
 
  return(CVEL/cp.h0*resultado);
}

void change_positions(int n)
{
  int ip, idim;
  RED("Inicio Change Positions\n");

  printf("xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  printf("ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  printf("zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);

  for(ip = 0; ip < n; ip++)
    for(idim = 0; idim < 3; idim++)
      P[ip].Pos[idim] -= pmin[idim];

  for(ip = 0; ip < ngrup; ip++)
    for(idim = 0; idim < 3; idim++)
      Gr[ip].Pos[idim] -= pmin[idim];

  cp.lbox = pmax[0] - pmin[0];
  for(idim = 1; idim < 3; idim++)
    if(cp.lbox < (pmax[idim] - pmin[idim])) cp.lbox = (pmax[idim] - pmin[idim]);

  cp.lbox *= 1.001;
  fprintf(stdout,"Changing cp.lbox %f....\n",cp.lbox);
  GREEN("Fin Change Positions\n");
}
