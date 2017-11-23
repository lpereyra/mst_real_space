#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "variables.h"
#include "cosmoparam.h"
#include "leesloanreal.h"
#include "timer.h"
#include "voronoi.h"
#include "mst_kruskal.h"
#include "colores.h"

#ifdef GAL_LUM
  std::vector<std::pair<float,int> > orden;
  #ifdef BRANCH_SURVIVE
    float mr_survive;
  #endif
#else
  std::vector<std::pair<int,int> > orden;
  #ifdef BRANCH_SURVIVE
    int N_part_survive;
  #endif
#endif

void Write_Segments(int *Padre, int *Rank, double *fof);

int main(int argc, char **argv)
{
  int    i;
  int    *Padre, *Rank;
  double start,end;
  std::vector<std::vector<int> > adjacency_list;
  std::vector<std::pair<float,std::pair<int,int> > > edges;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  #ifdef BRANCH_SURVIVE
  BLUE("********** Important *************\n");
    #ifdef GAL_LUM
     mr_survive = atof(argv[2]);
     if(mr_survive>=0) 
        RED("******  WARNING ******\n");
      sprintf(message,"Cut Mr Mag %g\n",mr_survive);
    #else
      N_part_survive = atoi(argv[2]);
      if(N_part_survive==0) 
        RED("******  WARNING ******\n");
      sprintf(message,"Cut Numpar %d\n",N_part_survive);
    #endif
  RED(message);
  BLUE("**********************************\n");
  fflush(stdout);
  #endif

  readoutsloanreal();
 
  #ifndef GAL_LUM
    read_grup_fof(fof);
  #endif 

  /// Importante sino los grupos quedan fuera del box ///
  change_positions(cp.npart);

  Voronoi_Grupos(fof[0],edges);

  fprintf(stdout,"%d NumEdges\n",(int)edges.size());
  fflush(stdout);

  Padre = (int *) malloc(ngrup*sizeof(int));
  Rank =  (int *) malloc(ngrup*sizeof(int));

  for(i=0;i<ngrup;i++)
  {
    Padre[i] = i;
    Rank[i] = 0;
    adjacency_list.push_back(std::vector<int>());
  }

  Kruskal(Padre,Rank,edges,adjacency_list);

  Podado(4,adjacency_list);
  
  #ifdef BRANCH_SURVIVE
    #ifdef GAL_LUM
      fprintf(stdout,"Sobreviven en ramas los nodos con magnitud %f\n",mr_survive);
    #else
      fprintf(stdout,"Sobreviven en ramas los nodos con %d\n",N_part_survive);
    #endif
  fflush(stdout);
  #endif

  for(i=0;i<ngrup;i++)
  {
    Padre[i] = -1;
    Rank[i]  = 0;

    if(adjacency_list[i].size()>0)
    #ifdef GAL_LUM
      orden.push_back(std::make_pair(fabs(Gr[i].mr),i));
    #else   
      orden.push_back(std::make_pair(Gr[i].NumPart,i));
    #endif
  }

  sort(orden.begin(),orden.end());
  
  for(i=(int)orden.size()-1;i>=0;i--)
  {
    int k = orden[i].second;

    if(Padre[k]==-1)
      DLU(k,Padre[k],adjacency_list,Padre,Rank);
  }

  fprintf(stdout,"Escribe\n");
  fflush(stdout);

  Write_Segments(Padre,Rank,fof);

  free(Gr);
  free(Padre);
  free(Rank);
  adjacency_list.clear();

  TIMER(end);
  fprintf(stdout,"Total time %f\n",end-start);
  fflush(stdout);

  return(EXIT_SUCCESS);
}

void Write_Segments(int *Padre, int *Rank, double *fof)
{
  char filename[200];
  char filename_ascii[200];
  FILE *pfout, *pfpropiedades;
  FILE *pfout_ascii, *pfpropiedades_ascii;
  int i,j,k,id;
  int *contador;
  float dx,dy,dz;
  float dux,duy,duz;
  float r,lenr,elong,rms;
  std::vector<int> aux;

  j = 0;
  #ifdef GAL_LUM

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_segmentos_%.2f_%2.2f_lum_%2.2f.bin",fabs(mr_survive),zcut,fof[0],fabs(mcut));
        sprintf(filename_ascii,"%.2f_segmentos_%.2f_%2.2f_lum_%2.2f.dat",fabs(mr_survive),zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"%.2f_segmentos_%2.2f_lum_%2.2f.bin",fabs(mr_survive),fof[0],fabs(mcut));
        sprintf(filename_ascii,"%.2f_segmentos_%2.2f_lum_%2.2f.dat",fabs(mr_survive),fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #else
      #ifdef LIM_VOL
        sprintf(filename,"segmentos_%.2f_%2.2f_lum_%2.2f.bin",zcut,fof[0],fabs(mcut));
        sprintf(filename_ascii,"segmentos_%.2f_%2.2f_lum_%2.2f.dat",zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"segmentos_%2.2f_lum_%2.2f.bin",fof[0],fabs(mcut));
        sprintf(filename_ascii,"segmentos_%2.2f_lum_%2.2f.dat",fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #endif   // close BRANCH_SURVIVE
    pfout=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfout);
   
    pfout_ascii = fopen(filename_ascii,"w");

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_propiedades_%.2f_%2.2f_lum_%2.2f.bin",fabs(mr_survive),zcut,fof[0],fabs(mcut));
        sprintf(filename_ascii,"%.2f_propiedades_%.2f_%2.2f_lum_%2.2f.dat",fabs(mr_survive),zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"%.2f_propiedades_%2.2f_lum_%2.2f.bin",fabs(mr_survive),fof[0],fabs(mcut));
        sprintf(filename_ascii,"%.2f_propiedades_%2.2f_lum_%2.2f.dat",fabs(mr_survive),fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #else
      #ifdef LIM_VOL
        sprintf(filename,"propiedades_%.2f_%2.2f_lum_%2.2f.bin",zcut,fof[0],fabs(mcut));
        sprintf(filename_ascii,"propiedades_%.2f_%2.2f_lum_%2.2f.dat",zcut,fof[0],fabs(mcut));
      #else
        sprintf(filename,"propiedades_%2.2f_lum_%2.2f.bin",fof[0],fabs(mcut));
        sprintf(filename_ascii,"propiedades_%2.2f_lum_%2.2f.dat",fof[0],fabs(mcut));
      #endif // close LIM_VOL
    #endif   // close BRANCH_SURVIVE
    pfpropiedades=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfpropiedades);

    pfpropiedades_ascii=fopen(filename_ascii,"w");

  #else

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_segmentos_%.2f_%2.2f_%2.2f.bin",(float)N_part_survive,zcut,fof[0],fof[1]);
        sprintf(filename_ascii,"%.2f_segmentos_%.2f_%2.2f_%2.2f.dat",(float)N_part_survive,zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"%.2f_segmentos_%2.2f_%2.2f.bin",(float)N_part_survive,fof[0],fof[1]);
        sprintf(filename_ascii,"%.2f_segmentos_%2.2f_%2.2f.dat",(float)N_part_survive,fof[0],fof[1]);
      #endif // close LIM_VOL
    #else
      #ifdef LIM_VOL
        sprintf(filename,"segmentos_%.2f_%2.2f_%2.2f.bin",zcut,fof[0],fof[1]);
        sprintf(filename_ascii,"segmentos_%.2f_%2.2f_%2.2f.dat",zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"segmentos_%2.2f_%2.2f.bin",fof[0],fof[1]);
        sprintf(filename_ascii,"segmentos_%2.2f_%2.2f.dat",fof[0],fof[1]);
      #endif // close LIM_VOL
    #endif   // close BRANCH_SURVIVE
    pfout=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfout);

    pfout_ascii=fopen(filename_ascii,"w");

    #ifdef BRANCH_SURVIVE    
      #ifdef LIM_VOL
        sprintf(filename,"%.2f_propiedades_%.2f_%2.2f_%2.2f.bin",(float)N_part_survive,zcut,fof[0],fof[1]);
        sprintf(filename_ascii,"%.2f_propiedades_%.2f_%2.2f_%2.2f.dat",(float)N_part_survive,zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"%.2f_propiedades_%2.2f_%2.2f.bin",(float)N_part_survive,fof[0],fof[1]);
        sprintf(filename_ascii,"%.2f_propiedades_%2.2f_%2.2f.dat",(float)N_part_survive,fof[0],fof[1]);
      #endif
    #else
      #ifdef LIM_VOL
        sprintf(filename,"propiedades_%.2f_%2.2f_%2.2f.bin",zcut,fof[0],fof[1]);
        sprintf(filename_ascii,"propiedades_%.2f_%2.2f_%2.2f.dat",zcut,fof[0],fof[1]);
      #else
        sprintf(filename,"propiedades_%2.2f_%2.2f.bin",fof[0],fof[1]);
        sprintf(filename_ascii,"propiedades_%2.2f_%2.2f.dat",fof[0],fof[1]);
      #endif
    #endif   // close BRANCH_SURVIVE
    pfpropiedades=fopen(filename,"w");
    fwrite(&j,sizeof(int),1,pfpropiedades);

    pfpropiedades_ascii=fopen(filename_ascii,"w");

  #endif // close GAL_LUM

  contador = (int *) calloc(3,sizeof(int));

  while(!orden.empty())
  {

    i = orden.back().second;
 
    if(Rank[i]==1 && Padre[i]>=0)
    {
      aux.push_back(i);
      id = Padre[i];
      Rank[i] *= -1;

      while(id>=0)
      {
        if(Rank[id]<0)
        {
          aux.push_back(id);
          break;
        }

        #ifdef BRANCH_SURVIVE    
          #ifdef GAL_LUM
          if(Gr[id].mr<mr_survive && Padre[id]!=-1)
          #else
          if(Gr[id].NumPart>N_part_survive && Padre[id]!=-1)
          #endif
        #else
          if(Rank[id]>2)
        #endif
        {
          aux.push_back(id);

          #ifdef SORT_DERECHA
            #ifdef GAL_LUM
            if(Gr[aux[0]].mr<Gr[aux.back()].mr) // las magnitudes son negativas recordar 
              reverse(aux.begin(),aux.end());
            #else
            if(Gr[aux[0]].NumPart>Gr[aux.back()].NumPart)
              reverse(aux.begin(),aux.end());
            #endif
          #endif
          
          k = 0;
          #ifdef BRANCH_SURVIVE    

            #ifdef GAL_LUM
              if(Gr[aux[0]].mr<mr_survive) k++;
              if(Gr[aux.back()].mr<mr_survive) k++;
            #else
              if(Gr[aux[0]].NumPart>N_part_survive) k++;
              if(Gr[aux.back()].NumPart>N_part_survive) k++;
            #endif

          #else

            if(abs(Rank[aux.back()])>2) k++;
            if(abs(Rank[aux[0]])>2) k++;

          #endif

          contador[k]++;

          fwrite(&k,sizeof(int),1,pfpropiedades);
          fprintf(pfpropiedades_ascii,"%d ",k);

          k = (int)aux.size();
          fwrite(&k,sizeof(int),1,pfout);
          fprintf(pfout_ascii,"%d ",k);
          fwrite(&k,sizeof(int),1,pfpropiedades);
          fprintf(pfpropiedades_ascii,"%d ",k);

          dux = Gr[aux.back()].Pos[0] - Gr[aux[0]].Pos[0];
          duy = Gr[aux.back()].Pos[1] - Gr[aux[0]].Pos[1];
          duz = Gr[aux.back()].Pos[2] - Gr[aux[0]].Pos[2];

          #ifdef PERIODIC
          dux = dux >= cp.lbox*0.5 ? dux-cp.lbox : dux;
          dux = dux < -cp.lbox*0.5 ? dux+cp.lbox : dux;

          duy = duy >= cp.lbox*0.5 ? duy-cp.lbox : duy;
          duy = duy < -cp.lbox*0.5 ? duy+cp.lbox : duy;

          duz = duz >= cp.lbox*0.5 ? duz-cp.lbox : duz;
          duz = duz < -cp.lbox*0.5 ? duz+cp.lbox : duz;
          #endif

          elong = dux*dux+duy*duy+duz*duz;

          lenr = rms = 0.0;

          for(k=0;k<(int)aux.size();k++)
          {

          #ifdef GAL_LUM
            fwrite(&Gr[aux[k]].id,sizeof(int),1,pfout);
            fprintf(pfout_ascii,"%d ",Gr[aux[k]].id);
          #else
            fwrite(&aux[k],sizeof(int),1,pfout);
            fprintf(pfout_ascii,"%d ",aux[k]);
          #endif

            if(k==0) continue;

            dx = Gr[aux[k]].Pos[0] - Gr[aux[k-1]].Pos[0];
            dy = Gr[aux[k]].Pos[1] - Gr[aux[k-1]].Pos[1];
            dz = Gr[aux[k]].Pos[2] - Gr[aux[k-1]].Pos[2];

            #ifdef PERIODIC
            dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
            dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

            dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
            dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

            dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
            dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
            #endif
            
            r = sqrt(dx*dx+dy*dy+dz*dz);

            lenr += r;

            if(k==(int)aux.size()-1) continue;

            dx = Gr[aux[k]].Pos[0] - Gr[aux[0]].Pos[0];
            dy = Gr[aux[k]].Pos[1] - Gr[aux[0]].Pos[1];
            dz = Gr[aux[k]].Pos[2] - Gr[aux[0]].Pos[2];

            #ifdef PERIODIC
            dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
            dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

            dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
            dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

            dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
            dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
            #endif
                
            r = pow(dy*duz-dz*duy,2);
            r += pow(dz*dux-dx*duz,2);
            r += pow(dx*duy-dy*dux,2);
            r /= elong;
            
            rms+=r;

          }

          r = sqrt(elong)/lenr;

          k = (int)aux.size();
          rms /= (float)k;
          rms = sqrt(rms);

          fwrite(&lenr,sizeof(float),1,pfpropiedades);
          fwrite(&r,sizeof(float),1,pfpropiedades);
          fwrite(&rms,sizeof(float),1,pfpropiedades);

          fprintf(pfpropiedades_ascii,"%f %f %f\n",lenr,r,rms);
          fprintf(pfout_ascii,"\n");

          aux.clear();
          j++;
        }

        aux.push_back(id);
        Rank[id] *= -1;
        id = Padre[id];               
      }

      #ifdef SORT_DERECHA
        #ifdef GAL_LUM
          if(Gr[aux[0]].mr<Gr[aux.back()].mr) // las magnitudes son negativas recordar 
            reverse(aux.begin(),aux.end());
        #else
          if(Gr[aux[0]].NumPart>Gr[aux.back()].NumPart)
            reverse(aux.begin(),aux.end());
        #endif
      #endif

      id = aux.back();

      k = 0;
      #ifdef BRANCH_SURVIVE    
        #ifdef GAL_LUM
          if(Gr[aux[0]].mr<mr_survive) k++;
          if(Gr[id].mr<mr_survive)     k++;
        #else
          if(Gr[aux[0]].NumPart>N_part_survive) k++;
          if(Gr[id].NumPart>N_part_survive)     k++;
        #endif
      #else
        if(abs(Rank[id])>2) k++;
        if(abs(Rank[aux[0]])>2) k++;
      #endif

      contador[k]++;
      fwrite(&k,sizeof(int),1,pfpropiedades);
      fprintf(pfpropiedades_ascii,"%d ",k);

      k = (int)aux.size();
      fwrite(&k,sizeof(int),1,pfout);
      fprintf(pfout_ascii,"%d ",k);
      fwrite(&k,sizeof(int),1,pfpropiedades);
      fprintf(pfpropiedades_ascii,"%d ",k);

      dux = Gr[id].Pos[0] - Gr[aux[0]].Pos[0];
      duy = Gr[id].Pos[1] - Gr[aux[0]].Pos[1];
      duz = Gr[id].Pos[2] - Gr[aux[0]].Pos[2];

      #ifdef PERIODIC
      dux = dux >= cp.lbox*0.5 ? dux-cp.lbox : dux;
      dux = dux < -cp.lbox*0.5 ? dux+cp.lbox : dux;
    
      duy = duy >= cp.lbox*0.5 ? duy-cp.lbox : duy;
      duy = duy < -cp.lbox*0.5 ? duy+cp.lbox : duy;
    
      duz = duz >= cp.lbox*0.5 ? duz-cp.lbox : duz;
      duz = duz < -cp.lbox*0.5 ? duz+cp.lbox : duz;
      #endif

      elong = dux*dux+duy*duy+duz*duz;

      lenr = rms = 0.0;

      for(k=0;k<(int)aux.size();k++)
      {

      #ifdef GAL_LUM
        fwrite(&Gr[aux[k]].id,sizeof(int),1,pfout);
        fprintf(pfout_ascii,"%d ",Gr[aux[k]].id);
      #else
        fwrite(&aux[k],sizeof(int),1,pfout);
        fprintf(pfout_ascii,"%d ",aux[k]);
      #endif

        if(k==0) continue;

        dx = Gr[aux[k]].Pos[0] - Gr[aux[k-1]].Pos[0];
        dy = Gr[aux[k]].Pos[1] - Gr[aux[k-1]].Pos[1];
        dz = Gr[aux[k]].Pos[2] - Gr[aux[k-1]].Pos[2];

        #ifdef PERIODIC
        dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

        dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

        dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif
        
        r = sqrt(dx*dx+dy*dy+dz*dz);
        
        lenr += r;

        if(k==(int)aux.size()-1) continue;

        dx = Gr[aux[k]].Pos[0] - Gr[aux[0]].Pos[0];
        dy = Gr[aux[k]].Pos[1] - Gr[aux[0]].Pos[1];
        dz = Gr[aux[k]].Pos[2] - Gr[aux[0]].Pos[2];

        #ifdef PERIODIC
        dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

        dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

        dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif
            
        r = pow(dy*duz-dz*duy,2);
        r += pow(dz*dux-dx*duz,2);
        r += pow(dx*duy-dy*dux,2);
        r /= elong;

        rms += r;

      }

      r = sqrt(elong)/lenr;

      k = (int)aux.size();
      rms /= (float)k;
      rms = sqrt(rms);

      fwrite(&lenr,sizeof(float),1,pfpropiedades);
      fwrite(&r,sizeof(float),1,pfpropiedades);
      fwrite(&rms,sizeof(float),1,pfpropiedades);

      fprintf(pfpropiedades_ascii,"%f %f %f\n",lenr,r,rms);
      fprintf(pfout_ascii,"\n");

      aux.clear();
      j++;
    }

    orden.pop_back();
  }

  rewind(pfout);
  fwrite(&j,sizeof(int),1,pfout);
  fclose(pfout);
  fclose(pfout_ascii);

  rewind(pfpropiedades);
  fwrite(&j,sizeof(int),1,pfpropiedades);
  fclose(pfpropiedades);
  fclose(pfpropiedades_ascii);

  for(i=0;i<3;i++)
    fprintf(stdout,"Segmentos %d %d\n",i,contador[i]);
  fprintf(stdout,"Segmentos %d\n",j);
  free(contador);

  return;

}
