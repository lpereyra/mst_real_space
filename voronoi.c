#include <assert.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <functional>
#include <algorithm>
#include <vector>

#include "cosmoparam.h"
#include "variables.h"
#include "grid.h"
#include "voronoi.h"
#include "leesloanreal.h" /// PARA EL RED2DIS
#include "voro++.hh"

bool dist_segment(double fof_2, int idg, float x, float y, float z)
{
  int i;
  int ixc, iyc, izc;
  int ixx, iyy, izz;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int xi, yi, zi;
  double fac;
  int ibox;
  double dx,dy,dz,r;

  fac = (double)grid.ngrid/(double)cp.lbox;

  ixc  = (int)((double)x*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)((double)y*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)((double)z*fac);
  izci = izc - 1;
  izcf = izc + 1;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= (int)grid.ngrid ) ixcf = (int)grid.ngrid - 1;
  if( iycf >= (int)grid.ngrid ) iycf = (int)grid.ngrid - 1;
  if( izcf >= (int)grid.ngrid ) izcf = (int)grid.ngrid - 1;
  #endif

  for(ixx = ixci; ixx <= ixcf; ixx++)
  {
    xi = ixx;
    #ifdef PERIODIC
    if(xi >= ngrid) xi -= grid.ngrid;
    if(xi < 0)      xi += grid.ngrid;
    #endif

    for( iyy = iyci ; iyy <= iycf ; iyy++)
    {

      yi = iyy;
      #ifdef PERIODIC
      if(yi >= ngrid) yi -= grind.ngrid;
      if(yi < 0)      yi += grind.ngrid;
      #endif

      for( izz = izci ; izz <= izcf ; izz++)
      {
        zi = izz;
        #ifdef PERIODIC
        if(zi >= ngrid) zi -= grid.ngrid;
        if(zi < 0)      zi += grid.ngrid;
        #endif

        ibox = (xi * grid.ngrid + yi) * grid.ngrid + zi;

        i = grid.llirst[ibox];

        while(i != -1)
        {

          if(P[i].sub==idg)
          {
            dx = P[i].Pos[0] - x;
            dy = P[i].Pos[1] - y;
            dz = P[i].Pos[2] - z;

            #ifdef PERIODIC
            dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
            dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

            dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
            dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

            dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
            dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
            #endif

            r = dx*dx+dy*dy+dz*dz;

            if(r<fof_2)
              return true;
          }

          i = grid.ll[i];

       }
      }
    }
  }

  return false;
}

void Voronoi_Grupos(double fof, std::vector<std::pair<float,std::pair<int,int> > > &edges)
{

  int i, j, idv, id, Ngrid, Tid, itera, N_threads;
  bool xbool, ybool, zbool;
  float  xp,yp,zp;
  double x_min, y_min, z_min;	  
  double x_max, y_max, z_max;	  
  double dx,dy,dz,r,r0,r0_2,frac;
  double rmablim,rmabmin;
  #ifdef LEN_MANUEL
  double rintlim;
  #endif
  std::vector<int>  vec;
  voro::voronoicell_neighbor cell;

  Ngrid = (int)pow((float)ngrup/5.0,1.0/3.0);
  #ifdef PERIODIC
    xbool = ybool = zbool = true;
  #else
    xbool = ybool = zbool = false;
  #endif
  r0 = 1.0e-10; // uso un pequeño gap
  x_min = y_min = z_min = 0.0-r0;	  
  x_max = y_max = z_max = cp.lbox+r0;	  

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  N_threads = NTHREADS > 16 ? 16 : NTHREADS;

  rmablim = rmaplim-25.0-5.0*log10(red2dis(0.1)*(1.+0.1));  // MAGNITUD ABSOLUTA LIMITE DEL CATALOGO  
  rmabmin = rmapmin-25.0-5.0*log10(red2dis(0.1)*(1.+0.1));  // MAGNITUD ABSOLUTA MINIMA DEL CATALOGO  
  #ifdef LEN_MANUEL
  rintlim = intfl(rmabmin,rmablim);                         // INTEGRAL ENTRE LAS MAGNITUDES LIMITES
  #endif

  fprintf(stdout, "INTEGRAL LIMITES\n");
  fprintf(stdout, "%f MABS correspondiente a %f MAPA\n",rmabmin, rmapmin);
  fprintf(stdout, "%f MABS correspondiente a %f MAPA\n",rmablim, rmaplim);
  fflush(stdout);
 
  #ifdef LEN_MANUEL
    r0 = 3.0/(4.0*M_PI*(fof+1)*(0.4*log(10.0)*flfia*rintlim));
    r0 = cbrt(r0);
  #else
    r0 = cbrt(cp.vol/cp.npart)*cbrt(1./(fof+1));
  #endif

  fprintf(stdout,"r0 %f\n",r0);
  fflush(stdout);

  /// INICIA LOS PARAMETROS DE LA GRILLA ///
  grid.nobj = cp.npart;
  grid.ngrid = (int)(1.1*cp.lbox/r0);
  grid.step = 1;

  fprintf(stdout,"r0 %f lbox %f\n",r0,cp.lbox);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }
  //////////////////////////////////////////
  grid_init();
  grid_build();

  r0_2 = r0*r0; // Para usar el cuadrado

  // Creo el container
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,Ngrid,Ngrid,Ngrid,xbool,ybool,zbool,8);

  for(i=0;i<ngrup; i++)
    con.put(i,Gr[i].Pos[0],Gr[i].Pos[1],Gr[i].Pos[2]);	  

  assert(ngrup==con.total_particles());

  voro::c_loop_all clo(con);
  std::vector<std::vector<std::pair<float,std::pair<int,int> > > > lados(N_threads);

  if(clo.start()) do if(con.compute_cell(cell,clo))
  {

    id = clo.pid();
    cell.neighbors(vec);

    #pragma omp parallel for num_threads(N_threads) schedule(static) default(none) \
    private(Tid,j,idv,dx,dy,dz,r,itera,frac,xp,yp,zp) shared(id,r0,r0_2,vec,cp,Gr,grid,lados,edges,stdout)
    for(j=0; j<(int)vec.size(); j++)
    {

      Tid = omp_get_thread_num();
      idv = vec[j];

      #ifndef PERIODIC
      if(idv<0) continue;
      #endif

      if(Gr[id].save!=Gr[idv].save) continue;

      if(Gr[id].id>Gr[idv].id)
      {

        dx = Gr[idv].Pos[0] - Gr[id].Pos[0];
        dy = Gr[idv].Pos[1] - Gr[id].Pos[1];
        dz = Gr[idv].Pos[2] - Gr[id].Pos[2];

        #ifdef PERIODIC
        dx = dx >= cp.lbox*0.5 ? dx-cp.lbox : dx;
        dx = dx < -cp.lbox*0.5 ? dx+cp.lbox : dx;

        dy = dy >= cp.lbox*0.5 ? dy-cp.lbox : dy;
        dy = dy < -cp.lbox*0.5 ? dy+cp.lbox : dy;

        dz = dz >= cp.lbox*0.5 ? dz-cp.lbox : dz;
        dz = dz < -cp.lbox*0.5 ? dz+cp.lbox : dz;
        #endif

        r = sqrt(dx*dx+dy*dy+dz*dz);

        if(r0<r)
        {

          itera = (int)(r/r0);

          while(itera>0)
          {
            frac = (double)itera*(r0/r);

            xp = Gr[id].Pos[0]+frac*dx;
            yp = Gr[id].Pos[1]+frac*dy;
            zp = Gr[id].Pos[2]+frac*dz;

            #ifdef PERIODIC
            xp = xp >= cp.lbox ? xp-cp.lbox : xp;
            xp = xp <  0.0     ? xp+cp.lbox : xp;

            yp = yp >= cp.lbox ? yp-cp.lbox : yp;
            yp = yp <  0.0     ? yp+cp.lbox : yp;

            zp = zp >= cp.lbox ? zp-cp.lbox : zp;
            zp = zp <  0.0     ? zp+cp.lbox : zp;
            #endif

            if(!dist_segment(r0_2,Gr[id].save,xp,yp,zp))
              break;
            else  
              itera--; 
          }

          if(itera!=0) continue;

        }

        #ifdef GAL_LUM
        // MAGsun_r = 4.71         // mag_r Hill et al 2010 - https://arxiv.org/pdf/1002.3788.pdf
          r = -(pow(10.0,-0.4*(Gr[id].mr-4.71))*pow(10.0,-0.4*(Gr[idv].mr-4.71)))/(r*r);

        #else
          r = -((float)Gr[id].NumPart*(float)Gr[idv].NumPart)/r;
        #endif

        lados[Tid].push_back(std::make_pair((float)r,std::make_pair(idv,id)));

      }
    }

    vec.clear();

   }while(clo.inc());

  con.clear();
  grid_free();
  free(P);

  for(i=0;i<N_threads;i++)
  {
    edges.insert(edges.end(),lados[i].begin(),lados[i].end());
    lados[i].clear();
  }

  return;
}


