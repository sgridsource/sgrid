/* Coordinates_set_bfaces.c */
/* Wolfgang Tichy 5/2017 */

#include "sgrid.h"
#include "Coordinates.h"


/* find box size L of smallest box */
double smallest_box_size(tGrid *grid)
{
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int i;
  double dx,dy,dz;
  double L=1e300;
  forallboxes(grid, i)
  {
    tBox *boxi = grid->box[i];
    double *pX = boxi->v[var_X];
    double *pY = boxi->v[var_Y];
    double *pZ = boxi->v[var_Z];
    double *px = boxi->v[var_x];
    double *py = boxi->v[var_y];
    double *pz = boxi->v[var_z];
    int n1 = boxi->n1;
    int n2 = boxi->n2;
    int n3 = boxi->n3;
    double Li;

    i = Index(n1,n2,n3);
    if(px!=NULL)
    {
      dx = px[i]-px[0];
      dy = py[i]-py[0];
      dz = pz[i]-pz[0];
    }
    else
    {
      dx = pX[i]-pX[0];
      dy = pY[i]-pY[0];
      dz = pZ[i]-pZ[0];
    }
    Li = sqrt(dx*dx + dy*dy + dz*dz);
    if(Li<L) L=Li;
  }
  return L;
}


/* determine which box faces touch or overlap other boxes, returns extface */
/* extface[f]=1 means face f is external, i.e. needs BC */
void find_external_faces_of_box(tBox *box, int *extface)
{
  tGrid *grid = box->grid;
  int pr = Getv("Coordinates_verbose", "yes");
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int b = box->b;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int d, dir;     /* we use dir = d+1 */ 
  int i,j,k, p, f;
  char str[1000];
  int *oblist = calloc(grid->nboxes,sizeof(int)); /* list that contains other boxes*/
  int nob; /* number of other boxes */
  int ob, oi;
  double oX,oY,oZ;
  double L, dL;
  double eps = dequaleps*1e3;

  /* mark faces in periodic dirs with not external, i.e. extface[f]=0  */
  for(f=0; f<6; f++)
  {
    if(box->periodic[1+f/2])  extface[f]=0;
    else                      extface[f]=1;
  }

  /* do nothing else if this is a Cartesian box */
  if(box->x_of_X[1]==NULL) return;

  /* find box size L of smallest box */
  L = smallest_box_size(grid);
  dL = L*eps;

  /* make oblist that contains all boxes */
  nob=0;
  forallboxes(grid, i) if(i!=b) nob=addto_boxlist(i, &oblist, nob);

  /* go over directions */
  for(d=0; d<3; d++)
  {
    dir = d+1;

    /* go over faces, but not edges */
    for(p=0; p<n1; p+=n1-1)
    {
      int s; /* sign s=+-1 */

      /* pick sign s of normal and pick face index */
      if(p==0) { s=-1;  f=2*d;   }
      else     { s=+1;  f=2*d+1; }

      /* mark other box as non-existent by default */
      ob = -1;

      /* do nothing if not an external face */
      if(!extface[f]) continue;

      /* look for points in other boxes just outside this box */
      forinnerplaneN(dir, i,j,k, n1,n2,n3, p)
      {
        int ijk = Index(i,j,k);
        double X = box->v[var_X][ijk];
        double Y = box->v[var_Y][ijk];
        double Z = box->v[var_Z][ijk];
        double x = box->v[var_x][ijk];
        double y = box->v[var_y][ijk];
        double z = box->v[var_z][ijk];
        double ox,oy,oz, Nx,Ny,Nz, Nmag, dx,dy,dz, dist;
        int li, ret;

        /* pick one of X,Y,Z on boundary */
        if(dir==1) X = box->bbox[f];
        if(dir==2) Y = box->bbox[f];
        if(dir==3) Z = box->bbox[f];

        /* use normal vector to find point ox,oy,oz slightly outside box */
        Nx = box->dx_dX[1][dir](box, ijk, X,Y,Z);
        Ny = box->dx_dX[2][dir](box, ijk, X,Y,Z);
        Nz = box->dx_dX[3][dir](box, ijk, X,Y,Z);
        Nmag = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
        dx = s*Nx*dL;
        dy = s*Ny*dL;
        dz = s*Nz*dL;
        if(Nmag!=0.0) { dx=dx/Nmag;  dy=dy/Nmag;  dz=dz/Nmag; }
        ox = x + dx;
        oy = y + dy;
        oz = z + dz;
//printf("b%d dir%d  p%d, f%d, i,j,k=%d,%d,%d\n", b, dir, p, f, i,j,k);
//printf(" Nx,Ny,Nz=%g,%g,%g Nmag=%g\n", Nx,Ny,Nz, Nmag);
//printf(" x,y,z=%g,%g,%g ox,oy,oz=%g,%g,%g\n", x,y,z, ox,oy,oz);
//printf(" dx,dy,dz=%g,%g,%g\n", dx,dy,dz);
        /* find point in other boxes */
        for(li=0; li<nob; li++)
        {
          int bi = oblist[li];
          tBox *obox = grid->box[bi];
          dist = nearestXYZ_of_xyz(obox, &oi, &oX,&oY,&oZ, ox,oy,oz);
          dist = sqrt(dist);
          if(dist<0.5*L)
          {
            ret=XYZ_of_xyz(obox, &oX,&oY,&oZ, ox,oy,oz);
            if(ret>=0)
            {
              if(!(obox->periodic[1]))
                if(dless(oX,obox->bbox[0]) || dless(obox->bbox[1],oX)) continue;
              if(!(obox->periodic[2]))
                if(dless(oY,obox->bbox[2]) || dless(obox->bbox[3],oY)) continue;
              if(!(obox->periodic[3]))
                if(dless(oZ,obox->bbox[4]) || dless(obox->bbox[5],oZ)) continue;
              ob=obox->b;
            }
          }
        }
        /* if we find one point in another box this face is external */
        if(ob>=0) 
        { i=n1; j=n2; k=n3; } /* leave plane loop if face is external*/
      }
      /* if ob<0, we found no other box face f is not external */
      if(ob<0) extface[f]=0;
    }
  } /* end loop over directions */

  free(oblist);
}


/* find and set all bfaces on an external box face f */
/* The idea is to loop over the box faces and then move out using the Cartesian
   normal vectors. Then we check if we are in an other box. */
int set_bfaces_on_boxface(tBox *box, int f)
{
  tGrid *grid = box->grid;
  int pr = Getv("Coordinates_verbose", "yes");
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  double *pX = box->v[var_X];
  double *pY = box->v[var_Y];
  double *pZ = box->v[var_Z];
  double *px = box->v[var_x];
  double *py = box->v[var_y];
  double *pz = box->v[var_z];
  int b = box->b;
  int nbfaces;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int d, dir;     /* we use dir = d+1 */ 
  int i,j,k, p, fi, li;
  int *oblist = calloc(grid->nboxes,sizeof(int)); /* list that contains other boxes*/
  int nob; /* number of other boxes */
  int ob, oi;
  double oX,oY,oZ;
  double eps = dequaleps*1e3;
  double L, dL;
  int s = 2*(f%2) - 1; /* direction of normal vector s=+-1 */

  d = f/2;
  dir = d+1;

  /* find box size L of smallest box */
  L = smallest_box_size(grid);
  dL = L*eps;

  /* make oblist that contains all boxes except b,
     and add one bface for each of the other boxes */
  nob=0;
  forallboxes(grid, i) if(i!=b)
  {
    nob=addto_boxlist(i, &oblist, nob);
    /* add one empty bface for each other box */
    fi = add_empty_bface(box, f);
    box->bface[fi]->ob = i;
  }

  /* look for points in other boxes just outside this box */
  p = ( (n1-1)*(dir==1) + (n2-1)*(dir==2) + (n3-1)*(dir==3) )*(f%2);
  //forinnerplaneN(dir, i,j,k, n1,n2,n3, p)
  forplaneN(dir, i,j,k, n1,n2,n3, p)
  {
    int ijk = Index(i,j,k);
    double X = pX[ijk];
    double Y = pY[ijk];
    double Z = pZ[ijk];
    double x,y,z;
    double ox,oy,oz, Nx,Ny,Nz, Nmag, dx,dy,dz;

    /* pick one of X,Y,Z on boundary */
    if(dir==1) X = box->bbox[f];
    if(dir==2) Y = box->bbox[f];
    if(dir==3) Z = box->bbox[f];

    if(box->x_of_X[1]==NULL) /* this is a Cartesian box */
    {
      x=X;
      y=Y;
      z=Z;
      /* normal vector */ 
      Nx = (dir==1);
      Ny = (dir==2);
      Nz = (dir==3);
    }
    else
    {
      x=px[ijk];
      y=py[ijk];
      z=pz[ijk];
      /* normal vector */
      if(box->dx_dX[1][dir]==NULL) errorexit("we need box->dx_dX[1][dir]");
      Nx = box->dx_dX[1][dir](box, ijk, X,Y,Z);
      Ny = box->dx_dX[2][dir](box, ijk, X,Y,Z);
      Nz = box->dx_dX[3][dir](box, ijk, X,Y,Z);
    }
    /* use normal vector to find point ox,oy,oz slightly outside box */
    Nmag = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
    dx = s*Nx*dL;
    dy = s*Ny*dL;
    dz = s*Nz*dL;
    if(Nmag!=0.0) { dx=dx/Nmag;  dy=dy/Nmag;  dz=dz/Nmag; }
    ox = x + dx;
    oy = y + dy;
    oz = z + dz;

    /* mark other box as non-existent by default */
    ob = -1;

    /* find point in other boxes */
    for(li=0; li<nob; li++)
    {
      int bi = oblist[li];
      int ret;
      tBox *obox = grid->box[bi];
      double dist = nearestXYZ_of_xyz(obox, &oi, &oX,&oY,&oZ, ox,oy,oz);
      dist = sqrt(dist);
//printf("  bi=%d dist=%g ob=%d oi=%d oX,oY,oZ=%g,%g,%g\n",
//bi,dist, ob, oi, oX,oY,oZ);
      if(dist<0.5*L)
      {
        ret=XYZ_of_xyz(obox, &oX,&oY,&oZ, ox,oy,oz);
        if(ret>=0)
        {
          if(!(obox->periodic[1]))
            if(dless(oX,obox->bbox[0]) || dless(obox->bbox[1],oX)) continue;
          if(!(obox->periodic[2]))
            if(dless(oY,obox->bbox[2]) || dless(obox->bbox[3],oY)) continue;
          if(!(obox->periodic[3]))
            if(dless(oZ,obox->bbox[4]) || dless(obox->bbox[5],oZ)) continue;
          ob=obox->b;
        }
      }
    }
    /* if we find a point we add it to the bface with the correct ob and f */
    if(ob>=0)
      for(fi=0; fi<box->nbfaces; fi++)
        if(box->bface[fi]->ob == ob && box->bface[fi]->f == f)
        {
          add_point_to_bface_inbox(box, fi, ijk, f);
          break;
        }
  }

  /* look for empty bfaces and remove them */
  remove_bfaces_with_NULL_ftps(box);

  free(oblist);
  return nbfaces;
}



/* set bfaces for each box on the grid */
int Coordinates_set_bfaces(tGrid *grid)
{
  int pr = Getv("Coordinates_verbose", "yes");
  int b;

  if(pr) printf("Coordinates_set_bfaces:\n");

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int extface[6];  /* extface[f]=1  means face f is external, i.e. needs BC */
    int f;

    /* find all external faces of box */
    find_external_faces_of_box(box, extface);
    if(pr)
    {
      printf("external faces on box%d: ", b);
      for(f=0; f<6; f++) if(extface[f]) printf(" %d", f);
      printf("\n");
    }
    /* set the bfaces for each face we found */
    for(f=0; f<6; f++)
      if(extface[f])  set_bfaces_on_boxface(box, f);

    if(pr)
    {
      printbfaces(box);
    }
  }
exit(88);
  return 0;
}
