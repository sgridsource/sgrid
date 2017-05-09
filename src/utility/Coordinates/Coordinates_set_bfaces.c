/* Coordinates_set_bfaces.c */
/* Wolfgang Tichy 5/2017 */

#include "sgrid.h"
#include "Coordinates.h"


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
  int periodic[3];
  int d, dir;     /* we use dir = d+1 */ 
  int i,j,k, p, f;
  char str[1000];
  int *oblist = calloc(grid->nboxes,sizeof(int)); /* list that contains other boxes*/
  int nob; /* number of other boxes */
  int ob, oi;
  double oX,oY,oZ;
  double eps = dequaleps*1e5;
  
  /* make oblist that contains all boxes */
  nob=0;
  forallboxes(grid, i) if(i!=b) nob=addto_boxlist(i, &oblist, nob);

  /* check for periodic dirs */
  for(d=0; d<3; d++)
  {
    snprintf(str, 999, "box%d_basis%d", b, d+1);
    if( Getv(str, "Fourier") )  periodic[d]=1;
    else                        periodic[d]=0;
  }
  /* mark faces in periodic dirs with not external, i.e. extface[f]=0  */
  for(f=0; f<6; f++)
  {
    if(periodic[f/2])  extface[f]=0;
    else               extface[f]=1;
  }

  /* go over directions */
  for(d=0; d<3; d++)
  {
    double dL = (box->bbox[2*d+1] - box->bbox[2*d])*eps;
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
        double ox,oy,oz, Nx,Ny,Nz, Nmag, dx,dy,dz;

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
        /* find point in other boxes */
        nearest_b_XYZ_of_xyz_inboxlist(grid, oblist,nob, 
                                       &ob, &oi, &oX,&oY,&oZ, ox,oy,oz);
printf("Nx,Ny,Nz=%g,%g,%g Nmag=%g\n", Nx,Ny,Nz, Nmag);
printf("b%d dir%d  p%d, f%d, i,j,k=%d,%d,%d\n", b, dir, p, f, i,j,k);
printf("b=%d x,y,z=%g,%g,%g ox,oy,oz=%g,%g,%g\n",
b, x,y,z, ox,oy,oz);
printf("ob=%d oi=%d oX,oY,oZ=%g,%g,%g\n", ob, oi, oX,oY,oZ);
        ob=b_XYZ_of_xyz_inboxlist(grid, oblist,nob, &oX,&oY,&oZ, ox,oy,oz);
printf("ob=%d oX,oY,oZ=%g,%g,%g\n", ob, oX,oY,oZ);

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
  int b = box->b;
  int nbfaces;
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int d, dir;     /* we use dir = d+1 */ 
  int i,j,k, p, fi;
  int *oblist = calloc(grid->nboxes,sizeof(int)); /* list that contains other boxes*/
  int nob; /* number of other boxes */
  int ob, oi;
  double oX,oY,oZ;
  double eps = dequaleps*1000.0;
  double dL;
  int s = 2*(f%2) - 1; /* direction of normal vector s=+-1 */

  d = f/2;
  dir = d+1;
  dL = (box->bbox[2*d+1] - box->bbox[2*d])*eps;

  /* mark other box as non-existent by default */
  ob = -1;

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
  forinnerplaneN(dir, i,j,k, n1,n2,n3, p)
  {
    int ijk = Index(i,j,k);
    double X = box->v[var_X][ijk];
    double Y = box->v[var_Y][ijk];
    double Z = box->v[var_Z][ijk];
    double x = box->v[var_x][ijk];
    double y = box->v[var_y][ijk];
    double z = box->v[var_z][ijk];
    double ox,oy,oz, Nx,Ny,Nz, Nmag, dx,dy,dz;

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

if(j==0 || j==n2-1) continue;

    /* find point in other boxes */
    nearest_b_XYZ_of_xyz_inboxlist(grid, oblist,nob, 
                                   &ob, &oi, &oX,&oY,&oZ, ox,oy,oz);
//printf("");

    ob=b_XYZ_of_xyz_inboxlist(grid, oblist,nob, &oX,&oY,&oZ, ox,oy,oz);
    /* if we find a point we add it to the bface with the correct ob */
    if(ob>=0)
      for(fi=0; fi<box->nbfaces; fi++)
        if(box->bface[fi]->ob == ob) 
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
    int extface[6]; /* extface[f]=1 means face f is external, i.e. needs BC */
    int f;

    /* find all external faces of box */
    find_external_faces_of_box(box, extface);
    if(pr)
    {
      printf("  external faces in box%d: ", b);
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

  return 0;
}
