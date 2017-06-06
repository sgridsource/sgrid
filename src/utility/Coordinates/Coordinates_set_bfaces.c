/* Coordinates_set_bfaces.c */
/* Wolfgang Tichy 5/2017 */

#include "sgrid.h"
#include "Coordinates.h"

#define EPS dequaleps*1e3

#define LESS(a, b) ((a) < (b) - EPS*L )


/* find box size L of box */
double find_box_size(tBox *box)
{
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int var_X = Ind("X");
  int var_Y = Ind("Y");
  int var_Z = Ind("Z");
  int i,j,n, ind[8];
  double dx,dy,dz;
  double sum, L;
  double *pX = box->v[var_X];
  double *pY = box->v[var_Y];
  double *pZ = box->v[var_Z];
  double *px = box->v[var_x];
  double *py = box->v[var_y];
  double *pz = box->v[var_z];
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int m1 = n1-1;
  int m2 = n2-1;
  int m3 = n3-1;

  ind[0] = Index( 0, 0, 0);
  ind[1] = Index(m1, 0, 0);
  ind[2] = Index( 0,m2, 0);
  ind[3] = Index(m1,m2, 0);
  ind[4] = Index( 0, 0,m3);
  ind[5] = Index(m1, 0,m3);
  ind[6] = Index( 0,m2,m3);
  ind[7] = Index(m1,m2,m3);

//for(i=0; i<8; i++)
//printf("i%d %4d x=%g y=%g z=%g  X=%g Y=%g Z=%g\n",
//i, ind[i], px[ind[i]], py[ind[i]], pz[ind[i]],
//pX[ind[i]], pY[ind[i]], pZ[ind[i]]);

  sum=0.0;
  n=0;
  for(i=1; i<8; i++)
  for(j=0; j<i; j++)
  {
    if(px!=NULL)
    {
      dx = px[ind[i]]-px[ind[j]];
      dy = py[ind[i]]-py[ind[j]];
      dz = pz[ind[i]]-pz[ind[j]];
    }
    else
    {
      dx = pX[ind[i]]-pX[ind[j]];
      dy = pY[ind[i]]-pY[ind[j]];
      dz = pZ[ind[i]]-pZ[ind[j]];
    }
    sum += dx*dx + dy*dy + dz*dz;
    n++;
  }

  L = sqrt(sum/n);
  return L;
}

/* find box size L of smallest box */
double smallest_box_size(tGrid *grid)
{
  int i;
  double L=1e300;
  forallboxes(grid, i)
  {
    tBox *box = grid->box[i];
    double Li = find_box_size(box);
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
  intList *obl = alloc_intList(); /* list that contains other boxes*/
  int ob, oi;
  double oX,oY,oZ;
  double L, dL;

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
  dL = L*EPS;

  /* make obl that contains all boxes */
  forallboxes(grid, i) if(i!=b) unionpush_intList(obl, i);

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
        Nx = box->dx_dX[1][dir](box, -1, X,Y,Z);
        Ny = box->dx_dX[2][dir](box, -1, X,Y,Z);
        Nz = box->dx_dX[3][dir](box, -1, X,Y,Z);
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

        /* find point in this box */
        dist = nearestinnerXYZ_of_xyz(box, &oi, &oX,&oY,&oZ, ox,oy,oz);
        dist = sqrt(dist);
        ret=XYZ_of_xyz(box, &oX,&oY,&oZ, ox,oy,oz);
        if(ret>=0)
        {
          if(!(box->periodic[1]))
            if(oX < box->bbox[0] || box->bbox[1] < oX) ret=-1;
          if(!(box->periodic[2]))
            if(oY < box->bbox[2] || box->bbox[3] < oY) ret=-1;
          if(!(box->periodic[3]))
            if(oZ < box->bbox[4] || box->bbox[5] < oZ) ret=-1;
        }
        dist=sqrt(x*x + y*y + z*z);
        if(ret>=0 && dist<1e60) /* point is in this box and x,y,z is not inf  */
        {
//printf("%g ", dist);
//printf("%d %d %d x,y,z=%g,%g,%g dx,dy,dz=%g,%g,%g\n", i,j,k, x,y,z, dx,dy,dz);
//printf("%d %d %d x,y,z=%g,%g,%g Nx,Ny,Nz=%g,%g,%g\n", i,j,k, x,y,z, Nx,Ny,Nz);
          extface[f]=0; /* mark face as not external */
          goto endplaneloop; /* break; does not work for nested loop */
        }

        /* find point in other boxes */
        for(li=0; li<obl->n; li++)
        {
          int bi = obl->e[li];
          tBox *obox = grid->box[bi];
          dist = nearestXYZ_of_xyz(obox, &oi, &oX,&oY,&oZ, ox,oy,oz);
          //dist = nearestinnerXYZ_of_xyz(obox, &oi, &oX,&oY,&oZ, ox,oy,oz);
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
      } endplaneloop:
      /* if ob<0, we found no other box face f is not external */
      if(ob<0) extface[f]=0;
    }
  } /* end loop over directions */

  free_intList(obl);
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
  int n1 = box->n1;
  int n2 = box->n2;
  int n3 = box->n3;
  int d, dir;     /* we use dir = d+1 */ 
  int i,j,k, p, fi, li;
  intList *obl = alloc_intList(); /* list that contains other boxes*/
  int ob, oi;
  double oX,oY,oZ;
  double L, dL;
  int s = 2*(f%2) - 1; /* direction of normal vector s=+-1 */

  d = f/2;
  dir = d+1;

  /* find box size L of smallest box */
  L = smallest_box_size(grid);
  dL = L*EPS;

  /* make obl that contains all boxes except b,
     and add one bface for each of the other boxes */
  forallboxes(grid, i) if(i!=b)
  {
    unionpush_intList(obl, i);
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
    int face[6];
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
      Nx = box->dx_dX[1][dir](box, -1, X,Y,Z);
      Ny = box->dx_dX[2][dir](box, -1, X,Y,Z);
      Nz = box->dx_dX[3][dir](box, -1, X,Y,Z);

      /* if we are also on another face add a bit of the normal as well */
      if(XYZ_on_face(box, face, X,Y,Z))
      {
        int ff;
        face[f]=0; /* remove current face */
        for(ff=0; ff<6; ff++)
          if(face[ff])
          {
            int dir = 1+ff/2;
            Nx += 0.1 * ( box->dx_dX[1][dir](box, -1, X,Y,Z) );
            Ny += 0.1 * ( box->dx_dX[2][dir](box, -1, X,Y,Z) );
            Nz += 0.1 * ( box->dx_dX[3][dir](box, -1, X,Y,Z) );
          }
      }
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
//if(i==7 && j==1 && k==1)
//printf("b%d %d %d %d x,y,z=%g,%g,%g dx,dy,dz=%g,%g,%g\n", b, i,j,k, x,y,z, dx,dy,dz);

    /* mark other box as non-existent by default */
    ob = -1;

    /* find point in other boxes */
    for(li=0; li<obl->n; li++)
    {
      int bi = obl->e[li];
      int ret;
      tBox *obox = grid->box[bi];
      double osize = find_box_size(obox);
      double dist = nearestinnerXYZ_of_xyz(obox, &oi, &oX,&oY,&oZ, ox,oy,oz);
      //double dist = guessXYZ_of_xyz(obox, &oi, &oX,&oY,&oZ, ox,oy,oz);
      dist = sqrt(dist);
//printf("o%dsize=%g ", obox->b, osize);
//printf("  bi=%d dist=%g ob=%d oi=%d oX,oY,oZ=%g,%g,%g\n",
//bi,dist, ob, oi, oX,oY,oZ);
      if(dist<osize)
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
          /*
          if(!(obox->periodic[1]))
            if(oX < obox->bbox[0] || obox->bbox[1] < oX) continue;
          if(!(obox->periodic[2]))
            if(oY < obox->bbox[2] || obox->bbox[3] < oY) continue;
          if(!(obox->periodic[3]))
            if(oZ < obox->bbox[4] || obox->bbox[5] < oZ) continue;
          */
          /*
          if(!(obox->periodic[1]))
            if(LESS(oX,obox->bbox[0]) || LESS(obox->bbox[1],oX)) continue;
          if(!(obox->periodic[2]))
            if(LESS(oY,obox->bbox[2]) || LESS(obox->bbox[3],oY)) continue;
          if(!(obox->periodic[3]))
            if(LESS(oZ,obox->bbox[4]) || LESS(obox->bbox[5],oZ)) continue;
          */
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
  remove_bfaces_with_NULL_fpts(box);

  free_intList(obl);
  return box->nbfaces;
}


/* figure out how each bface is connected to other bfaces */
int set_ofi_in_all_bfaces(tGrid *grid)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int fi, ob, ofi;

    for(fi=0; fi<box->nbfaces; fi++)
    {
      tBox *obox;
      ob = box->bface[fi]->ob;
      if(ob<0) continue; /* do nothing if there is no other box */
      obox = grid->box[ob];
      for(ofi=0; ofi<obox->nbfaces; ofi++)
      {
        /* check if other box ob refers to this box b */
        if(obox->bface[ofi]->ob == b)
        {
          /* set box->bface[fi]->ofi=ofi if nothing else is in there */
          if(box->bface[fi]->ofi == -1 || box->bface[fi]->ofi == ofi)
            box->bface[fi]->ofi=ofi;
          /* if one other was already there use bits */
          else if(box->bface[fi]->ofi >= 0)
          {
            int c = box->bface[fi]->ofi;
            unsigned int cbt = 1<<(c+1);
            int n = ofi;
            unsigned int nbt = 1<<(n+1);
            int r = cbt | nbt;
            box->bface[fi]->ofi=-r;
          }
          /* if several were there already use bits */
          else
          {
            unsigned int cbt = -box->bface[fi]->ofi;
            unsigned int nbt = 1<<(ofi+1);
            int r = cbt | nbt;
            box->bface[fi]->ofi=-r;
          }
        } /* end if */
      }/* end ofi loop */
    }
  } /* end box loop */
  return 0;
}


/* set bfaces for each box on the grid */
int Coordinates_set_bfaces(tGrid *grid)
{
  int pr = Getv("Coordinates_verbose", "yes");
  int maxits = Geti("Coordinates_newtMAXITS"); /* save par */
  int b;

  if(pr) printf("Coordinates_set_bfaces:\n");

  /* reduce iteration number in newton_linesrch_itsP
     and then loop over boxes */
  Seti("Coordinates_newtMAXITS", 75);
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int extface[6];  /* extface[f]=1  means face f is external, i.e. needs BC */
    int f;

    free_all_bfaces(box);
/*
tBox *obox=grid->box[5];
int oi;
double oX,oY,oZ, ox,oy,oz;
oX=0.08; oY=1e-10; oZ=0;
ox=25.839;
oy=0.0;
oz=0;//1e-7;
printf("L=%g\n", find_box_size(obox));
nearestinnerXYZ_of_xyz(obox, &oi, &oX,&oY,&oZ, ox,oy,oz);
printf("G ox,oy,oz=%g,%g,%g  oX,oY,oZ=%g,%g,%g\n",ox,oy,oz , oX,oY,oZ);
XYZ_of_xyz(obox, &oX,&oY,&oZ, ox,oy,oz);
printf("S ox,oy,oz=%g,%g,%g  oX,oY,oZ=%g,%g,%g\n",ox,oy,oz , oX,oY,oZ);
//exit(78);
*/
    //check_box_dx_dX(box, 0.1, 0.2, 0.3);

    /* When not Cartesian, check if dx/dX exists */
    if(box->x_of_X[1]!=NULL)
      if(box->dx_dX[1][1]==NULL)
        errorexiti("Coordinates_set_bfaces: dx_dX[1][1]==NULL in box%d", b);

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
  /* restore Coordinates_newtMAXITS */
  Seti("Coordinates_newtMAXITS", maxits);

  /* set ofi and bit fields for all bfaces */
  set_ofi_in_all_bfaces(grid);
  set_bits_in_all_bfaces(grid);

  if(pr) forallboxes(grid, b) printbfaces(grid->box[b]);

prPointList(grid->box[3]->bface[0]->fpts);
prPointList(grid->box[3]->bface[1]->fpts);
prPointList(grid->box[4]->bface[0]->fpts);
prPointList(grid->box[4]->bface[4]->fpts);
exit(77);
exit(88);
  return 0;
}


/* set bit fields in bfaces for each box on the grid */
int set_bits_in_all_bfaces(tGrid *grid)
{
  int iX = Ind("X");
  int iY = iX+1;
  int iZ = iX+2;
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    double *X = box->v[iX];
    double *Y = box->v[iY];
    double *Z = box->v[iZ];
    int fi;

    /* loop over bfaces */
    for(fi=0; fi<box->nbfaces; fi++)
    {
      tBface *bface = box->bface[fi];
      int ob  = bface->ob;
      int ofi = bface->ofi;
      int oXi = bface->oXi;
      int oYi = bface->oYi;
      int oZi = bface->oZi;
      tBox *obox = NULL;
      tBface *obface = NULL;
      int touch, sameX, sameY, sameZ;
      int pi, ind;
      double *oX;
      double *oY;
      double *oZ;
      int of, od;

      if(ob<0)  continue; /* do nothing if there is no other box at this bface */
      if(ofi<0) continue; /* do nothing if there not one face index in other box */

      /* other box and corresponding bface */
      obox = grid->box[ob];
      obface = obox->bface[ofi];

      if(obface->ofi < 0) continue; /* nothing if not one face index */

      /* if we get here there are two bfaces that may touch */

      /* mark as touching if bbox values agree */
      if( dequal(box->bbox[bface->f], obox->bbox[obface->f]) )
      {
        double bX,bY,bZ;
        int i,j,k, off_face;
        int f = bface->f;
        int dir = 1 + f/2;

        bface->touch = 1; /* we are touching */

        /* is point list in fpts on face? */
        i=j=k = off_face = 0;
        if(dir==1)
        {
          i = ( (n1-1) )*(f%2);
          bX = X[Index(i,j,k)];
          if( !dequal(bX, box->bbox[bface->f]) ) off_face=1;
        }
        else if(dir==2)
        {
          j = ( (n2-1) )*(f%2);
          bY = Y[Index(i,j,k)];
          if( !dequal(bY, box->bbox[bface->f]) ) off_face=1;
        }
        else /* dir=3 */
        {
          k = ( (n3-1) )*(f%2);
          bZ = Z[Index(i,j,k)];
          if( !dequal(bZ, box->bbox[bface->f]) ) off_face=1;
        }
        bface->fpts_off_face = off_face;

        /* should we set fields or their nomral derivs? */
        if(bface->setnormalderiv == 0)  obface->setnormalderiv = 1;
        else                            obface->setnormalderiv = 0;
      }

      /* do nothing if bface->fpts is NULL */
      if(bface->fpts == NULL) continue;

      /* do nothing if boxes have a different number of grid points */
      if(box->nnodes != obox->nnodes) continue;

      /* loop over bface, and check which coords agree */
      oX = obox->v[iX];
      oY = obox->v[iY];
      oZ = obox->v[iZ];
      of = obface->f;
      od = 1+of/2;
      sameX = sameY = sameZ = 1;
      forPointList_inbox(bface->fpts, box, pi, ind)
      {
        int ok = kOfInd_n1n2(ind, n1,n2);
        int oj = jOfInd_n1n2_k(ind, n1,n2, ok);
        int oi = iOfInd_n1n2_jk(ind, n1,n2, oj,ok);
        int oind;

        if(od==1)      oi = ( (obox->n1-1) )*(of%2);
        else if(od==2) oj = ( (obox->n2-1) )*(of%2);
        else if(od==3) ok = ( (obox->n3-1) )*(of%2);
        oind = Ind_n1n2(oi,oj,ok, obox->n1,obox->n2);

        if( !dequal(oX[oind], X[ind]) ) sameX = 0;
        if( !dequal(oY[oind], Y[ind]) ) sameY = 0;
        if( !dequal(oZ[oind], Z[ind]) ) sameZ = 0;
        if( !(sameX || sameY || sameZ) ) break;
      }
      bface->sameX = sameX;
      bface->sameY = sameY;
      bface->sameZ = sameZ;
    }
  }

  return 0;
}
