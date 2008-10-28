/* OpenMPsave_loops.h */
/* Wolfgang Tichy 10/2008 */

/* some loops needed in derivs.c and other places */
/* to use with OpenMP do e.g. something like this
   #pragma omp parallel for
   forLines_alloc2Lines(j,k, n2,n3, uline,duline,n1)
   {
     get_memline(u, uline, 1, j, k, n1, n2, n3);
     matrix_times_vector(box->DD1, uline, duline, n1);
     put_memline(du, duline, 1, j, k, n1, n2, n3);
   } end_forLines_free2Lines(uline,duline)
*/

/* loop over two directions */
#define forLines(j,k, n2,n3)\
  for((k)=0; (k)<(n3); (k)++ )\
  {\
    int (j);\
    for((j)=0; (j)<(n2); (j)++)

/* end the loop forLines */
#define end_forLines }

/* loop over two directions and allocate 2 lines in the third direction */
#define forLines_alloc2Lines(j,k, n2,n3, uline,duline,n1)\
  for((k)=0; (k)<(n3); (k)++ )\
  {\
    int (j);\
    double *(uline)  = (double*)  calloc((n1), sizeof(double));\
    double *(duline) = (double*)  calloc((n1), sizeof(double));\
    for((j)=0; (j)<(n2); (j)++)

/* end the loop forLines_alloc2Lines */
#define end_forLines_free2Lines(uline,duline) free((uline)); free((duline));}

/* loop over two directions and allocate 3 lines in the third direction */
#define forLines_alloc3Lines(j,k, n2,n3, uline,duline,dduline,n1)\
  for((k)=0; (k)<(n3); (k)++ )\
  {\
    int (j);\
    double *(uline)  = (double*)  calloc((n1), sizeof(double));\
    double *(duline) = (double*)  calloc((n1), sizeof(double));\
    double *(dduline)= (double*)  calloc((n1), sizeof(double));\
    for((j)=0; (j)<(n2); (j)++)

/* end the loop forLines_alloc3Lines */
#define end_forLines_free3Lines(uline,duline,dduline)\
  free((uline));  free((duline));  free((dduline));}
