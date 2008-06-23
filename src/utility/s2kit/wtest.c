/***************************************************************************
  **************************************************************************
  
                           S2kit 1.0

          A lite version of Spherical Harmonic Transform Kit

   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
   Copyright 2004 Peter Kostelec, Dan Rockmore

   This file is part of S2kit.

   S2kit is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   S2kit is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with S2kit; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


/*************************************************************************/

/* test_naive.c - top-level code for computing forward and inversee
   Legendre transforms using the naive algorithm; will return timing
   and error information

   m    - order of the problem
   bw   - bandwidth
   loops - number of loops thru timed portion of code.  Intended
           to reduce noise due to multiprocessing and 
	   discretization errors
   
   Sample calls:

   test_naive m bw loops

   test_naive 0 32 10

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "s2kit_pmls.h"
#include "s2kit_makeweights.h"
#include "s2kit_naive_synthesis.h"

#define mymax(a, b) ((a) > (b) ? (a) : (b))


int main(argc, argv)
     int argc;
     char **argv;
{
  int i, j, k;
  int bw, m;
  double *samples, *coeffs, *newcoeffs ;
  double *plm, *weights, *workspace ;

  if (argc < 3)
    {
      fprintf(stdout,"Usage: test_naive m bw\n");
      return(0);
    }

  m = atoi(argv[1]);
  bw = atoi(argv[2]);

  /* space for samples and coefficients */
  samples = (double *) malloc(sizeof(double) * 2 * bw);
  coeffs = (double *) malloc(sizeof(double) * (bw - m) );
  newcoeffs = (double *) malloc(sizeof(double) * (bw - m) );

  /* space for precomputed Plms */
  plm = (double *) malloc(sizeof(double) * 2 * bw * (bw - m) );

  /* for weights */
  weights = (double *) malloc(sizeof(double) * 4 * bw);

  /* workspace space */
  workspace = (double *) malloc(sizeof(double) * 18 * bw);

  /* precompute the Plms */
  PmlTableGen( bw, m, plm, workspace ) ;

  /* make the weights */
  makeweights( bw, weights );

  /* generate random coefficients */
  for( i = 0 ; i < (bw - m) ; i ++ )
  {
	/* coeffs[ i ] = 2.0 * ( drand48() - 0.5 ) ; */
	coeffs[ i ] = i/100.0;
	printf("coeffs[ i ] = %f\n", coeffs[ i ]);
  }
  printf("\n");

  /* do inverse naive transform */
  Naive_SynthesizeX(coeffs,
			bw,
			m,
			samples,
			plm ) ;

  /* now do forward naive transform */
  Naive_AnalysisX( samples,
		       bw,
		       m,
		       weights,
		       newcoeffs,
		       plm,
		       workspace );
  /*
	now tally up the error between the original
	coefficients and the new ones
  */

  /* generate random sample */
  for( i = 0 ; i < 2*bw; i ++ )
  {
    double th = 3.141592653589793*(2*i+1)/(4*bw);
    samples[i] = sin(th);
    printf("samples[%d] = %f\n", i, samples[ i ]);
  }
  printf("\n");

  /* now do forward naive transform */
  Naive_AnalysisX( samples,
		       bw,
		       m,
		       weights,
		       coeffs,
		       plm,
		       workspace );

  /* do inverse naive transform */
  Naive_SynthesizeX(coeffs,
			bw,
			m,
			samples,
			plm ) ;
  /* print new */
  for( i = 0 ; i < 2*bw; i ++ )
  {
	printf("new samples[%d] = %f\n", i, samples[ i ]);
  }
  printf("\n");


  fprintf(stdout,"bw = %d   m = %d\n",bw, m);

  /* now free up memory */
  free( workspace );
  free( weights );
  free( newcoeffs );
  free( coeffs );
  free( samples );

  return 0 ;
}










