#include "mex.h"
#include <math.h>
#include <time.h>

#define DEBUG 0

/*
 * Efficient computation of a similarity measure between the rows of two matrices
 
 in:
 k x m  - matrix A
 l x m  - matrix B
 string - measure
 
 out:
 k x l - similarity

 * algorithm:
 * loop over rows of A
 * - for each row of A, loop over rows of B
 *   - compute similarity between row of A and row of B
 *   - similarity:
 *     - correlation: (preprocess data to remove mean and standard deviation),
 *                     and this becomes dot product with scaling
 *     - cosine: (preprocess to make norm 1), becomes dot product with scaling
 *     - euclidean: add up squared distances
 *
 * notes:
 * - euclidean/cosine don't work yet
 * - no error handling - this is meant to be called O(10K) times, all handling is done in the caller
 * history:
 * 2009 Nov 23 - created - fpereira@princeton.edu
 *
*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /* input */
  double *A; double *Ap; mxArray *Ap_mx;
  double *B; double *Bp; mxArray *Bp_mx;
  int *mptr;

  /* output */
  double *S; /* similarity matrix */
  double *indexMostSimilar; /* index (starts at 1, MATLAB convention, of closest example)*/

  /* rest of variables pertaining to input/output contents */
  int nrowsA, ncolsA, nrowsB, ncolsB;
  int measurenum; /* numeric version of the measure */
  double rowmean, rowstdv;
  double *rowptr; double *rowptr2; double *rowptrA; double *rowptrB; double *rowptrS;
  int ia, ib, is;
  int ja, jb, js;
  double *valueMostSimilar;
  int noneighbourhack; /* use if voxel has no neighbors, so input matrices have 1 column */

  /* everything else */
  int e,i,ig,j,k,v;

  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=3) 
    mexErrMsgTxt("three inputs required");
  if(nlhs!=3) 
    mexErrMsgTxt("three outputs required");
  
  /* check to make sure the input arguments are double matrices */
  if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) ) {
    mexErrMsgTxt("first two inputs must be matrices");
  }
  
  /*  create pointers to the inputs */
  A = mxGetPr(prhs[0]);
  B = mxGetPr(prhs[1]);
  mptr = (int *) mxGetPr(prhs[2]); 
  measurenum = *mptr;

  nrowsA = mxGetN(prhs[0]);  ncolsA = mxGetM(prhs[0]); /* dimensions A */
  nrowsB = mxGetN(prhs[1]);  ncolsB = mxGetM(prhs[1]); /* dimensions B */

  if ((ncolsA==1) && (measurenum==0) ) {
    /* if a voxel has no neighbours, matrix A will have one column,
       compute euclidean distance instead and then turn to similarity */
    noneighbourhack = 1;
    measurenum = 1;
  } else {
    noneighbourhack = 0;
  }

  /* create output matrices */
  plhs[0] = mxCreateDoubleMatrix(nrowsB,nrowsA,mxREAL);
  S = mxGetPr(plhs[0]);

  plhs[1] = mxCreateDoubleMatrix(nrowsA,1,mxREAL);
  indexMostSimilar = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(nrowsA,1,mxREAL);
  valueMostSimilar = mxGetPr(plhs[2]);

  /** preprocess for each measure **/
  
  switch (measurenum) {
  case 0: /* correlation */
    /* temporary space for input matrices after preprocessing */
    Ap = (double *) mxMalloc( nrowsA*ncolsA*sizeof(double) );
    Bp = (double *) mxMalloc( nrowsB*ncolsB*sizeof(double) );

    for (ia=0; ia<nrowsA; ia++) { valueMostSimilar[ia] = -1000; };

    /* subtract the mean of each row and divide by its standard deviation */
    for (i=0; i<nrowsA; i++) {
      rowptr  = A  + i*ncolsA;
      rowptr2 = Ap + i*ncolsA;
      rowmean = 0; rowstdv = 0;
      
      /* compute row mean */
      for (j=0; j<ncolsA; j++) { rowmean = rowmean + rowptr[j]; }
      rowmean = rowmean / ncolsA;
      
      /* compute row stdv */
      for (j=0; j<ncolsA; j++) { rowstdv = rowstdv + pow(rowptr[j]-rowmean,2); }
      rowstdv = sqrt(rowstdv / (ncolsA-1));
      
      /* subtract mean and divide by standard deviation while copying*/
      for (j=0; j<ncolsA; j++) {
	  rowptr2[j] = (rowptr[j] - rowmean) / rowstdv;
      }
    }
    
    /* subtract the mean of each row and divide by its standard deviation */
    for (i=0; i<nrowsB; i++) {
      rowptr  = B  + i*ncolsB;
      rowptr2 = Bp + i*ncolsB;
      rowmean = 0; rowstdv = 0;
      
      /* compute row mean */
      for (j=0; j<ncolsB; j++) { rowmean = rowmean + rowptr[j]; }
      rowmean = rowmean / ncolsB;
      
      /* compute row stdv */
      for (j=0; j<ncolsB; j++) { rowstdv = rowstdv + pow(rowptr[j]-rowmean,2); }
      rowstdv = sqrt(rowstdv / (ncolsB-1));
      
      /* subtract mean and divide by standard deviation while copying*/
      for (j=0; j<ncolsB; j++) {
	rowptr2[j] = (rowptr[j] - rowmean) / rowstdv;
      }
    }
    break;
  
  case 1: /* euclidean */
    /* no need for temporary space, just set the pointer to the originals */
    Ap = A;
    Bp = B;

    for (ia=0; ia<nrowsA; ia++) { valueMostSimilar[ia] = DBL_MAX; };
    break;
  }

  /* debug, make sure the pointers are different
  printf("Ap=%p A=%p\n",Ap,A);
  printf("Bp=%p B=%p\n",Bp,B);
  fflush(stdout); return;
  */

  /* debug, print matrices */
  if (0) {
    printf("measure %d\n",measurenum);

    printf("matrix A:\n");
    for (ia=0; ia<nrowsA; ia++) {
      rowptrA = Ap + ia*ncolsA;
      for (ja=0; ja<ncolsA; ja++) {
	printf("%1.4f ",rowptrA[ja]);
      }
      printf("\n");
    }

    printf("matrix B:\n");
    for (ib=0; ib<nrowsB; ib++) {
      rowptrB = Bp + ib*ncolsB;
      for (jb=0; jb<ncolsB; jb++) {
	printf("%1.4f ",rowptrB[jb]);
      }
      printf("\n");
    }
  }


  /** compute measure between each example of A and all examples of B */

  for (ia=0; ia<nrowsA; ia++) {
    rowptrA = Ap + ia*ncolsA;

    rowptrS = S + ia*nrowsB;

    for (js=0; js<nrowsB; js++) { rowptrS[js] = 0; }

    for (ib=0; ib<nrowsB; ib++) {
      rowptrB = Bp + ib*ncolsB;

      rowptrS[ib] = 0;

      switch (measurenum) {
      case 0:
	/* correlation */
	for (js=0; js<ncolsB; js++) {
	  rowptrS[ib] = rowptrS[ib] + rowptrA[js]*rowptrB[js];
	}
	rowptrS[ib] = rowptrS[ib] / (ncolsB-1);

	if (rowptrS[ib] > valueMostSimilar[ia]) {
	  valueMostSimilar[ia] = rowptrS[ib];
	  indexMostSimilar[ia] = (double) ib+1;
	}
	
	break;
      case 1:
	/* euclidean */
	for (js=0; js<ncolsB; js++) {
	  rowptrS[ib] = rowptrS[ib] + pow(rowptrA[js]-rowptrB[js],2);
	}
	rowptrS[ib] = sqrt(rowptrS[ib]);

	if (rowptrS[ib] < valueMostSimilar[ia]) {
	  valueMostSimilar[ia] = rowptrS[ib];
	  indexMostSimilar[ia] = (double) ib+1;
	}

	break;
      }
    } /* for over B rows */
  } /* for over A rows */

  if (noneighbourhack == 1) {
    /* no need to fix anything here, as the only output used is indexMostSimilar[ia] */
    for (ia=0; ia<nrowsA; ia++) {
      rowptrA = Ap + ia*ncolsA;
    }
  }

  /* printf("aqui\n"); fflush(stdout); return; */
  
  /* get rid of all the matrices obtained with mxMalloc (in the order the mxMallocs appear in the code) */

  switch (measurenum) {
  case 0:
    mxFree(Ap); mxFree(Bp);
    break;
  case 1:
    /* nothing to free */
    break;
  }

}
