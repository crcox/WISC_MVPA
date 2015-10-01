#include "mex.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define DEBUG 0

/*
 * Efficient crossvalidation and permutation test searchlight GNB
 * (fpereira@princeton.edu)
 
 in:
 n x m  - examples;
 n x 1  - labels;
 n x 1  - labelsGroup;
 m x <radius^3 -1> - voxelsToNeighbours;
 m x 1  - numberOfNeighbours;
 
 out:
 1 x m  - error
 1 x m  - fraction
 
 assumes:
 - examples are ordered by group
 - group numbers are 1:#groups, all of them occurring
 - class labels are 1:#classes, all of them occurring
 
*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /* input */
  double *examples;
  int *labels;             double *plabels; /* original labels */
  int *labelsGroup;        double *plabelsGroup;
  int *voxelsToNeighbours; double *pvoxelsToNeighbours;
  int *numberOfNeighbours; double *pnumberOfNeighbours;
  int nPermutations;       double *pnPermutations;

  /* output */
  double *accuracy;
  double *accuracyCount;
  double *fraction;

  /* rest of variables pertaining to input/output contents */
  int nGroups;
  int itmp;
  int iprev,iseen;
  int group;
  int *groupStarts; double *pgroupStarts; mxArray *groupStarts_mx;
  int *groupEnds;   double *pgroupEnds;   mxArray *groupEnds_mx;
  int *groupSize;   double *pgroupSize;   mxArray *groupSize_mx;
  int *groupNclass; double *pgroupNclass; mxArray *groupNclass_mx;
  int **groupPermutations;
  int p;
  int nClasses;
  int n,m;
  int g,c,gidx;
  int *gptr;    /* group pointer */
  double *eptr, *e2ptr; /* example pointer */
  int ngc;
  double *sumX1; mxArray *sumX1_mx; double *x1ptr;
  double *sumX2; mxArray *sumX2_mx; double *x2ptr;
  double *mptr,*m2ptr;
  double *sptr;
  double *exampleLargestValue;
  double *examplesCorrect; mxArray *examplesCorrect_mx;
  double *examplesSquared;
  double *means;
  double *meansSquared;
  double *vars;
  double *sumX1train;
  double *sumX2train;
  int gtest, gtrain;
  int *nc;
  int ntest, ntrain, ntrainmo;
  int nn;
  int *nptr, *iptr;
  int neighbour;
  double voxelValue;
  double *cptr;
  int neighbourRadius; double *pneighbourRadius; int neighbourMax;
  
  /* everything else */
  int     e,i,ig,j,k,v;

  double *X,*M,*Ez,*Vz,*sigma;
  double *pEz,*pVz,*pX,*pM,*pL;
  double *pEzM, *pEzM2;
  double **mEz,**mVz,**mX,**mM,**mL;
  double fixed;
  double *log_px;
  double *log_px_pd;
  mxArray *tmp1_mx,*tmp2_mx,*sigma2_mx,*logsigma_mx,*pEzM_mx,*pEzM2_mx;
  double  *tmp1,*tmp2,*sigma2,*logsigma;
  int     status,mrows,ncols;
  int     cX,rX,cM,rM,cEz,rEz,cVz,rVz,cS,rS;
  int     idx;
  char    buffer[10];

  struct timeval tvStart; struct timeval tvEnd; struct timeval tvDuration;
  double  mduration;


  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=7) 
    mexErrMsgTxt("seven inputs required");
  if(nlhs!=3) 
    mexErrMsgTxt("three outputs required");
  
  /* check to make sure the input arguments are double matrices */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("inputs must be matrices");
  }
  
  /*  create a pointer to the input matrices */
  examples            = mxGetPr(prhs[0]);
  plabels             = mxGetPr(prhs[1]);
  plabelsGroup        = mxGetPr(prhs[2]);
  pneighbourRadius    = mxGetPr(prhs[3]);
  pvoxelsToNeighbours = mxGetPr(prhs[4]);
  pnumberOfNeighbours = mxGetPr(prhs[5]);
  pnPermutations      = mxGetPr(prhs[6]); nPermutations = round(*pnPermutations);
  n  = mxGetN(prhs[0]);  m = mxGetM(prhs[0]); /* dimensions */
  

  /* mexPrintf("#examples=%d m=%d\n",n,m);fflush(stdout); */
  /* mexPrintf("pvn = %d,%d pnn = %d,%d\n",mxGetN(prhs[3]),mxGetM(prhs[3]),mxGetN(prhs[4]),mxGetM(prhs[4]));fflush(stdout);*/

  /*  create output matrices */
  plhs[0] = mxCreateDoubleMatrix(1,m,mxREAL);
  accuracy = mxGetPr(plhs[0]);
  accuracyCount = (double *) mxMalloc( m*sizeof(double) );
  for (j=0; j<m; j++) { accuracyCount[j] = 0; /* no need to initialize accuracy */ };

  /* used both to return the binary correct/incorrect matrix but also to store intermediate classification results */
  plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
  examplesCorrect = mxGetPr(plhs[1]);
  for (i=0; i<n; i++) { cptr = examplesCorrect+i*m; for (j=0; j<m; j++) { cptr[j] = 0; }};

  /* used to return the fraction if a permutation test was asked for
     keeps counts during processing and starts at 1, since the correct labelling must be counted */
  plhs[2] = mxCreateDoubleMatrix(1,m,mxREAL);
  fraction = mxGetPr(plhs[2]);
  for (j=0; j<m; j++) { fraction[j] = 1; }


  /** figure out things and create temporary matrices */

  gettimeofday(&tvStart,NULL);

  /* some of the inputs are really matrices of integers, convert from double to a local integer copy */

  labels     = (int *) mxMalloc( n*sizeof(int) );
  for (i=0; i<n; i++) { labels[i] = (int) plabels[i]; }

  /* how many classes do we have (assume they are already 1:C), convert to 0:(C-1) */
  nClasses = 0;
  for (i=0; i<n; i++)
    {
      if (labels[i] > nClasses) { nClasses = labels[i]; }
      labels[i] = labels[i] - 1;
    }

  labelsGroup = (int *) mxMalloc( n*sizeof(int) );
  for (i=0; i<n; i++) { labelsGroup[i] = (int) plabelsGroup[i]; }
  
  /* how many groups do we have and where do they start/end */

  nGroups = 0; iprev = -1;
  for (i=0; i<n; i++){ if (labelsGroup[i] != iprev) { nGroups++; iprev = labelsGroup[i]; } }

  groupStarts = (int *) mxMalloc( nGroups*sizeof(int) );
  groupEnds   = (int *) mxMalloc( nGroups*sizeof(int) );
  groupSize   = (int *) mxMalloc( nGroups*sizeof(int) );

  /* for (i=0; i<n; i++) { mexPrintf("%d\t%d\t%d\n",i,labels[i], labelsGroup[i]);fflush(stdout); } */

  iprev = 0; group = 0;
  for (i=0; i<n; i++)
    {
      if (labelsGroup[i] != iprev)
	{ /* new group starts here */
	  /* mexPrintf("new group at %d\n",i);fflush(stdout);*/
	  groupStarts[group] = i;
	  if (group > 0) { groupEnds[group-1] = i-1; }
	  group++;
	  iprev = labelsGroup[i];
	}
    }
  groupEnds[group-1] = n-1;
  for (g=0; g<nGroups; g++) { groupSize[g] = groupEnds[g]-groupStarts[g]+1; }

  /* for (g=0; g<nGroups; g++) { mexPrintf("%d\t%d\t%d\n",groupStarts[g],groupEnds[g],groupSize[g]);}; fflush(stdout);*/

  /* count how many examples of each class in each group  */
  groupNclass = (int *) mxMalloc( nGroups*nClasses*sizeof(int) );
  
  for (g=0; g<nGroups; g++) {
    gptr = groupNclass + g*nClasses;
    for (c=0; c<nClasses; c++) { gptr[c] = 0; }
  
    /* scroll over examples in this group and tally the different classes */
    for (i=groupStarts[g]; i<=groupEnds[g]; i++) {
      gptr[labels[i]] = gptr[labels[i]] + 1;
    }
  }

  /* turns neighbour information arguments into arrays of ints, also switch to 0-based indexing */
  neighbourRadius = (int) *pneighbourRadius;
  neighbourMax    = pow( 2*neighbourRadius+1, 3)-1;

  numberOfNeighbours = (int *) mxMalloc( m*sizeof(int) );
  for (j=0; j<m; j++) { numberOfNeighbours[j] = (int) floor(pnumberOfNeighbours[j]); }

  voxelsToNeighbours = (int *) mxMalloc( m*neighbourMax*sizeof(int) );
  for (j=0; j<m; j++) {
    nptr = voxelsToNeighbours  + j*neighbourMax;
    mptr = pvoxelsToNeighbours + j*neighbourMax; 
    for (v=0; v<neighbourMax; v++) { nptr[v] = (int) floor(mptr[v]-1); }
  }
  

  /* debug examples
     for (i=0; i<10; i++) {
       eptr = examples + i*m;
       for (j=0; j<8; j++) { mexPrintf("%1.2f ",eptr[j]); }; mexPrintf("\n"); fflush(stdout);
     }
  */

  /* debug information about #groups/#classes
  for (g=0; g<nGroups; g++) {
    gptr = groupNclass + g*nClasses;
    for (c=0; c<nClasses; c++) { 
      mexPrintf("%d ",gptr[c]);
    }
    mexPrintf("\n"); fflush(stdout);
  }
  */
  mexPrintf("#examples=%d #voxels=%d #groups=%d #classes=%d\n",n,m,nGroups,nClasses); fflush(stdout);


  /*
  for (j=0; j<m; j++)
    {
      nn   = numberOfNeighbours[j];
      nptr = voxelsToNeighbours + j*neighbourMax;

      mexPrintf("neighbours %d:\t",j);fflush(stdout);
      for (v=0; v<nn; v++)
	{
	  neighbour = nptr[v];
	  mexPrintf("%d ",neighbour);fflush(stdout);
	}
      mexPrintf("\n");fflush(stdout);
    }
  */

  /*
    for (j=0; j<m; j++)
      {
        nptr = voxelsToNeighbours + j*neighbourMax;
        mexPrintf("neighbours %d:\t",j);fflush(stdout);
        for (v=0; v<neighbourMax; v++) { mexPrintf("%d ",nptr[v]);fflush(stdout); }; mexPrintf("\n");fflush(stdout);
      }
  */

  /* debug neighbour info
  for (j=0; j<m; j++)
    {
      nn   = numberOfNeighbours[j];
      nptr = voxelsToNeighbours + j*neighbourMax;

      mexPrintf("neighbours %d:\t",j);fflush(stdout);
      for (v=0; v<nn; v++)
	{
	  neighbour = nptr[v];
	  mexPrintf("%d ",neighbour);fflush(stdout);
	}
      mexPrintf("\n");fflush(stdout);
    }
  */


  /* get all the memory allocation out of the way here (only needs to be done once) */

  ngc = nGroups * nClasses;

  sumX1 = mxMalloc( m*ngc*sizeof(double));
  sumX2 = mxMalloc( m*ngc*sizeof(double));
  examplesSquared = mxMalloc( n*m*sizeof(double));

  nc = mxMalloc( nClasses*sizeof(double) );

  means = mxMalloc( nClasses*m*sizeof(double) );
  vars  = mxMalloc( 1*m*sizeof(double) );
  meansSquared = mxMalloc( nClasses*m*sizeof(double) );

  sumX1train = mxMalloc( nClasses*m*sizeof(double) );
  sumX2train = mxMalloc( nClasses*m*sizeof(double) );

  exampleLargestValue = mxMalloc( 1*m*sizeof(double) );

  gettimeofday(&tvEnd,NULL);

  timersub(&tvEnd,&tvStart,&tvDuration);
  mduration = (double) tvDuration.tv_usec / 1000000;
  mexPrintf("finished setting up in %lf second(s)",mduration);
  mexPrintf("\n");fflush(stdout);


  /**
   ** compute the sum of examples and examples^2 for each group/class
   ** (store in array which has #groups x #classes rows, all classes for each group)
   **/

  gettimeofday(&tvStart,NULL);

  for (g=0; g<nGroups; g++)
    {
      #if DEBUG 
      mexPrintf("processing group %d\n",g);fflush(stdout);
      #endif

      /* base of the matrix for this group (#classes x #voxels) */
      x1ptr = sumX1 + g*(nClasses*m);
      x2ptr = sumX2 + g*(nClasses*m);

      /* initialize the respective matrix portion to 0 */

      for (c=0; c<nClasses; c++) {
	mptr = x1ptr + c*m;
	for (j=0; j<m; j++) { mptr[j] = 0; }
      }

      for (c=0; c<nClasses; c++) {
	mptr = x2ptr + c*m;
	for (j=0; j<m; j++) { mptr[j] = 0; }
      }
      /* compute the sums of examples, squares of examples and sums of squares of examples */

      for (i=groupStarts[g]; i<=groupEnds[g]; i++) {
	eptr  = examples + i*m;
	e2ptr = examplesSquared + i*m;
	
	#if DEBUG
	mexPrintf("\texample %d\tclass %d\n",i,labels[i]);fflush(stdout);
	#endif

	/* add this example to the appropriate class tally */
	mptr = x1ptr + labels[i]*m;
	for (j=0; j<m; j++) { mptr[j] = mptr[j] + eptr[j]; }

	/* square it */
	for (j=0; j<m; j++) { e2ptr[j] = eptr[j] * eptr[j]; }

	/* add the squared of the example to the appropriate class tally */
	mptr = x2ptr + labels[i]*m;
	for (j=0; j<m; j++) { mptr[j] = mptr[j] + e2ptr[j]; }
      } 
    }

  /* debug tallies (X1 and X2, comment either x1ptr or x2ptr out)
    for (g=0; g<nGroups; g++) {
      x1ptr = sumX1 + g*(nClasses*m);
      x2ptr = sumX2 + g*(nClasses*m);
      mexPrintf("group %d\n",g);
      for (c=0; c<nClasses; c++) { 
        mptr = x1ptr + c*m;
        mptr = x2ptr + c*m;
        mexPrintf("\t"); for (j=0; j<8; j++) { mexPrintf("%1.2f ",mptr[j]); }; mexPrintf("\n");
      };fflush(stdout);
    }
  */

  /**
   ** Main cross-validation loop, in every loop
   ** - add up the x1 and x2 matrices for groups in the training set
   ** - derive class mean and standard deviations from that
   **/
  
  
  for (gtest=0; gtest<nGroups; gtest++)
    {
      #if DEBUG
      mexPrintf("testing group %d\n",gtest); fflush(stdout);
      #endif

      /* rezero the space for class means, stdv, sum X1, sum X2 and nc */      
      for (c=0; c<nClasses; c++) { mptr = means + c*m; for (j=0; j<m; j++) { mptr[j] = 0; }; }
      for (c=0; c<nClasses; c++) { mptr = meansSquared + c*m; for (j=0; j<m; j++) { mptr[j] = 0; }; }
      mptr = vars;  for (j=0; j<m; j++) { mptr[j] = 0; }
      for (c=0; c<nClasses; c++) { mptr = sumX1train + c*m; for (j=0; j<m; j++) { mptr[j] = 0; } }
      for (c=0; c<nClasses; c++) { mptr = sumX2train + c*m; for (j=0; j<m; j++) { mptr[j] = 0; } }
      for (c=0; c<nClasses; c++) { nc[c] = 0; }

      gidx = 0; /* counter of # groups being added */
      
      /* sum X1 and X2 over the other groups, producing a #classes x #voxels array for each */
      for (gtrain=0; gtrain<nGroups; gtrain++)
	{
	  if (gtrain != gtest)
	    {
	      /* add counts of each class in this group */
	      gptr = groupNclass + gtrain*nClasses;
	      for (c=0; c<nClasses; c++) { nc[c] = nc[c] + gptr[c]; }
	      
	      /* start of the X1 and X2 arrays for this group */
	      x1ptr = sumX1 + gtrain*(nClasses*m);
	      x2ptr = sumX2 + gtrain*(nClasses*m);
	      
	      /* add up X1 */
	      for (c=0; c<nClasses; c++)
		{
		  mptr = sumX1train + c*m;
		  for (j=0; j<m; j++) { mptr[j] = mptr[j] + x1ptr[j]; }
		  x1ptr = x1ptr + m;
		}
	      
	      /* add up X2 */
	      for (c=0; c<nClasses; c++)
		{
		  mptr = sumX2train + c*m;
		  for (j=0; j<m; j++) { mptr[j] = mptr[j] + x2ptr[j]; }
		  x2ptr = x2ptr + m;
		}
	
	      gidx = gidx + 1;
	    }
	} /* end of for tallying over training groups */
      ntrain = 0; for (c=0; c<nClasses; c++) { ntrain = ntrain + nc[c]; }; ntrainmo = ntrain - 1;

      /* compute means and meansSquared (one vector per class) */
      /* compute variances (one vector for all classes) */

      for (c=0; c<nClasses; c++)
	{
	  x1ptr = sumX1train + c*m;
	  x2ptr = sumX2train + c*m;
	  mptr  = means + c*m;
	  sptr  = meansSquared + c*m;

	  for (j=0; j<m; j++)
	    {
	      mptr[j] = x1ptr[j] / nc[c];
	      sptr[j] = mptr[j] * mptr[j]; /* compute the means squared to use below */
	      vars[j] = vars[j] + x2ptr[j] - 2*mptr[j]*x1ptr[j] + sptr[j] * nc[c];
	    }
	}

      for (j=0; j<m; j++) { vars[j] = vars[j] / ntrainmo; }

      /* add to the tally */

      /** 
       ** apply to test set (scroll over examples in the test group)
       ** 
       **
       **/
      
      /* create a buffer for each example where we store the current largest value,
	 class assignment is put in examplesCorrect
	 at the end check the label in examplesCorrect and set to 1 or 0
      */
      
      for (i=groupStarts[gtest]; i<=groupEnds[gtest]; i++)
	{
          #if DEBUG
	  mexPrintf("\ttesting example %d\n",i); fflush(stdout);
	  #endif

	  eptr  = examples        + i*m;
	  e2ptr = examplesSquared + i*m;	  
	  cptr  = examplesCorrect + i*m; /* while going over classes stores prediction, then binary 1/0 for right/wrong */
	  mptr  = exampleLargestValue; for (j=0; j<m; j++) { mptr[j] = -1000000000; }

	  /* main loop over voxels (within it, for every voxel, we loop over neighbours) */
	  for (j=0; j<m; j++)
	    {
	      nn   = numberOfNeighbours[j];
	      nptr = voxelsToNeighbours + j*neighbourMax;
      
	      for (c=0; c<nClasses; c++)
		{
		  mptr  = means + c*m;
		  m2ptr = meansSquared + c*m;

		  #if DEBUG
		  mexPrintf("class %d\n",c); fflush(stdout);
		  #endif

		  /* compute p(x,x_neighbours|c), minus parts that are same for all classes */
		  voxelValue = (2*eptr[j]*mptr[j] - m2ptr[j] - e2ptr[j])/vars[j];

		  #if DEBUG
		  mexPrintf("\tvoxel %d\n\t\tneighbours: ",j); fflush(stdout);
		  #endif
		  
		  for (v=0; v<nn; v++)
		    {
		      neighbour = nptr[v];
                      #if DEBUG
		      mexPrintf("%d ",neighbour);fflush(stdout);
		      #endif
		      voxelValue += (2*eptr[neighbour]*mptr[neighbour] - m2ptr[neighbour] - e2ptr[neighbour])/vars[neighbour];
		    }

                  #if DEBUG
		  mexPrintf("\n\tneighbours done\n"); fflush(stdout);
		  #endif

		  /* mexPrintf("%1.4f\t%1.4f\t",voxelValue,exampleLargestValue[j]); fflush(stdout); */

		  /* if p(x,x_neighbours|c) is larger than what we've seen so far, update prediction */
		  if (voxelValue > exampleLargestValue[j]) {
		    exampleLargestValue[j] = voxelValue;
		    cptr[j] = (double) c;
		    /* mexPrintf("better!\t");fflush(stdout); */
		  } else {
		    /* mexPrintf("worse!\t");fflush(stdout); */
		  }
		  /* mexPrintf("%1.2f\n",cptr[j]); */
		  
                  #if DEBUG
		  mexPrintf("\tclassified\n"); fflush(stdout);
                  #endif
		  
	        } /* end of loop over classes */
      

              /* transform to correct/incorrect */
              if ( cptr[j] != ((double) labels[i]) ) { cptr[j] = 0; } else { cptr[j] = 1; }
	      accuracyCount[j] = accuracyCount[j] + cptr[j];
	      
          } /* end of loop over voxels */

	  #if DEBUG
          for (j=0; j<20; j++) { mexPrintf("%d ",(int) cptr[j]);}; mexPrintf("\n");fflush(stdout);
          #endif

      } /* end of loop over test examples */

      /* debug: show the means
        for (c=0; c<nClasses; c++) {
	  mptr = means + c*m;
	  mexPrintf("%d:\t",c); for (j=0; j<8; j++) { mexPrintf("%1.2f ",mptr[j]); }; mexPrintf("\n");
        }
      */

      /* debug: show the variances
        mptr = vars;
        mexPrintf("var:\t"); for (j=0; j<8; j++) { mexPrintf("%1.2f ",mptr[j]); }; mexPrintf("\n");
      */

      /* debug examples squared 
        for (i=0; i<10; i++) {
	  sptr = examplesSquared + i*m;
	  for (j=0; j<8; j++) { mexPrintf("%1.2f ",sptr[j]); }; mexPrintf("\n"); fflush(stdout);
        }
      */
      
      /*mexPrintf("test %d:\t",gtest); for (c=0; c<nClasses; c++) { mexPrintf("%d ",nc[c]); };mexPrintf("\n");fflush(stdout); */

    }  /* end of for over test groups */
  

  /** finish tallying **/
  for (j=0; j<m; j++) { accuracy[j] = accuracyCount[j]; } /* divide later, keep counts for now */

  gettimeofday(&tvEnd,NULL);

  /* give the user an estimate of running time, so that they can abort if necessary */
  timersub(&tvEnd,&tvStart,&tvDuration);
  mduration = (double) tvDuration.tv_usec / 1000000;

  mexPrintf("finished running one iteration in %lf second(s)",mduration);
  if (nPermutations) { mexPrintf(", estimated running time for the permutation test (%d permutations) is %d seconds.",nPermutations,(int) round(mduration*nPermutations));}
  else { mexPrintf("\tall done!"); }
  mexPrintf("\n");fflush(stdout);

  if (nPermutations) {
    
    /**
     ** Permutation test part 
     ** The code is exactly the same as before, except where labels[i] appeared (when computing sums of examples).
     ** There, instead of using labels[i], where i is an example index in groupStarts[g]:groupEnds[g], we have
     ** to use the label of the example <i> got permuted to. Hence the loops that would be groupStarts[g]:groupEnds[g]
     ** are now 0:(groupSize[g]-1), and i comes from the permutation
     **/
  
    /* create space for the permutation for each group and initialize */
    groupPermutations = (int **) mxMalloc( nGroups*sizeof(int *) );
    for (g=0; g<nGroups; g++) {
      groupPermutations[g] =  (int *) mxMalloc( groupSize[g]*sizeof(int) ); 
      idx = 0; for (i=groupStarts[g]; i<=groupEnds[g]; i++) { groupPermutations[g][idx] = i; idx++; }
    }

    /*
      for (g=0; g<nGroups; g++) {
      for(i=0; i<groupSize[g]; i++) { mexPrintf("%d ",groupPermutations[g][i]); }; mexPrintf("\n"); fflush(stdout);
      }
    */

    for (p=0; p<nPermutations; p++)
      {
	/** compute the sum of examples and examples^2 for each group/class
	 ** (store in array which has #groups x #classes rows, all classes for each group) */
  
	for (j=0; j<m; j++) { accuracyCount[j] = 0; };

	for (g=0; g<nGroups; g++)
	  {
	    /* do the permutation *within* the group (loops up to groupSize[g]-2) */
	    for (i=0; i<(groupSize[g]-1); i++) {
	      /* a number between i and groupSize[g]-1 (the last index in the vector) */
	      j = i + (int) ((groupSize[g]-i) * (rand() / (RAND_MAX + 1.0)));
	      itmp = groupPermutations[g][i];
	      groupPermutations[g][i] = groupPermutations[g][j];
	      groupPermutations[g][j] = itmp;
	    }

	    /* base of the matrix for this group (#classes x #voxels) */
	    x1ptr = sumX1 + g*(nClasses*m);
	    x2ptr = sumX2 + g*(nClasses*m);

	    /* initialize the respective matrix portion to 0 */
	    for (c=0; c<nClasses; c++) { mptr = x1ptr + c*m; for (j=0; j<m; j++) { mptr[j] = 0; } }
	    for (c=0; c<nClasses; c++) { mptr = x2ptr + c*m; for (j=0; j<m; j++) { mptr[j] = 0; } }
	    
	    /* compute the sums of examples, squares of examples and sums of squares of examples */
	    /* (this is where normal code and permutation code differ, both loop over i for contiguous
	       examples but the label access uses an index different from i, ig, which belongs to
	       another example in the same group) */
	  
	    /* mexPrintf("%d\t",g); */
	    idx = 0;
	    for (i=groupStarts[g]; i<=groupEnds[g]; i++) {

	      eptr  = examples + i*m;
	      e2ptr = examplesSquared + i*m;

	      ig = groupPermutations[g][idx];

	      /*mexPrintf("%d|%d\t",i,ig); fflush(stdout); */
	      /*mexPrintf("%d|%d ",labels[i],labels[ig]); fflush(stdout);*/
	    
	      /* add this example to the appropriate class tally */
	      mptr = x1ptr + labels[ig]*m; for (j=0; j<m; j++) { mptr[j] = mptr[j] + eptr[j]; }
	    
	      /* square it */
	      for (j=0; j<m; j++) { e2ptr[j] = eptr[j] * eptr[j]; }
	    
	      /* add the squared of the example to the appropriate class tally */
	      mptr = x2ptr + labels[ig]*m; for (j=0; j<m; j++) { mptr[j] = mptr[j] + e2ptr[j]; }

	      idx++;
	    } 
	    /* mexPrintf("\n");fflush(stdout); */
	  
	  }

	/*
	  for (g=0; g<nGroups; g++) {
	  for(i=0; i<groupSize[g]; i++) { mexPrintf("%d ",groupPermutations[g][i]); }; mexPrintf("\n"); fflush(stdout);
	  }
	*/

	/**
	 ** Main cross-validation loop, in every loop
	 ** - add up the x1 and x2 matrices for groups in the training set
	 ** - derive class mean and standard deviations from that
	 **/

	for (gtest=0; gtest<nGroups; gtest++)
	  {
	    /* mexPrintf("aqui! %d\n",gtest);fflush(stdout); */
	  
	    /* rezero the space for class means, stdv, sum X1, sum X2 and nc */      
	    for (c=0; c<nClasses; c++) { mptr = means + c*m; for (j=0; j<m; j++) { mptr[j] = 0; }; }
	    for (c=0; c<nClasses; c++) { mptr = meansSquared + c*m; for (j=0; j<m; j++) { mptr[j] = 0; }; }
	    mptr = vars;  for (j=0; j<m; j++) { mptr[j] = 0; }
	    for (c=0; c<nClasses; c++) { mptr = sumX1train + c*m; for (j=0; j<m; j++) { mptr[j] = 0; } }
	    for (c=0; c<nClasses; c++) { mptr = sumX2train + c*m; for (j=0; j<m; j++) { mptr[j] = 0; } }
	    for (c=0; c<nClasses; c++) { nc[c] = 0; }
	  
	    gidx = 0; /* counter of # groups being added */

	    /* sum X1 and X2 over the other groups, producing a #classes x #voxels array for each */
	    for (gtrain=0; gtrain<nGroups; gtrain++)
	      {
		if (gtrain != gtest)
		  {
		    /* add counts of each class in this group */
		    gptr = groupNclass + gtrain*nClasses; for (c=0; c<nClasses; c++) { nc[c] = nc[c] + gptr[c]; }
	      
		    /* start of the X1 and X2 arrays for this group */
		    x1ptr = sumX1 + gtrain*(nClasses*m);
		    x2ptr = sumX2 + gtrain*(nClasses*m);
	      
		    /* add up X1 */
		    for (c=0; c<nClasses; c++)
		      {
			mptr = sumX1train + c*m;
			for (j=0; j<m; j++) { mptr[j] = mptr[j] + x1ptr[j]; }
			x1ptr = x1ptr + m;
		      }
	      
		    /* add up X2 */
		    for (c=0; c<nClasses; c++)
		      {
			mptr = sumX2train + c*m;
			for (j=0; j<m; j++) { mptr[j] = mptr[j] + x2ptr[j]; }
			x2ptr = x2ptr + m;
		      }
	
		    gidx = gidx + 1;
		  }
	      } /* end of for tallying over training groups */
	    ntrain = 0; for (c=0; c<nClasses; c++) { ntrain = ntrain + nc[c]; }; ntrainmo = ntrain - 1;

	    /* compute means and meansSquared (one vector per class) */
	    /* compute variances (one vector for all classes) */

	    for (c=0; c<nClasses; c++)
	      {
		x1ptr = sumX1train + c*m;
		x2ptr = sumX2train + c*m;
		mptr  = means + c*m;
		sptr  = meansSquared + c*m;

		for (j=0; j<m; j++)
		  {
		    mptr[j] = x1ptr[j] / nc[c];
		    sptr[j] = mptr[j] * mptr[j]; /* compute the means squared to use below */
		    vars[j] = vars[j] + x2ptr[j] - 2*mptr[j]*x1ptr[j] + sptr[j] * nc[c];
		  }
	      }

	    for (j=0; j<m; j++) { vars[j] = vars[j] / ntrainmo; }

	    /** 
	     ** apply to test set (scroll over examples in the test group)
	     **/
      
	    /* create a buffer for each example where we store the current largest value,
	       class assignment is put in examplesCorrect
	       at the end check the label in examplesCorrect and set to 1 or 0
	    */

	    /*
	      idx = 0;
	      for (i=groupStarts[g]; i<=groupEnds[g]; i++) {
	  
	    */


	    idx = 0;
	    for (i=groupStarts[gtest]; i<=groupEnds[gtest]; i++)
	      {
		eptr  = examples        + i*m;
		e2ptr = examplesSquared + i*m;	  
		cptr  = examplesCorrect + i*m; /* while going over classes stores prediction, then binary 1/0 for right/wrong */
		mptr  = exampleLargestValue; for (j=0; j<m; j++) { mptr[j] = -1000000000; }

		/* index of example label after permutation */
		ig = groupPermutations[gtest][idx];

		/* main loop over voxels (within it, for every voxel, we loop over neighbours) */
		for (j=0; j<m; j++)
		  {
		    nn   = numberOfNeighbours[j];
		    nptr = voxelsToNeighbours + j*neighbourMax;
      
		    for (c=0; c<nClasses; c++)
		      {
			mptr  = means + c*m;
			m2ptr = meansSquared + c*m;

			/* compute p(x,x_neighbours|c), minus parts that are same for all classes */
			voxelValue = (2*eptr[j]*mptr[j] - m2ptr[j] - e2ptr[j])/vars[j];

			for (v=0; v<nn; v++)
			  {
			    neighbour = nptr[v];
			    voxelValue += (2*eptr[neighbour]*mptr[neighbour] - m2ptr[neighbour] - e2ptr[neighbour])/vars[neighbour];
			  }

			/* mexPrintf("%1.4f\t%1.4f\t",voxelValue,exampleLargestValue[j]); fflush(stdout); */

			/* if p(x,x_neighbours|c) is larger than what we've seen so far, update prediction */
			if (voxelValue > exampleLargestValue[j]) {
			  exampleLargestValue[j] = voxelValue;
			  cptr[j] = (double) c;
			  /* mexPrintf("better!\t");fflush(stdout); */
			} else {
			  /* mexPrintf("worse!\t");fflush(stdout); */
			}
			/* mexPrintf("%1.2f\n",cptr[j]); */
		  
		      } /* end of loop over classes */
      

		    /* transform to correct/incorrect (uses permutation labels rather than original ones)*/
		    if ( cptr[j] != ((double) labels[ig]) ) { cptr[j] = 0; } else { cptr[j] = 1; }
		    accuracyCount[j] = accuracyCount[j] + cptr[j];
	      
		  } /* end of loop over voxels */

		idx++;
	      } /* end of loop over test examples */

	  }  /* end of for over test groups */
      
	/* now that we have the accuracy count, compare that with accuracy on original labelling (stores count for now) */
      
	for (j=0; j<m; j++) { if (accuracyCount[j]>=accuracy[j]) { fraction[j]++; } }
      
      } /* end of loop over permutations */
    
    /* turn fraction into p-values */
  
    /* turn accuracy into actual accuracy  */
    for (j=0; j<m; j++) { fraction[j] = fraction[j]/(nPermutations+1); }

  }  /* end of if on permutation testing */

  /* convert accuracy from counts to accuracy values (needs to be done here because p-tests use counts) */
  for (j=0; j<m; j++) { accuracy[j] = accuracy[j]/n; }


  
  /* get rid of all the matrices obtained with mxMalloc (in the order the mxMallocs appear in the code) */
  /* mxFree(accuracyCount); */
  mxFree(labels);
  mxFree(labelsGroup);
  mxFree(groupStarts);
  mxFree(groupEnds);
  mxFree(groupSize);  
  mxFree(groupNclass);
  mxFree(numberOfNeighbours);
  mxFree(voxelsToNeighbours);
  mxFree(sumX1);
  mxFree(sumX2);
  mxFree(examplesSquared);
  mxFree(nc);
  mxFree(means);
  mxFree(vars);
  mxFree(meansSquared);
  mxFree(sumX1train);
  mxFree(sumX2train);
  /* mxFree(exampleLargestValue); */
  /* for (g=0; g<nGroups; g++) { mxFree(groupPermutations[g]); } */
  /* mxFree(groupPermutations); */

  return;


  /*  mexPrintf("here\n");fflush(stdout); */

  /* create other temporary matrices and C pointers to them
  sigma2_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  sigma2 = mxGetPr(sigma2_mx);

  logsigma_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  logsigma = mxGetPr(logsigma_mx);

  tmp1_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  tmp1 = mxGetPr(tmp1_mx);

  pEzM_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  pEzM = mxGetPr(pEzM_mx);
  pEzM2_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  pEzM2 = mxGetPr(pEzM2_mx);
  */


  /*
   * do the computation
   */
  
  /* debug: access  matrices coming in */
  /* NO PROBLEM HERE
  for (e=0; e<n; e++) {
    pEz = Ez + e*l;
    for (i=0; i<l; i++) { mexPrintf("%1.4f ",pEz[i]); }; mexPrintf("\n"); fflush(stdout);
  }
  mexPrintf("\n");

  for (e=0; e<n; e++) {
    pVz = Vz + e*l;
    for (i=0; i<l; i++) { mexPrintf("%1.4f ",pVz[i]); }; mexPrintf("\n"); fflush(stdout);
  }

  mexPrintf("\n");

  for (e=0; e<n; e++) {
    pX = X + e*m;
    for (j=0; j<4; j++) { mexPrintf("%1.4f ",pX[j]); };
    for (j=(m-4); j<m; j++) { mexPrintf("%1.4f ",pX[j]); }; 
    mexPrintf("\n"); fflush(stdout);
  }
  mexPrintf("\n");

  for (i=0; i<l; i++) {
    pM = M + i*m;
    for (j=0; j<4; j++) { mexPrintf("%1.4f ",pM[j]); };
    for (j=(m-4); j<m; j++) { mexPrintf("%1.4f ",pM[j]); }; 
    mexPrintf("\n");
  }
  mexPrintf("\n");

  return;
  */

  /* initialize the result

  *log_px = 0;  
  for (i=0; i<l; i++) {
    pL = log_px_pd + i*m;
    for (j=0; j<m; j++) {
      pL[j] = 0;
    }
  }
  */


  /* precompute things outside the loop

  for (j=0; j<m; j++) { sigma2[j]   = sigma[j]*sigma[j]; }
  for (j=0; j<m; j++) { logsigma[j] = log(sigma[j]); }
  fixed = 0.5 * log(2*PI);
  */

  for (e=0; e<n; e++) {

    /** initialize pointers for this example */
    /*
    pX  = X  + e*m;
    pEz = Ez + e*l;
    pVz = Vz + e*l;
    */

    /** 1) work on the log_px computation for this example */

    /* compute E[zM] and E[(zM)^2]
    for (j=0; j<m; j++) { pEzM[j] = 0;}
    for (j=0; j<m; j++) { pEzM2[j] = 0;}

    for (i=0; i<l; i++) {
      pM = M + i*m;
      for (j=0; j<m; j++) {
	pEzM[j] = pEzM[j] + pEz[i] * pM[j];
      }
    }

    for (i=0; i<l; i++) {
      pM = M + i*m;
      for (j=0; j<m; j++) {
	pEzM2[j] = pEzM2[j] + pVz[i]*(pM[j]*pM[j]);
      }
    }
    for (j=0; j<m; j++) {
      pEzM2[j] = pEzM2[j] + (pEzM[j]*pEzM[j]);
    }
    

    for (j=0; j<m; j++) {
      *log_px = *log_px + (pX[j]*pEzM[j] -0.5*pEzM2[j] -0.5*(pX[j]*pX[j]))/sigma2[j] -logsigma[j] -fixed;
    }
    */


    /*
    for (j=0; j<5; j++) { mexPrintf("%1.3f ",pEzM[j]); };  mexPrintf("\n"); fflush(stdout);
    for (j=0; j<5; j++) { mexPrintf("%1.3f ",pEzM2[j]); }; mexPrintf("\n"); fflush(stdout);
    mexPrintf("%f\n",*log_px);
    mexPrintf("aqui!\n");fflush(stdout); scanf("%s",buffer);
    */

    /* 2) work on the log_px_pd computation for this example */

    /* for (j=0; j<5; j++) { mexPrintf("%1.3f ",pX[j]); }; mexPrintf("\n"); fflush(stdout); */
    /* for (i=0; i<l; i++) { mexPrintf("%1.3f ",pEz[i]); }; mexPrintf("\n"); fflush(stdout); */

    /* precompute Ez*M  for this example */
    
    /*
    for (j=0; j<m; j++) { tmp1[j] = 0; }
    for (i=0; i<l; i++) {
      pM = M + i*m;
      for (j=0; j<m; j++) {	
	tmp1[j] = tmp1[j] + pEz[i] * pM[j];
	idx++;
      }
    }
    */

    /*
    for (j=0; j<5; j++) { mexPrintf("%1.4f ",tmp1[j]); }; mexPrintf("\n"); fflush(stdout);
    mexPrintf("aqui!\n");fflush(stdout); scanf("%s",buffer);
    */
    /* printf("aqui!\n");fflush(stdout); scanf("%s",buffer); */

  } /* end of loop over examples */

  /*
  mxDestroyArray(tmp1_mx);
  mxDestroyArray(sigma2_mx);
  mxDestroyArray(logsigma_mx);
  */
  /* mxFree(tmp1); */
}
