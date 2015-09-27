#include <stdlib.h>		// (to include exit(1)) 
# include <stdio.h>
# include <math.h>
# include "integral.h"
//# include "data_types.h"



/****************************************************************
 * done_integral --  checks whether an integral has been 
 *                   previously calculated 
 * store_integral -  stores the values of previously calculated 
 *                   integrals
 *
 * Inputs --    y     : allowed intensities
 *              k     : clump number
 *              value : value of integral to be stored
 *
 * 
 * Returns --   done_integral() returns the value of the integral if it has
 *              been previously calculated and 0 otherwise     
 *               
 * 
 * Credits --   Numerical recipies in C - various routines
 *
 *
 *
 ***************************************************************/
      


double done_integral(int *y, long k)
{
	int i,j;
	
	int same;
	
	// simple checks all the stored values of y[] for a given clump,
	// if it finds one that matches the input then it returns the corresponding
	// value, if it can't find any match it returns 0;
	
	// integral[k].count is the number of integrals it has stored
	// integral[k].N is the length of y[]
	// integral[k].y[j][i] is the ith element of the jth stored y[]
	
	
	for(j=1;j<=integral[k].count;j++)
	{
    same=1;
    for(i=1;i<=integral[k].N;i++)
		{
			if(y[i]!=integral[k].y[j][i]) 
			{
				same=0;
				break; 
			}      
		}
    if(same==1) 
		{
			return integral[k].value[j];
		}
	}
	
	return 0.0;
}


// this simply stores y[] and values in the appropriate place

void store_integral(int *y, double value, long k)
{
	int i;
	
	if(integral[k].count<MAXINT)		
  {
    integral[k].count++;
		
    integral[k].value[integral[k].count] = value;
		
    for(i=1;i<=integral[k].N;i++)
			integral[k].y[integral[k].count][i] = y[i];
  }
	else
	{
		// MAXENT is the maximum no of integrals which can be stored for each clump.
		// MAXENT = 111 and it should be possible that integral[k].count exceed that
		// value. But in case this somehow happens

		printf("\n\n*******************************************************\n");
		printf("Extreme situation has occurred in function store_integral()\n");
		printf("Please contact the author of this program Anders J. Markvardsen\n");
		printf("Thank you\n");
		printf("*******************************************************\n\n");

		exit(1);
	}	
}
