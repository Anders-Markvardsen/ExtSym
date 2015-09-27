#include <stdlib.h>		// (to include exit(1)) 
#include <stdio.h>
#include <math.h>
#include "../lib/nrutil.h"
#include "../lib/nr.h"
#include "read_in_parameters.h"  // input_file_format


#include "data_types.h"
#include "integral.h"


/****************************************************************
* get_data --  a function that gets the data from the output
*              of a Pawly refinement code and turns it into 
*              something useful. It also has a secondary role of
*              initiating various arrays and matrices
*
* Inputs --    The file ~/space_data/*
*
* Returns --   returns the data in a usable structure 
*               : struct store_data data
*
* Comment --	The input file must be of a sequence of line with format
*
*				  h(int) k(int) l(int) Int(double) Sigma(double) Counter(long) ...
*				  CorToNext(double) CorToSecond(double) etc.
*				
*				where CorToNext is the correlation to the next reflection etc.
*
*				A double is here identified when the counter stays the same for the
*				next reflection. 				
*
***************************************************************/


store_data *data;			// where input data are kept

store_new_data *new_data;	// where symmetry modified data are kept

st_int *integral;		// st_int defined in integral.h 
       

// some haking : changed from # define PAWLEY_COR_ELEMENTS 15 to static long ... (22/2-2001)
// (also # undef PAWLEY_COR_ELEMENTS was commented out)

static long PAWLEY_COR_ELEMENTS;     // this is a constraint which orinates from the 
									// special format of the pawley input file.
									// Used for dynamical allocation of arrays.

// comment out for new
//extern struct store_data data[];             // (evidence.c)
//extern struct store_new_data new_data[];     // (evidence.c)

int total_clump;                             // number of clumps
// GLOBAL variable defined in this file 
//(used in contract.c and evidence.c)


typedef enum {FALSE, TRUE} my_boolean;		// define new type my_boolean


void get_information_infile(long *number_of_lines_ptr, long *num_extra_lines_ptr, 
							long *highest_multiplet_ptr, char infilename[], 
							long *input_correlation_element_ptr);

void convert_hkl_weight_to_hkl();
void save_data0_to_hkl_file();

void allocate_data (long clump_index, long position_first_peak, int clumpsize,
										long highest_multiplet);

int check_positive_definite(long position_first_peak, int clumpsize);

void print_out(long clump_index);

FILE *Fopen(char *fname, char *fmode); // from read_in_parameters.c 

//---------------------------------------------------------------------------------------
// get_data() first assigns pointer addresses to data and newdata. Next 
// get_data() loads the input data into data[0] first and then split data[0] into 
// smaller blocks data[1], data[2], ... , data[total_blocks]. 
// The necessary space is both allocated for data and new_data.
// Also the total number of blocks(clumps) are returned in total_blocks
//---------------------------------------------------------------------------------------
void get_data (long *total_blocks, char infilename[],
							 double sigma_multiply, double correlation_cutoff)
{
	long i,j,n;			// simple matrix indices
		
	int counter;			// the counter in the data file
	int last_counter=0;           // the last counter in the data file
	
	long num_refined_reflections;	// number of reflections refined in the Pawley 
									// refinement

	long num_extra_lines;				// number of extra lines from multiplets

	long highest_multiplet;    // highest multiplet in input data file
														 // highest_multiplet 0 = only singlets
														 // highest_multiplet 1 doublet and so on

	long number_of_lines;			// equal to num_refined_reflections + num_extra_lines


	yesno end_of_clump;					// flag to indicate when to end adding to 
															// correlation block

	int clumpsize;						// size of a correlation block

	long clump_index;					// to index the different clumps


#ifdef SAVE_PEAK_IGNORED
	FILE *peak_ignored;			// fp to peak_ignored.asc
  char filename[50];  // to add OUTPUT_DIRECTORY to path
#endif


	FILE *data_file;		             // where the identification of the data file is kept

	
	// checking and getting number of refined reflections etc from infile
	
	get_information_infile(&number_of_lines, &num_extra_lines, &highest_multiplet, 
		infilename, &PAWLEY_COR_ELEMENTS);
		

	// open infile which will be read below 
	
	data_file = Fopen(infilename, "r"); 

#ifdef SAVE_PEAK_IGNORED
	// open the output file peak_ignored.asc which will output all ignored peaks

  sprintf(filename, "%speak_ignored.asc", OUTPUT_DIRECTORY);
	peak_ignored = Fopen(filename, "w");
#endif

	// number of refined reflections

	num_refined_reflections = number_of_lines - num_extra_lines;
	
	
	// in the structure store_data there are pointers that must be 
	// initiallised. This is done here for data[0], which is used as a 
	// area. This is transfered to the correct data[k] in allocate_data().

	// this is a bit of an overkill - I may improve on this later 
	// The +1 below is because we index the blocks
	// from one rather than zero and upwards
	data = (store_data *) malloc((num_refined_reflections+1)*sizeof(store_data));
	new_data = (store_new_data *) 
		malloc((num_refined_reflections+1)*sizeof(store_new_data));
	integral = (st_int *) malloc((num_refined_reflections+1)*sizeof(st_int));
	
	
	data[0].N = num_refined_reflections;   
	
	// the +highest_multiplet extra element allocated for data[0].cov_matrix, data[0].ave_intensity,
	// data[0].h, data[0].k and data[0].l is required when reading in the inputfile 
	// for the special case where the inputfile ends with a multiplet. For that case
	// data are read into data[0].h[num_refined_reflections+highest_multiplet] etc. 
	// Not good programming

	data[0].cov_matrix    = dmatrix(1, num_refined_reflections+highest_multiplet, 1, PAWLEY_COR_ELEMENTS);	
	
	data[0].ave_intensity = dvector(1, num_refined_reflections+highest_multiplet);
	
	data[0].h = ivector(1, num_refined_reflections+highest_multiplet);	// purely for displaying purposes
	data[0].k = ivector(1, num_refined_reflections+highest_multiplet);
	data[0].l = ivector(1, num_refined_reflections+highest_multiplet);
	
	data[0].multiplet = ivector(1, num_refined_reflections);
	
	// the +1 in highest_multiplet+1 is included to make exist when highest_multiplet = 0
	// Not good programming: it would make more sense only to allocate space when
	// highest_multiplet > 0

	data[0].h_multiplet = imatrix(1, num_refined_reflections, 1, highest_multiplet+1);	// purely for displaying purposes
	data[0].k_multiplet = imatrix(1, num_refined_reflections, 1, highest_multiplet+1);
	data[0].l_multiplet = imatrix(1, num_refined_reflections, 1, highest_multiplet+1);
	
	
	
	// turn 'off' all multiplet flags and initialise the sparce multiplet matrices
	
	for(n=1;n<=num_refined_reflections;n++)
	{
		data[0].multiplet[n] = 0;
		for(i=1; i <= highest_multiplet; i++)
		{
			data[0].h_multiplet[n][i] = 0;
			data[0].k_multiplet[n][i] = 0;
			data[0].l_multiplet[n][i] = 0;
		}
	}
	

	// initialize last_counter to dummy number

	last_counter = -1;


	// n here played the role of a peak index (for the refined peaks)

	n = 1;   // pointing at first peak in infilename


	// read in infilename into the structure data[0]
	// Notice when infilename ends with a doublet data are read in to the 
	// num_refined_reflections+highest_multiplet. Not good programming

	for(i = 1; i <= number_of_lines; i++)
  {
		// read in h,k,l
		
		fscanf (data_file, "%d%d%d", &data[0].h[n], &data[0].k[n], &data[0].l[n]);
		
		
		// read in the likelihood expectation value of intensity i
		
		fscanf (data_file, "%le", &data[0].ave_intensity[n]);
		
		
		// read in the ith diagonal element of the covariance matrix
		// which is at this stage read in as a standard deviation
		
		fscanf (data_file, "%le", &data[0].cov_matrix[n][1]);


		// read in the data counter
		
		fscanf (data_file, "%d", &counter);
		
		
		// if the counter is unchanged then we have a multiplet
		
		if(counter==last_counter)
		{
			// step the refined reflection index one back 
			--n;  
			
			
			// increase multiplet by one
			
			(data[0].multiplet[n])++;      
			
			
			// put h,k,l values in the doublet storage arrays
			
			data[0].h_multiplet[n][data[0].multiplet[n]] = data[0].h[n+1];     
			data[0].k_multiplet[n][data[0].multiplet[n]] = data[0].k[n+1];     
			data[0].l_multiplet[n][data[0].multiplet[n]] = data[0].l[n+1]; 
			
			
			// since we have a multiplet there is no need to store its correlation 
			// properties again because it will be the same as it twin intentsity. 
			// However, since we may incounter an input with the format
			//     h k l ..... 32 100 4 0
			//     h k l ..... 32 4   0 0
			// rather then
			//     h k l ..... 32 4   0 0
			//     h k l ..... 32 4   0 0
			// then to take both these cases input account overwrite the line of
			// correlation matrix elements
			
			for(j=2;j<=PAWLEY_COR_ELEMENTS;j++)
      {
				fscanf(data_file,"%le",&data[0].cov_matrix[n][j]);

        // In principle the input file should not contain -100 or 100 correlation 
        // elements - since if two peaks are 100% correlated than better to treat
        // these as a one refineable. Further to avoid rule out the possible that
        // a matrix inversion struggle then do the stuff below

        if (data[0].cov_matrix[n][j] < -100.0001 || data[0].cov_matrix[n][j] > 100.0001)
        {
          printf("\n\n***********************************************************\n\n");
				  printf("ERROR IN INPUT FILE:\n");	
				  printf("  Correlation element for reflection %d %d %d\n", 
            data[0].h[n], data[0].k[n], data[0].l[n]);
				  printf("  is larger than 100 percent, which is clearly nonsense\n.");
				  printf("***********************************************************\n\n");

				  exit(1);
        }

        if (data[0].cov_matrix[n][j] < -99.001)
          data[0].cov_matrix[n][j] = -99.0;

        if (data[0].cov_matrix[n][j] > 99.001)
          data[0].cov_matrix[n][j] = 99.0;

      }

			// increase the refined peak index by one so we are ready to insert more data
			// into data[0]

			++n;


			// read in next line
			
			continue;
		}
		
		
		// update last counter
		
		last_counter = counter;
		
		
		// read in the remaining % covariance components for this intensity
		
		for(j=2;j<=PAWLEY_COR_ELEMENTS;j++)
    {
			fscanf(data_file,"%le",&data[0].cov_matrix[n][j]);

      // In principle the input file should not contain -100 or 100 correlation 
      // elements - since if two peaks are 100% correlated than better to treat
      // these as a one refineable. Further to avoid rule out the possible that
      // a matrix inversion struggle then do the stuff below

      if (data[0].cov_matrix[n][j] < -100.0001 || data[0].cov_matrix[n][j] > 100.0001)
      {
        printf("\n\n***********************************************************\n\n");
				printf("ERROR IN INPUT FILE:\n");	
				printf("  Correlation element for reflection %d %d %d\n", 
          data[0].h[n], data[0].k[n], data[0].l[n]);
				printf("  is larger than 100 percent, which is clearly nonsense\n.");
				printf("***********************************************************\n\n");

				exit(1);
      }

        //if (data[0].cov_matrix[n][j] < -90.001 + 5*(j-2))
        //  data[0].cov_matrix[n][j] = -90.0 + 5*(j-2);

        //if (data[0].cov_matrix[n][j] > 90.001 - 5*(j-2))
        //  data[0].cov_matrix[n][j] = 90.0 - 5*(j-2);

      if (data[0].cov_matrix[n][j] < -99.001)
        data[0].cov_matrix[n][j] = -99.0;

      if (data[0].cov_matrix[n][j] > 99.001)
        data[0].cov_matrix[n][j] = 99.0;
    }

		// increase the refined peak index by one

		++n;
  }


	// double check that infilename has been stored correctly in data[0]
	
	if (n-1 != num_refined_reflections)
	{
		printf("Error in reading %s into the structor data[0] in get_data.c\n",
			infilename);
		exit(1);
	}


  // For debugging you might what to convert an .hkl_weight into an .hkl. 
  // However notice this function will only work for reasonable behaving
  // .hkl_weight inputs

//  if (input_file_format == 'i')
//	{
//    convert_hkl_weight_to_hkl();
//  }


  // Also for debugging purpose you might want to save a data[0] as
  // a .hkl file named 'dummy.hkl'

//  save_data0_to_hkl_file();


  // Down-scale or up-scale the estimated errors on the intensities 
  // according to the sigma_multiply parameter in advanced.asc.
  //
  // Notice only the diagonal elements need to be multiplied (or divided)
  // by sigma_multiply since all the off-diagonal elemnents are read in
  // a normalised matrix elements

  if (input_file_format == 'i')
    for (i = 1; i <= data[0].N; i++)
		  data[0].cov_matrix[i][1] /= sigma_multiply;
  else
    for (i = 1; i <= data[0].N; i++)
		  data[0].cov_matrix[i][1] *= sigma_multiply;


	// which clump
	
	clump_index = 0;                       
	

	// cutting the correlation matrix into independent smaller correlation blocks.
	// i here plays the role of a refined peak index. 

	for(i = 1;i<=num_refined_reflections;)            
    {
		// initiate clumpsize

		clumpsize = 0;

		
		// Finding a correlated block in the correlation matrix. 
		// First add one to clumpsize and check if the correlation to the next peak is bigger than correlation_cutoff. 
		// Then if this is the case add one to clumpsize and then check if peak i is 
		// correlated with peak i+2 and if peak i+1 is correlated with i+2 and so on.
		// This way a triangular region of the correlation matrix is checked for any 
	    // elements which have a correlation higher than correlation_cutoff
		// At the end of the do loop we have the size of the correlation block, 
		// clumpsize, and of course, i, the first peak of the correlation block. 

		do 
		{
			// new member of clump
			++clumpsize;

			// make at check of input data file. Error if there is detected correlation to a 
			// peak with is not there

			if (clumpsize + i > num_refined_reflections + 1)
			{
				printf("\n\n***********************************************************\n\n");
				printf("ERROR IN INPUT FILE:\n");	
				printf("  Please double check your input file, since one or more correlation\n");
				printf("  matrix elements in the input data file predict correlation to\n");
				printf("  reflection number %d, which is not present in the input file!\n\n", num_refined_reflections + 1);
				printf("***********************************************************\n\n");

				exit(1);
			}

			// I do want to end clump

			end_of_clump = YES;


			// however not if there is correlation to the i+clump_size reflection

			if (clumpsize >= PAWLEY_COR_ELEMENTS)
			{
				for(j=clumpsize+2-PAWLEY_COR_ELEMENTS;j<=clumpsize;j++)
					if(fabs(data[0].cov_matrix[j+i-1][2+clumpsize-j]) >
						correlation_cutoff) 
						end_of_clump = NO;
			}
			else
			{
				for(j=1;j<=clumpsize;j++)
					if(fabs(data[0].cov_matrix[j+i-1][2+clumpsize-j]) >
						correlation_cutoff) 
						end_of_clump = NO;
			}

		} while (end_of_clump == NO);


		// increase clump-index by one

		++clump_index;


		// allocate the data to the correct storage arrays ONLY if the covariance is
		// positive definite. This extra criterion is included because the input data
		// might be crap, and so crap that the correlation matrix is no longer positive
		// definite
					
		if (check_positive_definite(i, clumpsize) == 0)
		{
			// rubbish data and in principle I should display and error message and stop the
			// problem, but for now instead just ignore such data by decrease clump_index
			// by one

			--clump_index;

#ifdef SAVE_PEAK_IGNORED
			fprintf(peak_ignored, "\n\n----------- WARNING: A NON-POSITIVE DEFINITE BLOCK ---------\n\n");
			fprintf(peak_ignored, "Intensity from %ld to %ld -- now ignored.\n\n", i, i+clumpsize-1);
#endif
		}
		else
		{
			allocate_data(clump_index, i, clumpsize, highest_multiplet);
		
			// diagnostics for this clump
					
			//print_out(clump_index);
		}						
		

		// let i point to first peak following the most recent registered clump

		i += clumpsize;
	}

#ifdef SAVE_PEAK_IGNORED
	// finish ignoring peaks for now

	fclose(peak_ignored);
#endif

	// a double check here. Because the last peak of infilename must be uncorrelated with
	// its next peak (since there is no next peak) then the index i must at this point 
	// have the value num_refined_reflections+1.

	if (i != num_refined_reflections + 1) 
	{
		printf("\n\n***********************************************************\n\n");
				printf("ERROR:\n");
				printf("  Error in block diagonalising the correlation matrix of\n");
				printf("  %s.\n", infilename);
				printf("  Please contact author of program Anders J. Markvardsen,\n");
				printf("  since I would be interested in knowing what caused this error\n");
				printf("  to be displayed. Thank you in advance and I apologize for the\n");
				printf("  inconvenience.\n\n");
				printf("***********************************************************\n\n");

				exit(1);
	}


	// global variable used in contract.c and evidence.c 
	
	total_clump = clump_index;	

	
	// return total_blocks -- equal to number of clumps(blocks)

	*total_blocks = total_clump;
	


	// force all isolated negative peaks to zero

	for(i = 1;i<=total_clump;i++)            
    {
		if(data[i].N == 1 && data[i].ave_intensity[1] < 0.0)
		{
			data[i].ave_intensity[1] = 0.0;
		}
	}


	
	// free up the data structure which is used to originally load in the input data
	// - before the data are distributed into clumps

	free_dmatrix(data[0].cov_matrix,   1, num_refined_reflections+highest_multiplet, 1, PAWLEY_COR_ELEMENTS);	
	
	free_dvector(data[0].ave_intensity,1, num_refined_reflections+highest_multiplet);
	
	free_ivector(data[0].h, 1, num_refined_reflections+highest_multiplet);
	free_ivector(data[0].k, 1, num_refined_reflections+highest_multiplet);
	free_ivector(data[0].l, 1, num_refined_reflections+highest_multiplet);
	
	free_ivector(data[0].multiplet,1, num_refined_reflections);
	
	free_imatrix(data[0].h_multiplet, 1, num_refined_reflections, 1, highest_multiplet+1);
	free_imatrix(data[0].k_multiplet, 1, num_refined_reflections, 1, highest_multiplet+1);
	free_imatrix(data[0].l_multiplet, 1, num_refined_reflections, 1, highest_multiplet+1);
}


//--------------------------------------------------------------------
// temporary function for checking if correlation matrix is positive definite
//     returns   0 -- not positive definite
//               1 -- then positive definite
//----------------------------------------------------------------------

int check_positive_definite(long position_first_peak, int clumpsize)
{
	int i, j;
	double **v,*d;
	int nrot;

	int return_value;
	
	double **Cov;

	Cov = dmatrix(1, clumpsize, 1, clumpsize);
	

	// reading in diagonal elements first without squaring them
	
	for(i=1;i<=clumpsize;i++) { 
		Cov[i][i] = data[0].cov_matrix[position_first_peak-1+i][1];
	}
	
	
	// filling in the off-diagonal elements

	for(i=1;i<=clumpsize;i++)
	{
		for(j=i+1;j<=clumpsize;j++)
		{
			if (j-1+1 <= PAWLEY_COR_ELEMENTS)
			{
				Cov[i][j] =  
					data[0].cov_matrix[position_first_peak-1+i][j-i+1]
					* Cov[i][i]
					* Cov[j][j] / 100.0;
			}
			else
			{
				// this is because I only read in correlation matrix element up to the 15th next intensity

				Cov[i][j] = 0.0;
			}

			Cov[j][i] =  Cov[i][j];
		}
	}
		

	// square the diagonal elements to get variances

	for (i = 1; i <= clumpsize; i++)
		Cov[i][i] = Cov[i][i] * Cov[i][i];


	d = dvector(1,clumpsize);
	v = dmatrix(1,clumpsize,1,clumpsize);

	// jacobi() returns the eigenvalues of invin in d[1...ndim] and its eigenvectors in 
	// v[1...ndim][1...ndim]

	jacobi(Cov,clumpsize,d,v,&nrot);


	// check if any of the eigenvectors are negative

	return_value = 1;
	for(i=1; i<=clumpsize; i++)
	{
		if (d[i] <= 0.0)
		{
			return_value = 0;
		}
	}

	
	free_dvector(d,1,clumpsize);	
	free_dmatrix(v,1,clumpsize,1,clumpsize);

	free_dmatrix(Cov, 1, clumpsize, 1, clumpsize);

	return return_value;
}







// ********************************************************************
// this allocates the data, but also initiates various storage arrays 
// used by other functions
// ********************************************************************

void allocate_data (long clump_index, long position_first_peak, int clumpsize,
										long highest_multiplet)
{
	int i,j;                    // simple counters
	
	double **dummy_matrix;               // local storage matrix
	
	
	
	// initiate   struct store_data data[]
	
	data[clump_index].reference_point = position_first_peak;	// for displaying purposes

	data[clump_index].N                 = clumpsize;

	data[clump_index].monte_ave     = dvector(1,data[clump_index].N); 

	data[clump_index].extra_mean_values = dvector(1,data[clump_index].N);
	
	data[clump_index].ave_intensity     = dvector(1, data[clump_index].N);
	
	data[clump_index].inv_cov           = dmatrix(1, data[clump_index].N, 1, data[clump_index].N);
	data[clump_index].cov_matrix				= dmatrix(1, data[clump_index].N, 1, data[clump_index].N);
	
	// purely for displaying purposes
	data[clump_index].h                 = ivector(1, data[clump_index].N);
	data[clump_index].k                 = ivector(1, data[clump_index].N);
	data[clump_index].l                 = ivector(1, data[clump_index].N);
	
	data[clump_index].multiplet           = ivector(1, data[clump_index].N);	// = 1 doublet
																			// = 0 singlet	= 2 triplet and so on

	// the +1 in highest_multiplet+1 not good programming see above allocation of
	// data[0].h_multiplet

	data[clump_index].h_multiplet                = imatrix(1, data[clump_index].N, 1, highest_multiplet+1);
	data[clump_index].k_multiplet                = imatrix(1, data[clump_index].N, 1, highest_multiplet+1);
	data[clump_index].l_multiplet                = imatrix(1, data[clump_index].N, 1, highest_multiplet+1);
	
	
	
	
	
	// initiate   struct store_new_data new_data[]
	
	new_data[clump_index].monte_ave     = dvector(1,data[clump_index].N); 

	new_data[clump_index].extra_mean_values = dvector(1,data[clump_index].N);
	
	new_data[clump_index].cov_matrix    = dmatrix(1, data[clump_index].N, 1, data[clump_index].N);
	
	new_data[clump_index].ave_intensity = dvector(1, data[clump_index].N);
	
	new_data[clump_index].inv_cov       = dmatrix(1, data[clump_index].N, 1, data[clump_index].N);
	
	
	new_data[clump_index].allowed       = ivector(1, data[clump_index].N);		// = 0 forbiden
															// = 1 singlet       
															// = 2 doublet   
	
	
	new_data[clump_index].multiplet       = ivector(1, data[clump_index].N);		// = 0 singlet
															// = 1 doublet

	// purely for displaying purposes
	new_data[clump_index].h                 = ivector(1, data[clump_index].N);
	new_data[clump_index].k                 = ivector(1, data[clump_index].N);
	new_data[clump_index].l                 = ivector(1, data[clump_index].N);
															// = 0 singlet	
	new_data[clump_index].h_multiplet                = imatrix(1, data[clump_index].N, 1, highest_multiplet+1);
	new_data[clump_index].k_multiplet                = imatrix(1, data[clump_index].N, 1, highest_multiplet+1);
	new_data[clump_index].l_multiplet                = imatrix(1, data[clump_index].N, 1, highest_multiplet+1);
	
	
	// intitiate   struct st_int integral[]

  integral[clump_index].N      = clumpsize; 

  integral[clump_index].count  = 0;

  integral[clump_index].y      = imatrix(1, MAXINT, 1, integral[clump_index].N); // MAXINT is defined in 
																																								 // integral.h
	integral[clump_index].value = dvector(1, MAXINT);


	// because of the form of the input data, the covariance matrix is in 
	// a tricky form. This bit changes the upper off-diagonal elements from 
	// %-values to true values and also copies them into the lower diagonal.
	// It then squares the diagonal elements to make sd's into variances
	
	// to be used in get inverse

  dummy_matrix = dmatrix(1, data[clump_index].N, 1, data[clump_index].N);
	

	// reading in diagonal elements first without squaring them
	
	for(i=1;i<=data[clump_index].N;i++) { 
		data[clump_index].cov_matrix[i][i] = data[0].cov_matrix[position_first_peak-1+i][1];
	}
	
	
	// filling in the off-diagonal elements

	for(i=1;i<=data[clump_index].N;i++)
	{
		for(j=i+1;j<=data[clump_index].N;j++)
		{
			if (j-1+1 <= PAWLEY_COR_ELEMENTS)
			{
				data[clump_index].cov_matrix[i][j] =  
					data[0].cov_matrix[position_first_peak-1+i][j-i+1]
					* data[clump_index].cov_matrix[i][i]
					* data[clump_index].cov_matrix[j][j] / 100;
			}
			else
			{
				// this is because I only read in correlation matrix element up to the 15th next intensity

				data[clump_index].cov_matrix[i][j] =  0.0;				
			}

			data[clump_index].cov_matrix[j][i] =  data[clump_index].cov_matrix[i][j];
		}
	}
		

	// square the diagonal elements to get variances

	for (i = 1; i <= clumpsize; i++)
		data[clump_index].cov_matrix[i][i] = 
		data[clump_index].cov_matrix[i][i] * data[clump_index].cov_matrix[i][i];


	// Cov is used in inverse_jac

	for (i = 1; i <= clumpsize; i++)
		for (j = 1; j<= clumpsize; j++)
			dummy_matrix[i][j] = data[clump_index].cov_matrix[i][j];


	// assign values the clump data arrays

	for(i=1;i<=data[clump_index].N;i++)
	{
		data[clump_index].ave_intensity[i] = data[0].ave_intensity[position_first_peak-1+i];
		
		data[clump_index].h[i]       = data[0].h[position_first_peak-1+i];
		data[clump_index].k[i]       = data[0].k[position_first_peak-1+i];
		data[clump_index].l[i]       = data[0].l[position_first_peak-1+i];
		
		data[clump_index].multiplet[i] = data[0].multiplet[position_first_peak-1+i];
		
		// as data[0].h2, data[clump_index].h2 will also be a sparce array etc
		// the +1 in j+1 is not good programming see above allocation of
		// data[0].h_multiplet

		for (j = 0; j <= highest_multiplet; j++)
		{
			data[clump_index].h_multiplet[i][j+1] = data[0].h_multiplet[position_first_peak-1+i][j+1];	 
			data[clump_index].k_multiplet[i][j+1] = data[0].k_multiplet[position_first_peak-1+i][j+1];	
			data[clump_index].l_multiplet[i][j+1] = data[0].l_multiplet[position_first_peak-1+i][j+1]; 
		}
	}
	
	
	// the matrix dummy_matrix[][] is now inverted to get data[clump_index].inv_cov[][] or 
	// data[clump_index].cov_matrix[][].
	// Notice dummy_matrix is destroyed when inverse_jac is called

	inverse_jac(dummy_matrix, data[clump_index].inv_cov, 
		&data[clump_index].log_det_cov, data[clump_index].N);


	// in case input_file_format = inverse correlation matrix

	if (input_file_format == 'i')
	{
		for(i=1;i<=data[clump_index].N;i++)
		{
			for(j=1;j<=data[clump_index].N;j++)
			{
				dummy_matrix[i][j] = data[clump_index].cov_matrix[i][j];
				data[clump_index].cov_matrix[i][j] = data[clump_index].inv_cov[i][j];
				data[clump_index].inv_cov[i][j] = dummy_matrix[i][j];
			}
		}
		data[clump_index].log_det_cov = - data[clump_index].log_det_cov;
	}




	// finished used local variable Cov	

	 free_dmatrix(dummy_matrix, 1, data[clump_index].N, 1, data[clump_index].N);
}



// ***********************************************************************
// This function checks if the file infilename has the format
//		h(int) k(int) l(int) Int(double) Sigma(double) Counter(long) ...
//		CorToNext(double) CorToSecond(double) etc. 	
// And returns the number of lines and the number of doubles in infilename through
// its first two arguments
// ***********************************************************************

void get_information_infile(long *number_of_lines_ptr, long *num_extra_lines_ptr, 
							long *highest_multiplet_ptr, char infilename[], 
							long *input_correlation_element_ptr)
{
	long i, j;

	long highest_multiplet_dummy;   // used in determining highest_multiplet

	FILE *infile_fp;		// where the file infilename is kept
	
    int dummy_h, dummy_k, dummy_l;		// used to read in numbers from the
	double dummy_int, dummy_sigma;		// file infilename
	long dummy_counter;					//	
		
	int previous_counter;				// used to check for doubles
	
	char line[5000];						// to read in first line of infilename
	
	
	// open infilename
	
	infile_fp = Fopen(infilename,"r");

	
	// read in first line	
	
	if (fgets(line, 5000, infile_fp) == NULL) 
	{
		printf("\nError reading first line of %s\n", infilename);

		exit(1);
	}

	
	// counting the number of numbers in top line of input file

	i = 0; // plays the role of stepping through the string line
	j = 0; // plays the role of counting the number of numbers in line
	do
	{
		i++;
		// new element if line[i] is either a space or tab and the element before is not

		if ((line[i] == ' ' || line[i] == '\t') && (line[i-1] != ' ' && line[i-1] != '\t')) 
		{
			j++;
		}
		if (line[i] == '\0') break;
	} while (line[i] != '\n');
	
	if (line[i-1] != ' ' && line[i-1] != '\t') j++; // check if a number at the end of line

	*input_correlation_element_ptr = j - 5;  // this is the no of cor elem read in
																					 // including the diagonal element


	// checks that the format of the input file is not completely crazy

	if ( 6 == sscanf(line, "%d%d%d%lf%lf%ld", &dummy_h, &dummy_k, &dummy_l, 
		&dummy_int, &dummy_sigma, &dummy_counter) )
	{

	}
	else
	{
		printf("The input format of %s is invalid\n", infilename);
		printf("Please see manual.doc explaining the required format of the input file\n");

		exit(1);
	}
	
	
	// move the file pointer associated with infile_fp back to the beginning of the file
	
	fseek(infile_fp, 0L, SEEK_SET);
	
	
	// initiate previous_counter to dummy number
	
	previous_counter = -1;
	
	
	// initiate number of lines and multiplet
	
	*number_of_lines_ptr = 0;
	*num_extra_lines_ptr = 0;
	*highest_multiplet_ptr = 0;
	highest_multiplet_dummy = 0;
	
	// find the number of lines etc
	
	while ( fgets(line, 5000, infile_fp) != NULL )
	{
		
		i = 0; // plays the role of stepping through the string line
		j = 0; // plays the role of counting the number of numbers in line
		do
		{
			i++;
			if ((line[i] == ' ' || line[i] == '\t') && (line[i-1] != ' ' && line[i-1] != '\t'))
			{
				j++;
			}
			if (line[i] == '\0') break;
		} while (line[i] != '\n');
		if (line[i-1] != ' ' && line[i-1] != '\t') j++; 
		
		if (j == 0)
		{
			break;
		}
		if(j-5 != *input_correlation_element_ptr)
		{
			printf("\n\n***********************************************************\n\n");
			printf("ERROR IN INPUT FILE:\n");
			printf("  Line %d of input file does not contain the same number of\n", *number_of_lines_ptr+1);
			printf("  input data as the line above\n\n");
			printf("***********************************************************\n\n");

			exit(1);
			
		}
		
		
		// increase number of lines by one
		
		++(*number_of_lines_ptr);
		
		
		// check for the presence of a multiplet
		
		sscanf(line, "%d%d%d%lf%lf%ld", &dummy_h, &dummy_k, &dummy_l, 
			&dummy_int, &dummy_sigma, &dummy_counter);
		
		if (previous_counter == dummy_counter)
		{
			// increase number of doublets by one
			
			++(highest_multiplet_dummy);
			++(*num_extra_lines_ptr);
		}
		else
		{
			if (highest_multiplet_dummy > *highest_multiplet_ptr)
				*highest_multiplet_ptr = highest_multiplet_dummy;

			// put highest_multiplet_dummy back to zero again

			highest_multiplet_dummy = 0;

			// set previous_counter to current counter
		
			previous_counter = dummy_counter;
		}
	}
    

	// close file
	
	fclose(infile_fp);
}

// ********************************
// Used when a .hkl_weight file is read into the program. 
// Then the populated data[0].cov_matrix contains normalized
// inverse weight matrix elements and the purpose of this 
// function is to convert these into normalised covariance
// matrix element of the format as obtained when reading a hkl 
// file
//
// 
// ********************************
void convert_hkl_weight_to_hkl()
{
  double **weight_matrix; 
  double **cov_matrix;

  int i, j, N;

  double dummy;

//	double **v,*d;  // to get eigenvalues
//	int nrot;       // to get eigenvalues
// int no_negative; // to judge eigenvalues

//  for debugging

  	FILE *fp;
  char filename[50];  // used to add OUTPUT_DIRECTORY path

  sprintf(filename, "%sfisse_weight.hkl", OUTPUT_DIRECTORY);

  fp = Fopen(filename, "w");

  // end for debugging





  N = data[0].N;

  weight_matrix = dmatrix(1, N, 1, N);
  cov_matrix = dmatrix(1, N, 1, N);


  // initialising weight_matrix

  for (i = 1; i <= N; i++)
  {
    for (j = i+1; j <= N; j++)
    {
      weight_matrix[i][j] = 0.0;  
    }
  }

  // expand data[0].cov_matrix into a proper weight matrix

  // first read in diagonal elements without squaring these

  for (i = 1; i <= N; i++)
  {
    weight_matrix[i][i] = data[0].cov_matrix[i][1];
  }

  // filling in the off diagonal elements

  for (i = 1; i <= N; i++)
  {
    for (j = i+1; j <= N; j++)
    {  
      if (j-i+1 <= PAWLEY_COR_ELEMENTS)
      { 
        weight_matrix[i][j] = data[0].cov_matrix[i][j-i+1]
        * weight_matrix[i][i] * weight_matrix[j][j] / 100.0;

        weight_matrix[j][i] = weight_matrix[i][j];
      }
    }
  }


  // check to see if weight_matrix is positive definite
  // jacobi() returns the eigenvalues in d[1...ndim] and its eigenvectors in 
	// v[1...ndim][1...ndim]

 // d = dvector(1,N);
	//v = dmatrix(1,N,1,N);

	//jacobi(weight_matrix,N,d,v,&nrot);


	//// check if any of the eigenvectors are negative

	//no_negative = 1;
	//for(i=1; i<=N; i++)
	//{
	//	if (d[i] <= 0.0)
	//	{
	//		no_negative = 0;
	//	}
 //   //printf("%g\n", d[i]);
	//}

	//free_dvector(d,1,N);	
	//free_dmatrix(v,1,N,1,N);


  // squaring diagonal elements

  for (i = 1; i <= N; i++)
  {
    weight_matrix[i][i] = weight_matrix[i][i] * weight_matrix[i][i];
  }

  
  inverse_jac(weight_matrix, cov_matrix, &dummy, N);

    // for debuggin

  for (i = 1; i <= N-15; i++)
  {
    for (j = i; j <= i+14; j++)
    {
      fprintf(fp, "%5g ", cov_matrix[i][j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  // end for debugging

  // wipe data[0].cov_matrix

  for (i = 1; i <= N; i++)
  {
    for (j = 1; j <= PAWLEY_COR_ELEMENTS; j++)
    {
      data[0].cov_matrix[i][j] = 0.0;
    }
  }


  // now convert the proper cov_matrix into data[0].cov_matrix format

  for (i = 1; i <= N; i++)
  {
    if (cov_matrix[i][i] <= 0.0)
      data[0].cov_matrix[i][1] = 0.03;
    else
      data[0].cov_matrix[i][1] = sqrt(cov_matrix[i][i]);
  }


  // fill in off-diagonal elements

  for (i = 1; i <= N; i++)
  {
    // be aware that j = 2 refers to the correlation between the ith
    // and i+1th reflection

    for (j = 2; j <= PAWLEY_COR_ELEMENTS; j++)
    {
      if (i+j-1 <= N)
      {
        data[0].cov_matrix[i][j] = cov_matrix[i][i+j-1] /
          (data[0].cov_matrix[i][1] * 0.01 *data[0].cov_matrix[i+j-1][1]);
      }

      if (data[0].cov_matrix[i][j] > 99.001) data[0].cov_matrix[i][j] = 99.0;
      if (data[0].cov_matrix[i][j] < -99.001) data[0].cov_matrix[i][j] = -99.0;
    }
  }

}


// For debugging. Saving to file a data[0] structure to .hkl file. 
// Initially meant to be used when a .hkl_weight is used as input
// and after data[0].cov_matrix has been through the convert_hkl_weight_to_hkl
// function

void save_data0_to_hkl_file()
{
	long i,j,k;

	FILE *fp;
  char filename[50];  // used to add OUTPUT_DIRECTORY path

  sprintf(filename, "%sdummy.hkl", OUTPUT_DIRECTORY);

  fp = Fopen(filename, "w");

  for (i = 1; i <= data[0].N; i++)
  {
    // singlet when multiplet[i]=0
    for (j = 0; j <= data[0].multiplet[i]; j++)
    {
      if (j == 0)
      {
        fprintf(fp, "%5d%5d%5d", data[0].h[i], data[0].k[i], data[0].l[i]);
      }
      else
      {
        fprintf(fp, "%5d%5d%5d", data[0].h_multiplet[i][j], data[0].k_multiplet[i][j], data[0].l_multiplet[i][j]);
      }

      fprintf(fp, "%15g  %15g  %7d  ", data[0].ave_intensity[i], data[0].cov_matrix[i][1], i);

      for(k = 2; k <= PAWLEY_COR_ELEMENTS; k++)
      {
        fprintf(fp, "%5g ", data[0].cov_matrix[i][k]);
      }

      fprintf(fp, "\n");
    }
  }
  fclose(fp);
}


// ********************************
// printout some diagnostics
// ********************************

void print_out(long clump_index)
{
	long i,j;

	FILE *fp;
  char filename[50];  // used to add OUTPUT_DIRECTORY path

  sprintf(filename, "%slogfile.temp", OUTPUT_DIRECTORY);

	if(clump_index == 1) 
	{
		fp = Fopen(filename,"w");
	} 
	else
	{
    fp = Fopen(filename,"a");
	}

	fprintf(fp,"N = %i,\t clump_index = %ld \n\n",data[clump_index].N,clump_index);
	
	fprintf(fp,"**** h,k,l values and multiplet ****\n");
	for(i=1;i<=data[clump_index].N;i++)
		fprintf(fp,"%i\t%i\t%i\t%i\n",data[clump_index].h[i],data[clump_index].k[i],
		data[clump_index].l[i], data[clump_index].multiplet[i]);
	fprintf(fp,"\n");
	

	fprintf(fp,"**** covariance values ****  log_det_cov %g\n", 
		data[clump_index].log_det_cov);
	for(i=1;i<=data[clump_index].N;i++)
	{
		for(j=1;j<=data[clump_index].N;j++)
			fprintf(fp,"%g\t",data[clump_index].cov_matrix[i][j]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");


	fprintf(fp,"**** inv_cov values ****\n");
	for(i=1;i<=data[clump_index].N;i++)
	{
		for(j=1;j<=data[clump_index].N;j++)
			fprintf(fp,"%g\t",data[clump_index].inv_cov[i][j]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	

	fprintf(fp,"**** intensity values ****\n");
	for(i=1;i<=data[clump_index].N;i++)
		fprintf(fp,"%f\n",data[clump_index].ave_intensity[i]);
	fprintf(fp,"\n\n");

	fclose(fp);
}
