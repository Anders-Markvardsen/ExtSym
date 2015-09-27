void get_data (long *total_blocks, char infilename[],
							 double sigma_multiply, double correlation_cutoff);
void initialize_sym_struct(char laue_class_symbol[]);
void read_in_parameters(double *sigma_multiply, 
												double *correlation_multiply,
												char *infilename,
												char *laue_class_symbol,
												char *prior_mean_value_option,
												double *prior_mean_value,
												char *special_sym_option,
												char special_sym_name[]);


struct description
{
    char* symmetry;
    char* name;
};


typedef struct              // definition of the store_data structure
  {
  int    *h;                   // pointer to an array of h's
  int    *k;                     
  int    *l;

  int *multiplet;									//		 triplet = 2
																	//     doublet = 1
																	//     singlet = 0 

  int    **h_multiplet;           // the additional hkl values when multiplet
  int    **k_multiplet;
  int    **l_multiplet;

  double *ave_intensity;       // pointer to an array of the likelihood 
                               // expectation values of the intensities


  double **inv_cov;            // pointer to inverse of the likelihood
                               // covariance matrix of the intensities

  double **cov_matrix;

  int N;   // the number of intensities

  double base_integral_value;  // integral value for the base extinction group

  double log_det_cov;           // log of determinant of covariance matrix
								// used to calculate base_integral_value

  double percent_error;			// percentage error in monte carlo integral

  double *monte_ave;			// used in error-checking monte carlo algorithm

  double *extra_mean_values;	// used in error-checking monte carlo algorithm
								// Ones I now that monte carlo algorithm is optimized
								// I can delete monte_ave and extra_mean_values from 
								// this structure.
  char highway_name[80];

	long reference_point;		// which number in the 6th column of the .hkl file is the starting
													// intensity of the block. This option is only used for printing
													// purposes
} store_data;


typedef enum {NO, YES} yesno;

typedef struct 
  {
  int    *h;                   // pointer to an array of h's
  int    *k;                   // used purely for displaying purposes  
  int    *l;

  int    **h_multiplet;           // the hkl values of multiplets highter than singlets
  int    **k_multiplet;									// used purely for displaying purposes 
  int    **l_multiplet; 
  

  double *ave_intensity;       // pointer to an array of the likelihood 
                               // expectation values of the intensities

  double **cov_matrix;         // pointer to a matrix of the likelihood
                               // covariance matrix of the intensities

  double **inv_cov;            // pointer to inverse of the likelihood
                               // covariance matrix of the intensities

	// this one below can be made redundant by using multiplet where allowed 
	// is used

  int *allowed;                // pointer to an array of allowed flags
                               //
                               //    not allowed = 0 
                               //    singlet     = 1
                               //    doublet     = 2 

  int *multiplet;									//		 triplet = 2
																	//     doublet = 1
																	//     singlet = 0 

  int N;                       // the number of intensities

  double new_integral_value;   // integral value of the new extinction group

  double log_det_cov;			// log of deterninant of covariance matrix
								// used to calculate new_integral_value

  yesno any_changes_of_block; // if the block has exactly the same number
											 // of singlets and doublet then 
											 //		any_changes_of_block = NO
											 // otherwise
											 //		any_changes_of_block = YES

//  long num_delta_priors;		// number of delta priors

  double percent_error;			// percentage error in monte carlo integral

  double *monte_ave;			// used in error-checking monte carlo algorithm

  double *extra_mean_values;	// used in error-checking monte carlo algorithm
								// Ones I now that monte carlo algorithm is optimized
								// I can delete monte_ave and extra_mean_values from 
								// this structure.

  char highway_name[80]; // which method of MC-integration turns to be best
} store_new_data;


// changed 21/2-2001 from # define to extern double

extern double TOL;          // the fractional accuracy that the estimate is
                           // calculated




