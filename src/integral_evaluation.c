# include <stdlib.h>		// (to include exit(1)) 
# include <stdio.h>
# include <math.h>
# include "../lib/nrutil.h"
# include "../lib/nr.h"

# include "data_types.h"   
# include "integral.h"      


extern double Ave_Int; //(three.c)

double sample  (int N, double **Cov, double *Ave, int *multiplet, 
				double **inv_Cov, double *percent_error, double *monte_ave,
				double *extra_mean_values, char highway_name[]); //(sample.c)

double direct_integration(double Variance, double Average, int multiplet);		//(sample.c)

void symm_hkl(char symmetry[], store_data data, store_new_data *new_data_ptr); //(symmetry.c)

void no_contract(store_data data, store_new_data *new_data);

void contract(store_data data, store_new_data *new_data, long clump_number);

void reduce_matrix_and_vector(int from_N, int to_N, double **matrix, double *vector,
							  int *project_vec);

double eval_integral(int N, double **cov_matrix, double *ave_intensity, int *multiplet, 
			  double **inv_cov, double log_det_cov, double *percent_error, 
			  double *monte_ave, double *extra_mean_values, char highway_name[], char store_option,
				int *allowed, long clump_number); 

//------------------------------------------------------------------------------------------
// Do the base work
//------------------------------------------------------------------------------------------

void base_integral_value(store_data *data)
{
	int dummy; // used as dummy argument for allowed

	long clump_dummy = 0; // used as dummy argument for clump_number

	data->base_integral_value = eval_integral(data->N, data->cov_matrix, 
		data->ave_intensity,
		data->multiplet, 
		data->inv_cov, 
		data->log_det_cov, 
		&(data->percent_error),
		data->monte_ave,
		data->extra_mean_values,
		data->highway_name,
		'n',
		&dummy,
		clump_dummy);
}


//----------------------------------------------------------------------------------
// Distributing work
// The argument clump_number is used when store integral value
//----------------------------------------------------------------------------------

void integral_value(store_data data, store_new_data *new_data, char *symmetry, long clump_number)
{
	long i,j;	// dummy indices


	// find new_data.allowed[], new_data.any_changes_of_block, new_data.N, 
	// new_data.h,k,l, new_data.h_multiplet,k_multiplet,l_multiplet and also
	// new_data.multiplet[] is assigned in symm_khl().
	
	symm_hkl(symmetry, data, new_data);


	// 

	if (new_data->any_changes_of_block == YES  && data.percent_error <= TOL*100.0)
	{
		if (new_data->N != data.N)
		{
			// contract
			contract(data, new_data, clump_number);
		}
		else if (new_data->N == data.N)  
		{
			// when contraction not needed just do

			new_data->new_integral_value = eval_integral(new_data->N, 
				data.cov_matrix, 
				data.ave_intensity,
				new_data->multiplet, 
				data.inv_cov, 
				data.log_det_cov,
				&(new_data->percent_error),
				new_data->monte_ave,
				new_data->extra_mean_values,
				new_data->highway_name,
				'y',
				new_data->allowed,
				clump_number); 

			// for displaying purposes

			for (i = 1; i <= data.N; ++i)
			{
				new_data->ave_intensity[i] = data.ave_intensity[i];
			}
			for (i = 1; i <= data.N; ++i)
			{
				for (j = 1; j <= data.N; ++j)
				{
					new_data->cov_matrix[i][j] = data.cov_matrix[i][j];
				}
			}
		}
		else
		{
			// error 
			printf("Error in integral_value\n");
			exit(1);
		}
	}
	else
	{
		// no change of block hence don't do anything
	}
}



//--------------------------------------------------------------------------------------
// Constract
// (generation of new multiplet vector is dealt with within symm_hkl())
// The argument clump_number is adding in include the store integral option
//--------------------------------------------------------------------------------------

void contract(store_data data, store_new_data *new_data, long clump_number)
{
	int i,j;                    // simple matrix indices and counters
	
	double *reduced_temp_vector;	// a vector required in the calculation
									// of the averages in the new space
		
	double x;                     // general storage variable
		
	long allowed_index;				// number of allowed peaks control index

	double gaus_left_over_term = 0.0;	// initialized here for convenience

	double **temp_inv_cov;

//	const double my_PI = acos(-1);

	/* comment out double check 2/11-99
	double dummy_ave_intensity; // used to double check that the new ave_intensity has 
								// calculated correctly	
	*/
	
	// new_data.inv_cov is first initialised to data.inv_cov, then cut down to its
	// required size and finally then also used to calculate new_data.cov_matrix
	
	for(i=1;i<=data.N;i++)
		for(j=1;j<=data.N;j++)
			new_data->inv_cov[i][j] = data.inv_cov[i][j];
		
	// this is used when new_data.ave_intensity is cut down

	for(i=1;i<=data.N;i++)
		new_data->ave_intensity[i] = data.ave_intensity[i];
		
			
	// The vector reduced_vector is created and calculated here for later use

	if (new_data->N) {
		reduced_temp_vector = dvector(1,new_data->N);
	
		allowed_index = 0;
		for (i = 1; i <= data.N; ++i)
		{
			if (new_data->allowed[i] != 0)		// Means column not cut out of covariance
			{									// matrix. 
				++allowed_index;				// Multiply this column by the vector ((I-P)x_0) 
												// and thereby is the allowed_index component
												// of reduced_temp_vector calculated
				reduced_temp_vector[allowed_index] = 0.0;
				for (j = 1; j <= data.N; ++j)
				{
					if (new_data->allowed[j] == 0) // only when allowed = 0 is ((I-P)x_0) 
					{							   // non-zero
						// the allowed_index component of reduced_temp_vector is calculated
						reduced_temp_vector[allowed_index] += data.inv_cov[j][i] *
							data.ave_intensity[j];
					}
				}
			}
		}		
	}
		

	// the input inverse covariance matrix is cut down by removing all
	// the elements corresponding to the unwanted intensities and rearranging
	// the remaining elements in the correct manner. (This is now the inverse
	// of the required matrix).		
	// The same type of thing is done for the averages intensities
	
	reduce_matrix_and_vector(data.N, new_data->N, new_data->inv_cov, 
			new_data->ave_intensity, new_data->allowed);


	// 

	if (new_data->N > 0)		// at least one allowed reflection in clump
	{
		// prepare for the use of inverse_jac which destroy it input matrix, therefore
		// copy the content of new_data->inv_cov to temp_inv_cov and use that matrix 
		// as input argument for inverse_jac

		temp_inv_cov = dmatrix(1, new_data->N, 1, new_data->N);

		for (i = 1; i <= new_data->N; ++i)
			for (j = 1; j <= new_data->N; ++j)
				temp_inv_cov[i][j] = new_data->inv_cov[i][j];
		
		// the cut down NxN matrix new_data.inv_cov is inverted to give the 
		// required covariance matrix. 
		// Notice temp_inv_cov is destroyed by inverse_jac
	
	
		inverse_jac(temp_inv_cov, 
			new_data->cov_matrix, 
			&new_data->log_det_cov, 
			new_data->N);
		// inverse_jac returned to log of the determinant of the input matrix, which is
		// inv_cov, therefore to get the log of the determinant of cov then do
		new_data->log_det_cov *= -1;
	
		// de-allocate temp_inv_cov

		free_dmatrix(temp_inv_cov, 1, new_data->N, 1, new_data->N);
	


	
		// calculation new averages equal to x_0^T + ((I-P)x_0)^T*(C^{-1})_spec*C_reduced	
		// C_reduced is of dimensions new_data.N*new_data.N. (C^{-1})_spec is is a cut down
		// version of C^{-1} where the column for which allowed = 0 have been taken out, 
		// therefore (C^{-1})_spec is of dimension data.N*new_data.N. However (C^{-1})_spec
		// is multiplied on the left by a vector of dimensions 1*data.N, hence 
		// ((I-P)x_0)^T*(C^{-1})_spec is of dimensions 1*new_data.N.
	
		for(i=1; i<=new_data->N; i++)
		{
			for(j=1; j<=new_data->N;j++)
			{
				new_data->ave_intensity[i] 
					+= reduced_temp_vector[j]*new_data->cov_matrix[i][j];
			}
		}
	
/* I don't need this double check any more - 2/11-99 
	
		// make another double check here. This is to check that new ave_intensities
		// have been calculated correctly. For the double check we use that the vector
		// x_0^T + ((I-P)x_0)^T*(C^{-1})_spec*C_reduced is equal to the vector 
		// x_0^T*(C^{-1})_spec*C_reduced. 
		// Hence calculate here first the 1*new_data.N vector, x_0^T*(C^{-1})_spec, by 
		// overwriting reduced_temp_vector and then multiply it with C_reduced to 
		// calculate the new ave_intensities a second time and compare these with the previous 
		// calcuated new ave_intensities. 

		allowed_index = 0;
		for (i = 1; i <= data.N; ++i)
		{
			if (new_data->allowed[i] != 0)  // if true i'th column allowed in correlation matrix
			{
				++allowed_index;
				reduced_temp_vector[allowed_index] = 0.0;
				for (j = 1; j <= data.N; j++)
				{
					reduced_temp_vector[allowed_index] += data.inv_cov[j][i] *
							data.ave_intensity[j]; 
					// reduced_temp_vector now contains x_0^T*(C^{-1})_spec
				}
			}
		}

		for(i=1; i<=new_data->N; i++)
		{
			dummy_ave_intensity = 0.0;
			for(j=1; j<=new_data->N;j++)
			{
				dummy_ave_intensity 
					+= reduced_temp_vector[j]*new_data->cov_matrix[i][j];
			}
			if (fabs( dummy_ave_intensity/new_data->ave_intensity[i]-1.0 ) > 1e-8)
			{
				printf("Error in calculating new ave_intensities in contract()\n");
				printf("fabs( dummy_ave_intensity/new_data->ave_intensity[i]-1.0 ) = %g\n",
					fabs( dummy_ave_intensity/new_data->ave_intensity[i]-1.0 ));
				exit(1);
			}
		}

     endoff double check */


		// free reduced_temp_vector

		free_dvector(reduced_temp_vector, 1, new_data->N);


		// calculate zero order Gaussian left over term

		// gaus_left_over_term is initialize to zero at the beginning of function

		for (i = 1; i <= new_data->N; i++)  // subtract x_new^T*(C^{-1})_reduced*x_new
		{
			x = 0.0;
			for (j = 1; j <= new_data->N; j++)
			{
				x += new_data->inv_cov[i][j] * new_data->ave_intensity[j];
			}
			gaus_left_over_term -= new_data->ave_intensity[i] * x;
		}

	}
	else						// no allowed reflection in clump
	{
		new_data->log_det_cov = 0.0;
	}


	// calculating zero order Gaussian left over term -- continued if new_data.N>0
	
	for (i = 1; i <= data.N; i++)  // add first x_0^T*C^{-1}*x_0
	{
		x = 0.0;
		for (j = 1; j <= data.N; j++)
		{
			x += data.inv_cov[i][j] * data.ave_intensity[j];
		}
		gaus_left_over_term += data.ave_intensity[i] * x;
	}

	
	// finally multiply by -0.5 to get the log of the zero order Gaussian left over term

	gaus_left_over_term *= -0.5;


	// start to add to the integral evaluation
	//

	// add the gaus_left_over_term to the integral_value

	new_data->new_integral_value = gaus_left_over_term;

	// add remaining stof

	new_data->new_integral_value += eval_integral(new_data->N, new_data->cov_matrix, 
				new_data->ave_intensity,
				new_data->multiplet, 
				new_data->inv_cov, new_data->log_det_cov,
				&(new_data->percent_error),
				new_data->monte_ave,
				new_data->extra_mean_values,
				new_data->highway_name,
				'y',
				new_data->allowed,
				clump_number); 
}


//-------------------------------------------------------------------------------------
// Project matrix and vector down from from_N to to_N in accordence with the projection
// vector project_vec, which has elements to one for dimensions which are kept and zero
// otherwise.
//
// NOTICE: this functions overwrites matrix and vector.
//-------------------------------------------------------------------------------------

void reduce_matrix_and_vector(int from_N, int to_N, double **matrix, double *vector,
							  int *project_vec)
{
	int allowed_countdown;  // counts down number of allowed dimensions

	int i,j,n;

	
	// inializing allowed_countdown

	allowed_countdown = from_N;
	
	for(n=from_N; n>=1; n--)
	{
		if(project_vec[n]==0) 
		{
			if (n == allowed_countdown)		// then simply do not consider that 
			{								// reflection. This condition will be valid
											// until the first allowed intensity starting
											// from the bottom, in which case from there
											// unwards n will be smaller then 
											// allowed_countdown
				allowed_countdown--; 	
			}
			else						    // the intensities and covariance elementes
			{								// below are all copied one step up to fill
											// the gap left be the missing reflection

				// cut down inverse covariance matrix
				
				for(i=n; i<=allowed_countdown-1; i++)
				{
					for(j=1; j<=allowed_countdown; j++)
					{
						matrix[i][j]	= matrix[i+1][j];
					}
				}
				
				for(i=n; i<=allowed_countdown-1; i++)
				{
					for(j=1; j<=allowed_countdown-1; j++)
					{
						matrix[j][i] = matrix[j][i+1];     
					}
				}	
				
				// cut down the average intensities
				// (used when calculating final new ave_intensities later)
				
				for(i=n;i<=allowed_countdown-1;i++)
				{
					vector[i] = vector[i+1];
				}
								
				// every time the elements corresponding to an 
				// unwanted intensity are removed allowed_countdown decreases by one
				
				allowed_countdown--; 
			}
		}
	}
		
	
	// run double check

	if (allowed_countdown != to_N)
	{
		printf("Error in reduce_matrix_to_vector\n");
		exit(1);
	}

}


//-----------------------------------------------------------------------------------------
// store_option, allowed and clump_number are added to include the store integral option
//-----------------------------------------------------------------------------------------

double eval_integral(int N, double **cov_matrix, double *ave_intensity, int *multiplet, 
			  double **inv_cov, double log_det_cov, double *percent_error, 
			  double *monte_ave, double *extra_mean_values, char highway_name[], char store_option,
				int *allowed, long clump_number) 
{
	double return_value;			// value to be returned

	double int_value;					// used to temporary store integral value

	const double my_PI = acos(-1);

	// start to add to the integral evaluation
	//
	// The first term comes from writing the not normalised probability as
	// (Please notice the original Gaussian likelihood normalization constant is never
	//  included since it is constant throughout the evaluation of all the probabilities.
	//  The sqrt(2*pi)^N * sqrt(|det(Cov)|) below refers to the normalization constant
	//  of the Gaussian sampler in <prior>_gaus).
	//
	//	       [sqrt(2*pi)^N * sqrt(|det(Cov)|) * mu^{-N}]*[mu^N * <prior>_gaus]  ,
	//
	// where then data.N*0.5*log(2.0*M_PI) below is the log of sqrt(2*pi)^N.
	//
	// The second term: this may not be an intuiative way of doing it, but mu^{-N} has
	// been brought outside of <prior>_gaus! 
	// In case of a multiplet we still have to multiply by an extra 1/mu for each multiplet 
	// (or add -log(mu)) - this is done in doublet_integration and within the 
	// function generate() in sample.c. 
	// mu is the prior mean value (or Ave_Int here).

	return_value = ((double) N) * ( 0.5*log(2*my_PI) - log(Ave_Int) );
	
	
	// next add the square root of log of the determinant of the covariance matrix which is 
	// calculated in get_data.c
	
	return_value += 0.5 * log_det_cov;

	
	// finally add log( mu^{N}*<prior>_guas ) 

	switch(N)
	{
        // goto next clump
		
	case 0: 
		{
			// then add nothing
			break;
		}
		

        // use erfc function to calculate mu*<prior>_gaus
		
	case 1:
		{
			// we have an isolated multiplet so we call direct_integration()
			
			return_value += direct_integration(cov_matrix[1][1], ave_intensity[1], multiplet[1]);
			
			break;
		}
		

        // use Monte Carlo integration
		
	default:
		{
			// for base integral do not store integral value

			if (store_option == 'n')
			{
				return_value += sample(N, cov_matrix, ave_intensity, multiplet, inv_cov,
				percent_error, monte_ave, extra_mean_values, highway_name); 
				break;
			}
			else
			{
				// check if integral value is already calculated 

				int_value = done_integral(allowed, clump_number);

				// if it is not then int_value == 0

				if(int_value==0)
	      {
					int_value = sample(N, cov_matrix, ave_intensity, multiplet, inv_cov,
					percent_error, monte_ave, extra_mean_values, highway_name); 
		
					// store integral value
					
					store_integral(allowed, int_value, clump_number);
				}
				
				return_value += int_value;
				break;
			}
		}
	}

	return return_value;
}
