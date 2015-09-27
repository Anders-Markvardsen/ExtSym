# include <stdlib.h>		// (to include exit(1)) 
# include <stdio.h>
# include <string.h>

# include <math.h>
# include "../lib/nr.h"
# include "../lib/nrutil.h"

#include "read_in_parameters.h"

# include "main_header.h"
# include "display.h"

/***********************************************************************
*
*  main_control -- Contains the main function for ExtSym.
*                  See http://www.markvardsen.net/projects/ExtSym/main.html
*                  for more info on ExtSym
*
*  Author : Anders J. Markvardsen 
*
*  Contributor : John Johnston
*
*  Credit and dependency : Numerical recipies in C
*
*************************************************************************/							

extern double Ave_Int;			// mean prior intensity
							
extern store_data *data;			// where input data are kept (assigned in get_data.c)
extern store_new_data *new_data;	// where symmetry modified data are kept (get_data.c)

extern struct description *symmetries;  // (initialize_sym_struct.c)
extern int num_extinction_groups; // (initialize_sym_struct.c)


// function declarations

void print_diagnosis_of_sym(store_data *data, store_new_data *new_data, long total_clump,
														char name[]);	 

double get_mean_prior_value(store_data *data, long total_blocks,
														long *num_uncorrelated_counter); 
void printout_MC_information(store_data *data, store_new_data *new_data, long total_blocks,
														 char name[]);

FILE *Fopen(char *fname, char *fmode); // from read_in_parameters.c 

void main(void)
{
	long i;                         // simple indices

	long num_uncorrelated_counter;  // count the number of uncorrelated peaks

	double mean_prior_value;		// mean prior value I will end up using

	long base_num_acceptable;		// counts number of acceptable integrals

	char special_sym_option;  // for debugging and investigation of specific sym
	char special_sym_name[20]; // (the code for the specific sym option is commented out 1/2-2001, 
                              // but I may want to use this code latter for debugging purposes) 	

	char infilename[1000];   // input file name

	double sigma_multiply;
	double correlation_cutoff;
	char laue_class_symbol[2];
	char prior_mean_value_option;

	long total_blocks;				// total number of blocks(or clumps)


	// open file   (probably put this somewhere else)	

#ifdef SAVE_PEAK_IGNORED
  
  char filename[50];
	FILE *peak_ignored;  


  sprintf(filename, "%speak_ignored.asc", OUTPUT_DIRECTORY);
  peak_ignored = Fopen(filename,"w");	
#endif

	//printf("%g\n\n", for_large_negative(10, -10));
	//exit(1);

	// read in input file

	read_in_parameters(&sigma_multiply, &correlation_cutoff, 
		infilename, &laue_class_symbol,
		&prior_mean_value_option, &mean_prior_value, &special_sym_option, special_sym_name);


	// get data

	get_data(&total_blocks, infilename, sigma_multiply, 
		correlation_cutoff);
	

	// initialize symmetries

	initialize_sym_struct(laue_class_symbol);


	//
	// mean prior value
	//

  // initialize number of uncorrelated peaks
	
	num_uncorrelated_counter = 0; 


	// decide on how to obtain mean prior value

//	if (prior_mean_value_option == 'y')
//	{
		Ave_Int = get_mean_prior_value(data, total_blocks, &num_uncorrelated_counter);
//	}
//	else
//	{
		// set Ave_Int to value assigned to local variable mean_prior_value assigned in read_in_parameters

//		Ave_Int = mean_prior_value; 
//	}


  //
	// calculate the base integral values
	//
  printf("\n**********************************************************\n");
	printf("ExtSym Version 2008 (update V1)\n\n");


	printf("Please wait while base integrals are calculated....");

	for (i = 1; i<=total_blocks; i++) 
	{	
		base_integral_value(&data[i]);
	}

	printf(" Thank-you\nPlease wait while remaining integrals are calculated\n");


	//
	// for debugging and investigating special symmetry
	// (probably put in some separate function)

	// commented out 1/2-2001 together with some lines in read_in_parameters.c
/*
//	if (special_sym_option == 'y')
//	{
		valid(special_sym_name);

		for (i = 1; i<=total_blocks; i++) 
		{	
			integral_value(data[i], &new_data[i], special_sym_name, i);
		}
	
		// save files related to this specific sym

		print_diagnosis_of_sym(data, new_data, total_blocks, special_sym_name); 
		printout_MC_information(data, new_data, total_blocks, special_sym_name); //exit(1);
//	} */


	//
	// calculate Bayesian probability table
	//

	cal_bayes_table(data, new_data, total_blocks);


	//
	// print out relevant information (probably put this somewhere else)
	//

	printf("\n----------------------------------------------------------\n");
	printf("Diagnosis of calculation : \n\n");
	//printf("---------------------------\n\n");
/*	if (prior_mean_value_option == 'y' && num_uncorrelated_counter != 0)
	{ */
		printf("The mean prior value was calculated to = %f\n\n", Ave_Int);

/*		printf("(This value was calculated from the mean of %ld uncorrelated peaks)\n\n", num_uncorrelated_counter);
	}
	else if (prior_mean_value_option == 'y' && num_uncorrelated_counter == 0)
	{
		printf("Mean (prior) value used was calculated to be = %f\n", Ave_Int);
		printf("(Taken as mean of all the intensities in input data file)\n\n");
	}
	else
	{
		printf("Mean (prior) value used was = %f\n", Ave_Int);
		printf("(This value was taken from the input parameter file)\n\n");
	} */


	printf("The total number of base integrals calculated was = %ld,\n", total_blocks);
	base_num_acceptable = 0;
	for (i = 1; i <= total_blocks; ++i)
	{
		if (data[i].percent_error <= TOL*100)
		{
			base_num_acceptable++;
		}
	}
	printf("where %ld of these have a relative error bigger than %g percent\n", total_blocks - base_num_acceptable, TOL*100);
	printf("and these %ld integrals are ignored.\n", total_blocks - base_num_acceptable);
	printf("----------------------------------------------------------\n\n");
//	printf("The total number of base integrals\n    with relative error less than %g percent is %ld\n", 
	//	TOL*100, base_num_acceptable);

#ifdef SAVE_PEAK_IGNORED
	// (save peak_ignored in some separate function)

	fprintf(peak_ignored, "\n-----------------------------------------------------------------------------\n\n");
	fprintf(peak_ignored, "\n--------- Peaks ignored because relative error of integral too high ---------\n\n\n");
	for (i = 1; i <= total_blocks; ++i)
	{
		if (data[i].percent_error > TOL*100)
		{
			fprintf(peak_ignored, "Percent error = %f -- therefore ignore intensity %ld to %ld\n\n", 
				data[i].percent_error, data[i].reference_point, 
				data[i].reference_point + data[i].N - 1);
		}
	}

	fclose(peak_ignored);
#endif  
	printf("\n** Probability table is saved in table.asc **");
  printf("\n**********************************************************\n");
	//printf("-------------------------------------------\n");

}			


//----------------------------------------------------------------------------
// print out diagnosis of a particular symmetry and compare it calculated information
// in new_data with the base data in data
//------------------------------------------------------------------------------ 

void print_diagnosis_of_sym(store_data *data, store_new_data *new_data, long total_blocks,
							char name[])
{
	long i, j, k;

	
	// open diagnosis_of_sym.temp
	
  char filename[50];
	FILE *diag_fp;  

  sprintf(filename, "%sdiagnosis_of_sym.temp", OUTPUT_DIRECTORY);
  diag_fp = Fopen(filename,"w");	
	
	
	//	printf("**************   ****************\n\n");
	
	fprintf(diag_fp, "Compare BASE with SYM = %s\n\n", name);

	fprintf(diag_fp,"Clump No.  Base reflections  multiplet reflection   ");
	fprintf(diag_fp,"remaining reflections  multiplet reflection\n\n");
	for (i = 1; i <= total_blocks; ++i)
	{
		for (j = 1; j <= data[i].N; j++)
		{
			if (j == 1)
			{
				fprintf(diag_fp,"   %ld            ", i);
			}
			else
			{
				fprintf(diag_fp,"                ");
			}
			fprintf(diag_fp,"%d  %d  %d         ", data[i].h[j], data[i].k[j], data[i].l[j]);
			if (data[i].multiplet[j] >= 1)
			{
				for (k = 1; k <= data[i].multiplet[j]; k++)
				{
					fprintf(diag_fp,"   %d  %d  %d               ", data[i].h_multiplet[j][k], data[i].k_multiplet[j][k], 
						data[i].l_multiplet[j][k]);
				}
			}
			else
			{
				fprintf(diag_fp,"                         ");
			}
			if (j <= new_data[i].N && new_data[i].N > 0)
			{
				fprintf(diag_fp,"%d  %d  %d              ", new_data[i].h[j], new_data[i].k[j], 
					new_data[i].l[j]);
				if (new_data[i].multiplet[j] >= 1)
				{
					for (k = 1; k <= data[i].multiplet[j]; k++)
					{
						fprintf(diag_fp,"%d  %d  %d         ", new_data[i].h_multiplet[j][k], 
							new_data[i].k_multiplet[j][k], new_data[i].l_multiplet[j][k]);
					}
				}
				else
				{
					// new line next
				}
			}
			else 
			{
				fprintf(diag_fp,"                                   ");
			}
			fprintf(diag_fp, "\n");
		}
		if (new_data[i].any_changes_of_block == YES && data[i].percent_error <= TOL*100)
		{
			fprintf(diag_fp,"\n * log(Bas_Int) = %g  log(Sym_Int) = %g Sym_int/Bas_int = %g", 
				data[i].base_integral_value, new_data[i].new_integral_value,
				new_data[i].new_integral_value-data[i].base_integral_value);
		}
		else
		{
			fprintf(diag_fp,"\n * Bas_Int = %g", 
				data[i].base_integral_value);
		}
		fprintf(diag_fp,"\n\n\n");
	}
	fprintf(diag_fp,"\n");
	
	
	fclose(diag_fp);
}



//---------------------------------------------------------------------------------
// return a value for the mean prior value
//---------------------------------------------------------------------------------
double get_mean_prior_value(store_data *data, long total_blocks, long *num_uncorrelated_counter)
{
	long i,j;

	double mean_prior_value_local = 0.0;		// mean prior value	
	
	long dummy_counter;        // used for the case of no uncorrelated peaks

#ifdef SAVE_UNCORRELATED_PEAK
	// open diagnosis_of_sym.temp
	
  char filename[50];
	FILE *uncor_fp;  

  sprintf(filename, "%suncorrelated_peaks.asc", OUTPUT_DIRECTORY);
  uncor_fp= Fopen(filename,"w");


	fprintf(uncor_fp, "List of reflections of the input data file\n");
	fprintf(uncor_fp, "which were catagorized as uncorrelated intensities\n\n");
	fprintf(uncor_fp, "h    k    l     Intensity\n\n");
#endif

	// save the uncorrelated intensities and calculated mean prior value
	// To be present h, k, l must all be even!!

	for (i=1; i<=total_blocks; ++i)
	{
		// if peak uncorrelated then 

		if (data[i].N == 1 && data[i].multiplet[1] == 0 &&
			1 == t0(data[i].h[1], data[i].k[1], data[i].l[1], 
			symmetries[num_extinction_groups].symmetry) )
		{

#ifdef SAVE_UNCORRELATED_PEAK
			// print out intensity to file

			fprintf(uncor_fp, "%d    %d    %d    %10.3g\n", data[i].h[1], data[i].k[1],
				data[i].l[1], data[i].ave_intensity[1]);
#endif

			++(*num_uncorrelated_counter);  // one more uncorrelated peak


			// sum up mean prior value

			mean_prior_value_local += data[i].ave_intensity[1];
		}
	}

	
	// to get mean devide by number of uncorrelated peaks

	if (*num_uncorrelated_counter < 5)
	{
		// no uncorrelated peaks -- display warning message and calculate alternative
		// prior mean value

		mean_prior_value_local = 0.0;
		dummy_counter = 0;

		for (i = 1; i <= total_blocks; i++)
		{
			for (j = 1; j <= data[i].N; j++)
			{
				mean_prior_value_local += data[i].ave_intensity[j];
				dummy_counter++;
			}
		}

	//	printf("lfsa %d\n", dummy_counter);

		mean_prior_value_local /= dummy_counter;

/*		printf("\n\n***********************************************************\n\n");
		printf("Warning -- no uncorrelated peaks to calculate mean prior\n");
		printf("           value. Therefore the mean prior value is calculated\n");
		printf("           as the mean value of all the intensities of the input data\n\n");
		printf("***********************************************************\n\n"); */

#ifdef SAVE_UNCORRELATED_PEAK
		// also write to file

		fprintf(uncor_fp, "\n\n***********************************************************\n\n");
		fprintf(uncor_fp, "Less than 5 uncorrelated peaks to calculate mean prior value.\n");
		fprintf(uncor_fp, "Therefore the mean prior value used was the mean value\n");
		fprintf(uncor_fp, "of all the intensities of the input data which was\n");
		fprintf(uncor_fp, "calculated to %g\n\n", mean_prior_value_local);
		fprintf(uncor_fp, "***********************************************************\n\n");
#endif

		return mean_prior_value_local;
	}
/*    else if(*num_uncorrelated_counter < 3)
	{
		// display warning

		printf("\n\n***********************************************************\n\n");
		printf("Warning -- less than 3 uncorrelated peaks to be used for\n");
		printf("           determining a value for the mean prior value\n\n");
		printf("***********************************************************\n\n");

		// also write to file

		fprintf(uncor_fp, "\n\n***********************************************************\n\n");
		fprintf(uncor_fp, "Warning -- less than 3 uncorrelated peaks to be used for\n");
		fprintf(uncor_fp, "           determining a value for the mean prior value\n\n");
		fprintf(uncor_fp, "***********************************************************\n\n");
	} */
	else 
	{
		// else ok
	}

	mean_prior_value_local /= ((double) *num_uncorrelated_counter);
 
#ifdef SAVE_UNCORRELATED_PEAK
	// write mean value to uncor_fp

	fprintf(uncor_fp, "\n\nThe mean value of these uncorrelated intensities is = %g\n", mean_prior_value_local);


  // finished writing file

	fclose(uncor_fp);
#endif

	// return mean prior value

	return mean_prior_value_local;
}




void printout_MC_information(store_data *data, store_new_data *new_data, long total_blocks,
							 char name[])
{
	long i, j;


	// open MC_information.temp

  char filename[50];
	FILE *mc_fp;  

  sprintf(filename, "%sMC_information.temp", OUTPUT_DIRECTORY);
  mc_fp= Fopen(filename,"w");			


	fprintf(mc_fp, "BASE INTEGRALS\n\n");


	for (i = 1; i <= total_blocks; i++)
	{
		if (data[i].N > 1) 
		{
			fprintf(mc_fp,"Clump No = %ld        percent error = %f     N = %d\n", i, data[i].percent_error,
				data[i].N);
			fprintf(mc_fp, "name of highway : %s\n", data[i].highway_name);
			fprintf(mc_fp,"extra_mean_values[i]        monte_ave[i]            Cov[i][i]\n");
			for (j = 1; j <= data[i].N; ++j)
			{
				fprintf(mc_fp, "%g                           %g                   %g\n", data[i].extra_mean_values[j], data[i].monte_ave[j],
					data[i].cov_matrix[j][j]);
			}
			fprintf(mc_fp, "\n\n");
		}
	}


    fprintf(mc_fp, "SYM INTEGRALS with sym = %s\n\n", name);
	for (i = 1; i <= total_blocks; i++)
	{
		if (new_data[i].N > 1 && new_data[i].any_changes_of_block) 
		{
			fprintf(mc_fp,"Clump No = %ld        percent error = %f     N = %d\n", i, new_data[i].percent_error,
				new_data[i].N);
			fprintf(mc_fp,"extra_mean_values[i]        monte_ave[i]            Cov[i][i]\n");
			for (j = 1; j <= new_data[i].N; ++j)
			{
				fprintf(mc_fp, "%g                           %g                   %g\n", new_data[i].extra_mean_values[j], new_data[i].monte_ave[j],
					new_data[i].cov_matrix[j][j]);
			}
			fprintf(mc_fp, "\n\n");
		}
	}

	fclose(mc_fp);
}



