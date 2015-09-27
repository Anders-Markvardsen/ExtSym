# include <stdlib.h>		// (to include exit(1)) 
# include <stdio.h>
# include <string.h>
# include <math.h>

#include "read_in_parameters.h"

# include "main_header.h"
# include "display.h"



typedef struct              
{
	char name[20];
	double value;
	int num_integrals_changed;
} table; 


static table *bayes_table;
static table *sorted_table;

extern struct description *symmetries;  // (initialize_sym_struct.c)
extern int num_extinction_groups; // (initialize_sym_struct.c)

void save_bayes_table(void);

FILE *Fopen(char *fname, char *fmode); // from read_in_parameters.c 

//----------------------------------------------------------------------------------
// function takes the external parameters: num_extinction_groups, symmetries[].
//----------------------------------------------------------------------------------

void cal_bayes_table(store_data *data, store_new_data *new_data, long total_blocks)
{
	long i;
	long sym_index; 

	long num_integrals_changed;

	double total_base_integral_value;
	double total_sym_integral_value;
	double ratio_sym_base;


	bayes_table = (table *) malloc((num_extinction_groups+1)*sizeof(table));
	sorted_table = (table *) malloc((num_extinction_groups+1)*sizeof(table));


	// main loop

	for (sym_index = 1; sym_index <= num_extinction_groups; ++sym_index)
	{
		// calculate new_data for symmetry with index number sym_index

		for (i = 1; i <= total_blocks; ++i)
		{
			integral_value(data[i], &new_data[i], symmetries[sym_index].symmetry, i);
		}
	

		// initialize integral_values

		total_base_integral_value = 0.0;
		total_sym_integral_value = 0.0;


		// calculating the total integral values

		num_integrals_changed = 0;
		for (i = 1; i <= total_blocks; ++i)
		{
			if (new_data[i].any_changes_of_block == YES && data[i].percent_error <= TOL*100)
			{
				total_base_integral_value += data[i].base_integral_value;
				total_sym_integral_value += new_data[i].new_integral_value;
				num_integrals_changed++;
			}
		}


		// and the (log of the) sym/base ratio for this sym is

		ratio_sym_base = total_sym_integral_value - total_base_integral_value;

		strcpy(bayes_table[sym_index].name, symmetries[sym_index].name);
		bayes_table[sym_index].value = ratio_sym_base;
		bayes_table[sym_index].num_integrals_changed = num_integrals_changed;
		strcpy(sorted_table[sym_index].name, symmetries[sym_index].name);
		sorted_table[sym_index].value = ratio_sym_base;
		sorted_table[sym_index].num_integrals_changed = num_integrals_changed;

	}

	save_bayes_table();
}


void save_bayes_table(void)
{
	int i;

	int a, b; // for bubble sort
	double double_temp; // for bubble sort
	char string_temp[20]; // for bubble sort
	int int_temp; // for bubble sort

	double prob_const;			// to obtain procentage 

//	FILE *syms_fp;	
	FILE *table_fp;				// to generate table ascii file
	
	char filename[50];

	// open diagnosis_of_sym.temp
	
//	if((syms_fp = fopen("ratio_information.asc","w")) == NULL)
//	{
//		printf("Error open ratio_information.asc\n");
//		printf("Now exiting program\n");
//		exit(1);
//	}

	// open table.asc
  
  sprintf(filename, "%stable.asc", OUTPUT_DIRECTORY);
	table_fp = Fopen(filename, "w");


	// sort table using bubble sort

	for (a = 2; a <= num_extinction_groups; ++a)
	{
		for (b = num_extinction_groups; b >= a; --b)
		{
			if (sorted_table[b-1].value < sorted_table[b].value)
			{
				// exchange elements
				
				double_temp = sorted_table[b-1].value;
				sorted_table[b-1].value = sorted_table[b].value;
				sorted_table[b].value = double_temp;

				strcpy(string_temp, sorted_table[b-1].name);
				strcpy(sorted_table[b-1].name, sorted_table[b].name);
				strcpy(sorted_table[b].name, string_temp);

				int_temp = sorted_table[b-1].num_integrals_changed;
				sorted_table[b-1].num_integrals_changed = sorted_table[b].num_integrals_changed;
				sorted_table[b].num_integrals_changed = int_temp;
			}
		}
	}	

	// find probability constant. (Under the assumpsion that 
	// the Laue class form a complete sample space then it
	// posssible to list the relative probabilies as absolute
	// probabilities.) 

	prob_const = 0.0;
	for (i = 1; i <= num_extinction_groups; ++i)
	{
		prob_const += exp(sorted_table[i].value - sorted_table[1].value);
	}

  // print out ralevant information


	fprintf(table_fp, " When reading the table below keep in mind that more than one space group may\n");
	fprintf(table_fp, " have the same Extinction Symbol. For more info about this see the International\n"); 
	fprintf(table_fp, " Tables for Crystallography Vol A. For further info about how to interpret this file\n");
	fprintf(table_fp, " visit: www.markvardsen.net/projects/ExtSym/explain-output.html\n\n");
  fprintf(table_fp, " The extinction symbols are ranked in order of probability with the most probable\n");
  fprintf(table_fp, " symbol being listed first\n\n");

    fprintf(table_fp, " Extinction Symbol      log-probability score\n\n");

	for (i = 1; i <= num_extinction_groups; ++i)
	{
		// print out information to table.asc

		fprintf(table_fp, "%12s               %10.6g\n", sorted_table[i].name, sorted_table[i].value);
	}


//	fclose(syms_fp);
	fclose(table_fp);
}