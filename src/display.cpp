//# include <iostream>
//# include <string>
//# include <vector>

# include "main_header.h"
# include "display.h"



struct table             
{
//	string name;
	double value;
	int num_integrals_changes;
}; 

extern struct description *symmetries;  // (initialize_sym_struct.c)
extern int num_extinction_groups; // (initialize_sym_struct.c)



//----------------------------------------------------------------------------------
// function takes the external parameters: num_extinction_group, symmetries[].
//----------------------------------------------------------------------------------

void compare_many_syms(store_data *data, store_new_data *new_data, long total_blocks)
{
	long i;
	long sym_index; 

	long num_integrals_changed;

	double total_base_integral_value;
	double total_sym_integral_value;
	double ratio_sym_base;

//	FILE *syms_fp;	
//	FILE *table_fp;				// to generate table ascii file
	
	
	// if I want to change default value

	//num_extinction_groups = 50;


	// open diagnosis_of_sym.temp
	
//	if((syms_fp = fopen("ratio_information.asc","w")) == NULL)
//	{
//		printf("Error open ratio_information.asc\n");
//		printf("Now exiting program\n");
//		exit(1);
//	}

	// open table.asc

//	if((table_fp = fopen("table.asc","w")) == NULL)
//	{
//		printf("Error open table.asc\n");
//		printf("Now exiting program\n");
//		exit(1);
//	}


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


		// print out ralevant information

//		fprintf(syms_fp, "Symmetry %s(%s):  ", symmetries[sym_index].name, 
//			symmetries[sym_index].symmetry);
//		fprintf(syms_fp, "Number integrals changed %ld\n", num_integrals_changed);
//		fprintf(syms_fp, "  base = %g  sym = %g  sym/base = %g\n\n", 
//			total_base_integral_value, total_sym_integral_value, ratio_sym_base);

		// print out information to table.asc

//		fprintf(table_fp, "%s , %g\n", symmetries[sym_index].name, 
//			ratio_sym_base);
	}

//	fclose(syms_fp);
//	fclose(table_fp);
}
