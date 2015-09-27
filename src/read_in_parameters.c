#include <stdio.h>
#include <stdlib.h>		// (to include exit(1)) 
#include <stdarg.h> 
#include "read_in_parameters.h"
#include <stdio.h>

/****************************************************************
* read_in_parameters -- reads in input parameters  
*
* Input --		parameter_input.asc
*
* Returns --	a set of global constants in read_in_parameters.h
* 
***************************************************************/

// some haking here (21/2-2001): also find TOL in data_types.h as extern double TOL
// TOL is read from advanced.asc below

double TOL;



//  crap definition and declarations

#define SAFEFSCANF(FP,FMT,VP,TEXT) \
  if (fscanf(FP,FMT,VP) != 1) Complain(TEXT "\n"); \
  else kill_line(FP);
void Complain(char *fmt, ...);
FILE *Fopen(char *fname, char *fmode);
void kill_line(FILE *fp);


//------------------------------------------------------------
// reads in parameters 
//------------------------------------------------------------
void read_in_parameters(double *sigma_multiply, 
												double *correlation_cutoff,
												char *infilename,
												char *laue_class_symbol,
												char *prior_mean_value_option,
												double *prior_mean_value,
												char *special_sym_option,
												char special_sym_name[]) 
{
  char filename[50];

  FILE *mepar; 
	FILE *fp_av;

  sprintf(filename, "%sparameter_input.asc", INPUT_DIRECTORY);
  mepar = Fopen(filename, "r");

  sprintf(filename, "%sadvanced.asc", INPUT_DIRECTORY);
  fp_av = Fopen(filename, "r"); 


  // scan in parameters
  
	kill_line(mepar);
	kill_line(mepar);
	kill_line(mepar);
	kill_line(mepar);
	kill_line(mepar);


	// reading in the input file name taking into account that this
	// name may contain spaces like "Program Files"

  {
    int i, ch;

    char hkl_ext[5];
    char hkl_ext_compare[] = ".hkl";
    char hkl_weight_ext[12];
    char hkl_weight_ext_compare[] = ".hkl_weight";
    int ii, j;

		
	  // find i value which currently 'ends' infilename

	  ch = fgetc( mepar );
	  for( i=0; (i < 996 ) && ( feof( mepar ) == 0 ); i++ )
	  {
		  infilename[i] = (char)ch;

		  if (infilename[i] == '\0' || infilename[i] == '\n') break;

      ch = fgetc( mepar );
	  }

	  // step down infilename until first letter and insert '\0'

	  for (;i >= 1;i--)
	  {
		  if (infilename[i] != ' ' && infilename[i] != '\t' && infilename[i] != '\n'
			  && infilename[i] != '\0')
		  {
			  infilename[i+1] = '\0';
			  break;
		  }
	  }
	
    // at this point i point to last letter of filename

    // check if filename ends with .hkl

    input_file_format = 'n';  // 'n' means here not a valid format
    ii = 0;
    if (i >= 4)
    {
      for (j = i-3; j <= i+1; j++)
      {
        hkl_ext[ii] = infilename[j];
        ii++;
      }

      if (strcmp(hkl_ext, hkl_ext_compare) == 0) input_file_format = 'c';
    }

    if (input_file_format == 'n')
    {
      ii = 0;
      if (i >= 11)
      {
        for (j = i-10; j <= i+1; j++)
        {
          hkl_weight_ext[ii] = infilename[j];
          ii++;
        }

        if (strcmp(hkl_weight_ext, hkl_weight_ext_compare) == 0) 
          input_file_format = 'i';        
      }
    }

    
    if (input_file_format == 'n')
    {
			printf("\n\n***********************************************************\n\n");
			printf("ERROR IN PARAMETER_INPUT.ASC FILE:\n");	
			printf("  The input powder data file must have either the extension\n");
			printf("  .hkl or .hkl_weight.\n");
			printf("***********************************************************\n\n");

      exit(1);
    }

  } // end of scope

	kill_line(mepar);
	kill_line(mepar);

  SAFEFSCANF(mepar, "%s", laue_class_symbol, 
		"Couldn't read laue_class.");

	fclose(mepar);	


	kill_line(fp_av);
	kill_line(fp_av);
	kill_line(fp_av);
	kill_line(fp_av);
	kill_line(fp_av);
	kill_line(fp_av);
	kill_line(fp_av);

  SAFEFSCANF(fp_av, "%lg", sigma_multiply, 
		"Couldn't read sigma_multiply.");

	kill_line(fp_av);
	kill_line(fp_av);
	kill_line(fp_av);
	kill_line(fp_av);

  SAFEFSCANF(fp_av, "%lg", correlation_cutoff, 
		"Couldn't read correlation_cutoff.");
  SAFEFSCANF(fp_av, "%ld", &max_monte_iteration, 
		"Couldn't read max_monte_iteration.");
	SAFEFSCANF(fp_av, "%lg", &TOL, 
		"Couldn't read TOL.");
	TOL /= 100.0;  // read in as a percentage

	// commented out 21/2-2001 
/*  SAFEFSCANF(mepar, "%c", prior_mean_value_option, 
		"Couldn't read prior_mean_value_option.");
  SAFEFSCANF(mepar, "%lg", prior_mean_value, 
		"Couldn't read prior_mean_value."); */

	// commented out 1/2-2001 together with some lines in mail_control.c
/*
  SAFEFSCANF(mepar, "%c", special_sym_option, 
		"Couldn't read special_sym_option."); 
  SAFEFSCANF(fp_av, "%s", special_sym_name, 
		"Couldn't read special_sym_name."); */

	// commented out 1/8-2001
/*	SAFEFSCANF(fp_av, "%c", &input_file_format, 
		"Couldn't read input_file_format."); */
	// instead the two lines below are substituted with



  // close streams

  fclose(fp_av);
}


//----------------------------------------------------------
//  crap functions below
// ---------------------------------------------------------

void Complain(char *fmt, ...) { 
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);

  exit(1);
}


FILE *Fopen(char *fname, char *fmode) {
  FILE *fp;
  
  if ((fp=fopen(fname, fmode)) == NULL)
  {
    Complain("Couldn't open file \"%s\" in mode \"%s\".\n",
	     fname, fmode);
  }
  return(fp);
}

void kill_line(FILE *fp) {
  int c=fgetc(fp);
  while((c != EOF) && (c != '\n'))
    c=fgetc(fp);
}