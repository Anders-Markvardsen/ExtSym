# include <stdlib.h>		// (to include exit(1)) 
# include <stdio.h>
# include <math.h>
# include "data_types.h"


/***********************************************************************
*
*  valid -- a routine which check if a seven letter orthorhombic extinction symbol 
*			is valid         
*
*  Usage :	if the symbol is not valid the program is abrupted and an error
*			message is displayed
*
*************************************************************************/

int strcmp_eight_char(char st1[],char st2[]);	// used by valid()

extern struct description *symmetries;  // (initialize_sym_struct.c)

extern int num_extinction_groups; // (initialize_sym_struct.c)

void valid(char st[])
{
	int i;              // simple index
	
	enum yes_no {NO, YES} found_one;	// found_one state if symmestry is valid
	
	
	// initialise string not to be a valid extinction group

	found_one = NO;


	// tests if the string is the same as each of the allowed strings
	// if it find one then set found_one = YES;
	
	for(i = 1; i <= num_extinction_groups; i++)	// nough element of symmetries is blanc
	{
		if(strcmp_eight_char(symmetries[i].symmetry,st)==0) 
		{
			//printf("%s   :    ",symmetries[i].name);
			found_one = YES;
			break;
		}
	}
	

	// if the string is not a valid extinction class then exit program

	if (found_one == NO) {
		printf("\n%s is not a valid extinction class", st);
		printf(".... Now exit program\n");
		exit(1);
	}
}

// strcmp_seven_char is used instead of strcmp incase the strings are not correctly
// initialised. This is probably unneccessary

int strcmp_eight_char(char st1[],char st2[])
{
	int i;
	
	for(i=0;i<=7;i++) if(st1[i]!=st2[i]) return 1;   // returns 1 if they are different
	
	return 0;         // returns 0 if they are the same
}





