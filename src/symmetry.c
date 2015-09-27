# include <stdio.h>
# include <math.h>

# include "data_types.h"


/************************************************************
 *
 *  symm_hkl -- takes a set of h,k,l values from a clump and sets y[i] to 1 or 2 if
 *                  h,k,l is allowed and 0 if it isn't
 *
 *  Returns:    returns the number of refined peaks which are allowed by the 
 *				input symmetry, string[], and vector y[], where each element
 *				states whether an element of the clump is either absent, singles	
 *				, doublet etc
 *
 ************************************************************/


int t0(int h, int k, int l, char string[]);
int t_special(int h, int k, int l, char string[]);
int t1(int h, int k, int l, char string[]);
int t2(int h, int k, int l, char string[]);
int t3(int h, int k, int l, char string[]);
int t4(int h, int k, int l, char string[]);
int t5(int h, int k, int l, char string[]);
int t6(int h, int k, int l, char string[]);


// ************************************************************
// does the set up and starts off the chain of calls
// ************************************************************

void symm_hkl(char string[], store_data data,  store_new_data *new_data_ptr)
{
	int i, j;								// simple index
	int absent_dummy;					// check if peak absent when doublet
	
	new_data_ptr->any_changes_of_block = NO;
	//	new_data_ptr->num_delta_priors = 0;
	new_data_ptr->N = 0;			// number of allowed peaks
	
	
	for(i=1;i<=data.N;i++)
	{
		new_data_ptr->allowed[i]  = t0(data.h[i], data.k[i], data.l[i] ,string);
				
		// check if peak absent
		
		if (new_data_ptr->allowed[i] == 0) {  // if peak absent then
			new_data_ptr->any_changes_of_block = YES; // block changes and
			//			++(new_data_ptr->num_delta_priors); // one more delta function is present in block
		}
		
			
		// generate new_data.h,k,l for displaying purposes
		
		if (new_data_ptr->allowed[i] == 1)
		{
			new_data_ptr->h[new_data_ptr->N+1] = data.h[i];
			new_data_ptr->k[new_data_ptr->N+1] = data.k[i];
			new_data_ptr->l[new_data_ptr->N+1] = data.l[i];
		}
		
			
		if(data.multiplet[i] >= 1) 
		{
			for (j = 1; j <= data.multiplet[i]; j++)
			{
				absent_dummy = t0(data.h_multiplet[i][j], data.k_multiplet[i][j],
					data.l_multiplet[i][j], string);
				new_data_ptr->allowed[i] += absent_dummy;
			
				// check if peak absent
			
				if (absent_dummy == 0) {  // if peak absent then
					new_data_ptr->any_changes_of_block = YES; // block changes and
					//				++(new_data_ptr->num_delta_priors); // one more delta function is present in block
				}
			
			
				// generate new_data.h2,k2,l2 for displaying purposes
			
				if (absent_dummy == 1)
				{	
					if (new_data_ptr->allowed[i] >= 2)
					{
						new_data_ptr->h_multiplet[new_data_ptr->N+1][new_data_ptr->allowed[i]-1] = data.h_multiplet[i][j];
						new_data_ptr->k_multiplet[new_data_ptr->N+1][new_data_ptr->allowed[i]-1] = data.k_multiplet[i][j];
						new_data_ptr->l_multiplet[new_data_ptr->N+1][new_data_ptr->allowed[i]-1] = data.l_multiplet[i][j];
					}
					else
					{
						new_data_ptr->h[new_data_ptr->N+1] = data.h_multiplet[i][j];
						new_data_ptr->k[new_data_ptr->N+1] = data.k_multiplet[i][j];
						new_data_ptr->l[new_data_ptr->N+1] = data.l_multiplet[i][j];
					}
				}
			}
		}
				
		
		// if new_data_ptr->allowed[i] non-zero then add one to the number of allowed peaks 
		// in clump
		//
		//		new_data_ptr->allowed[i] = 1 (singlet)     -  peak allowed
		//		new_data_ptr->allowed[i] = 2 (doublet)     -  peak allowed
		//      new_data_ptr->allowed[i] = 0			   -  peak not allowed
		
		if (new_data_ptr->allowed[i]!=0) 
		{
			++(new_data_ptr->N);
			
			
			// make multiplet array 
			// the new_data.multiplet array will of course only be different from 
			// data.doublet array when any_changes_of_block is YES
			
			if (new_data_ptr->allowed[i] >= 2)
			{
				new_data_ptr->multiplet[new_data_ptr->N] = new_data_ptr->allowed[i] - 1;
			}
			else
			{
				new_data_ptr->multiplet[new_data_ptr->N] = 0;
			}
		}	
	}
}


// *******************************************************************
// tn() checks its bit and returns 0 if its bit is not valid, or calls 
// tn+1 which checks the next until t6 is reached. If t6 passes the 
// values then 1 is returned 
// *******************************************************************

int t0(int h, int k, int l, char string[])
{
	
	switch(string[0])
	{
	case 'C':
		if((h+k)%2 != 0) return 0;         // (h+k)%2 = 0 if h+k is even
		break;
	case 'B':
		if((h+l)%2 != 0) return 0;
		break;		
	case 'A':
		if((k+l)%2 != 0) return 0;
		break;
	case 'I':
		if((h+k+l)%2 !=0) return 0;
		break;
	case 'R':   // test for both reverse and obverse hexagonal settings 
		if((h-k+l)%3 != 0 && (-h+k+l)%3 != 0) return 0;
		break;
	case 'F':
		if((h+k)%2 != 0 || (h+l)%2 != 0 || (k+l)%2 != 0) return 0;
		break;
	default:
		break;
	}
	
	// if it passes the check we go onto t1
	
	return t_special(h,k,l,string);
}


//-------------------------------------------------------------------------------
// this functions treat special symmetry absence rules
//-------------------------------------------------------------------------------
int t_special(int h, int k, int l, char string[])
{
	// hexa		cubic		cubic		hexa
	// trigo									trigo
	// cubic
	// tetra
	if(h==k || h==l || k==l || h==-k)
    {
		switch(string[7])
		{
		case 'l':
			if(l%2 != 0 && h==k) return 0;  // hexa, tetra and trigo
			break;
		case '4':
			if(((2*h+l)%4 != 0 || l%2 != 0) && h==k) return 0;  // tetra
			break;
		case 'i':
			if(l%2 != 0 && h==-k) return 0;	// trigo and hexa
			break;
		case 'e':
			if(l%2 != 0 && (h==k || h==-k)) return 0;		// hexa
			break;
		case 'c':
			if(h==k || h==l || k==l) 
			{
				if(h==k)
					if(l%2 != 0) return 0;		// cubic
				if(h==l)
					if(k%2 != 0) return 0;
				if(k==l)
					if(h%2 != 0) return 0;
			}
			break;
		case 'b':
			if(h==k || h==l || k==l) 
			{
				if(h==k)
					if((2*h+l)%4 != 0 || l%2 != 0) return 0;		// cubic
				if(h==l)
					if((2*h+k)%4 != 0 || k%2 != 0) return 0;
				if(k==l)
					if((2*k+h)%4 != 0 || h%2 != 0) return 0;
			}
			break;
    default:
			break;
		}
		return t1(h,k,l,string);
    }
	else	return t1(h,k,l,string);	
}



//-------------------------------------------------------------------------------
// 0kl
//-------------------------------------------------------------------------------
int t1(int h, int k, int l, char string[])
{
	if(h==0)
    {
		switch(string[1])
		{
		case 'l':
			if(l%2 != 0) return 0;
			break;
		case 'k':
			if(k%2 != 0) return 0;
			break;
		case 's':
			if((k+l)%2 != 0) return 0;
			break;
		case '4':
			if(k%2 != 0 || (k+l)%4 != 0) return 0;  // special rule k+l=4n;k,l
			break;
        default:
			break;
		}
		return t2(h,k,l,string);
    }
	else	return t2(h,k,l,string);	
}


//-------------------------------------------------------------------------------
// h0l
//-------------------------------------------------------------------------------
int t2(int h, int k, int l, char string[])
{
	if(k==0)
    {
		switch(string[2])
		{
		case 'l':
			if(l%2 != 0) return 0;
			break;
		case 'h':
			if(h%2 != 0) return 0;
			break;
		case 's':
			if((h+l)%2 != 0) return 0;
			break;
		case '4':
			if(h%2 != 0 || (h+l)%4 != 0) return 0;	// special rule h+l=4n;h,l
			break;
        default:
			break;
		}
		return t3(h,k,l,string);
    }
	else	return t3(h,k,l,string);	
}


//-------------------------------------------------------------------------------
// hk0
//-------------------------------------------------------------------------------
int t3(int h, int k, int l, char string[])
{
	if(l==0)
    {
		switch(string[3])
		{
		case 'h':
			if(h%2 != 0) return 0;
			break;
		case 'k':
			if(k%2 != 0) return 0;
			break;
		case 's':
			if((h+k)%2 != 0) return 0;
			break;
		case '4':
			if(h%2 != 0 || (h+k)%4 != 0) return 0;	// ortho and cubic
			break;
        default:
			break;
		}
		return t4(h,k,l,string);
    }
	else	return t4(h,k,l,string);	
}


//-------------------------------------------------------------------------------
// h00
//-------------------------------------------------------------------------------
int t4(int h, int k, int l, char string[])
{
	if(k==0 && l==0)
	{
		switch(string[4])
		{
		case 'h':
			if(h%2 != 0) return 0;
			break;
		case '4':    
			if(h%4 != 0) return 0;		// cubic
			break;
		default:
			break;
		}
		return t5(h,k,l,string);
	}	
	else	return t5(h,k,l,string);	
}


//-------------------------------------------------------------------------------
// 0k0
//-------------------------------------------------------------------------------
int t5(int h, int k, int l, char string[])
{
	if(h==0 && l==0)
	{
		switch(string[5])
		{
		case 'k':
			if(k%2 != 0) return 0;
			break;
		case '4':
			if(k%4 != 0) return 0;  // cubic 
			break;
		default:
			break;
		}
		return t6(h,k,l,string);
	}	
	else	return t6(h,k,l,string);	
}


//-------------------------------------------------------------------------------
// 00l
//-------------------------------------------------------------------------------
int t6(int h, int k, int l, char string[])
{
	if(h==0 && k==0)
	{
		switch(string[6])
		{
		case 'l':
			if(l%2 != 0) return 0;
			break;
		case '4':
			if(l%4 != 0) return 0;  // tetra, cubic
			break;
		case '3':
			if(l%3 != 0) return 0;	// trigo, hexa
			break;
		case '6':
			if(l%6 != 0) return 0;	// trigo, hexa
			break;
		default:
			break;
		}
		return 1;
	}
	else	return 1;	
}















