# define MAXINT 111  // is equal to the maximum number of integral which can be stored for
                     // each clump. It should not be possible that this exceed 111 (which is
                     // the no or ort. extinction groups) 


typedef struct {

int N;    // no of intensities in clump

int count; // counts no of integrals stored for each clump

int **y; // for each clump allowed vectors are stored in matrix y

double *value; 

} st_int;


extern st_int *integral;	// (get_data.c)

double done_integral(int *y, long clump_number); // (integral.c)
void store_integral(int *y, double value, long clump_number); // (integral.c)