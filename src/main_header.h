# include "data_types.h"
# include "integral.h"


void valid(char string[]);  // (valid.c)

void base_integral_value(store_data *data); // (integral_evaluation.c)
void integral_value(store_data data, store_new_data *new_data, 
										char *symmetry, long clump_number); // (integral_evaluation.c)

int t0(int h, int k, int l, char string[]); // (symmetry.c)


