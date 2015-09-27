long max_monte_iteration;
char input_file_format;


#if !defined(IN_READ_IN_PARAMETERS)
#define IN_READ_IN_PARAMETERS

//#define EXTSYM_PRODUCTION
#if defined(EXTSYM_PRODUCTION)
static char INPUT_DIRECTORY[] = "";
static char OUTPUT_DIRECTORY[] = "";
#else  // to be able to conveniently run ExtSym from visual studio
static char INPUT_DIRECTORY[] = "input/";
static char OUTPUT_DIRECTORY[] = "output/";
#endif

// Uncomment the #defines below to save peak_ignored.asc or
// uncorrelated_peaks.asc or both

//#define SAVE_PEAK_IGNORED
//#define SAVE_UNCORRELATED_PEAK 0

#endif

