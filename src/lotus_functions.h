#include <stdbool.h>

int calc_least_sq_fit (int, int, double [], double [], double [], double *);
int calculate_cross_correlation_offset (double [], double [], int, int, int *, double *);
int check_file_exists (char []);
int check_key_to_omit (char *, char *, char *, bool *);
int find_centroid_parabolic (double [], int [], int, double [], bool);
int find_peaks (int, double [], int [], int *, int, int, int, int, bool);
int find_peaks_contiguous (int, int, double **, int **, int *, int, int, int, int, int, bool);
int find_time (char []);
int flip_array_dbl (double [], int);
int interpolate (char [], double [], double [], int, double, double, double, double []);
int iterative_sigma_clip (double [], int, float, int [], double, double, double *, double *, int *, bool);
int lsearch_int (int [], int, int);
int median_filter (double [], double [], int, int);
int populate_img_parameters (char [], fitsfile *, int, int *, int *, long [], int *, char []);
int populate_env_variable (char [], char []);
int print_file (char *);
char * strdup(char *);
int write_additional_keys_file_to_header(char [], fitsfile *, char [], int, int *);
int write_additional_key_to_file_dbl (char *, char *, char *, double, char *, char *);
int write_additional_key_to_file_str (char *, char *, char *, char *, char *, char *);

char ADDITIONAL_KEYS_INITIAL_FILE_WRITE_ACCESS [2]	= "w";			// r readonly; w overwrite; a+ append
char ADDITIONAL_KEYS_FILE_WRITE_ACCESS [2]	 	= "a+";			// r readonly; w overwrite; a+ append

char ADDITIONAL_KEYS_FILE [100] 			= "additional_keys";	// location of file to write additional keys to



