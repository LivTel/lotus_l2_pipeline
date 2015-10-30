int get_error_description(FILE, int, char *, char *, char *);
int write_error_codes_file_to_header(char [], fitsfile *, int *);
int write_key_to_file(char *, char *, char *, int, char *, char *);

char ERROR_CODES_INITIAL_FILE_WRITE_ACCESS [2]	= "a+";			// r readonly; w overwrite; a+ append
char ERROR_CODES_FILE_WRITE_ACCESS [2]	 	= "a+";			// r readonly; w overwrite; a+ append

char REF_ERROR_CODES_FILE [100];					// variable to hold location of error codes reference file
char ERROR_CODES_FILE [100] 			= "error_codes";	// location of file to write error codes to

