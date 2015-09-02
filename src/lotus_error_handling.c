/************************************************************************

 File:				lotus_error_handling.c
 Last Modified Date:     	02/09/15

************************************************************************/

#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

/************************************************************************

 Function:		get_error_description
 Last Modified Date:    26/01/11
 Purpose:		gets an error description from the ERROR_CODES 
			file given an error value
 Required By:		all
 Additional Notes:	None

************************************************************************/

int get_error_description(FILE *REF_ERROR_CODES_FILE, int fits_key_value, char fits_key [], char verbose_comment [], char severity []) {

	char input_string [300];

	char file_fits_key [81];
	int file_fits_key_value;
	char file_verbose_comment [200];		

	bool found_key = FALSE;		// error checking param
	bool found_code = FALSE;	// error checking param

	while(!feof(REF_ERROR_CODES_FILE)) {

		memset(input_string, '\0', sizeof(char)*300);
		memset(file_fits_key, '\0', sizeof(char)*81);
		memset(file_verbose_comment, '\0', sizeof(char)*200);

		fgets (input_string, 300, REF_ERROR_CODES_FILE);
	
		sscanf(input_string, "%s\t%d\t%[^\n]", file_fits_key, &file_fits_key_value, file_verbose_comment);	// last argument in sscanf means read characters as one string up to newline char

		if (input_string[0] == 'L') {	// check the line begins with a valid key

			if (strcmp(file_fits_key, fits_key) == 0) {

				found_key = TRUE;

				if (file_fits_key_value == fits_key_value) {

					found_code = TRUE;

					strcpy(verbose_comment, file_verbose_comment);

					if (fits_key_value >= 1) {

						strcpy(severity, "WARNING");

					} else if (fits_key_value <= -1) {

						strcpy(severity, "CRITICAL");
	
					} else {

						strcpy(severity, "SUCCESS");

					}

				}

			}

		} 

	}

	if (found_code == FALSE && found_key == TRUE) {			// error check - found key, no code

		printf("\nWARNING:\tUnable to write error code to file. Can't find a key %s with error code %d in the reference file.\n", fits_key, fits_key_value);
		printf("\t\tRoutine may or may not continue depending on the severity of the error.\n\n");
		return 1;

	} else if (found_code == FALSE && found_key == FALSE) {		// error check - no key, no code

		printf("\nWARNING:\tUnable to write error code to file. Key %s doesn't exist in the reference file.\n", fits_key);
		printf("\t\tRoutine may or may not continue depending on the severity of the error.\n\n");
		return 1;

	} else {

		return 0;

	}

}

/************************************************************************

 Function:		write_error_codes_file_to_header
 Last Modified Date:    11/09/14
 Purpose:		writes an error codes file to the header
 Required By:		lotus_red_reformat.c
 Additional Notes:	None

************************************************************************/

int write_error_codes_file_to_header(char ERROR_CODES_FILE [], fitsfile *f_ptr, int *status) {

	FILE *file_ptr;
	file_ptr = fopen(ERROR_CODES_FILE, "r");

	char input_string [300];

	char this_keyname [300];
	int this_int_value;
	char this_comment [300];

	if (file_ptr) {
	
		while(!feof(file_ptr)) {

			memset(input_string, '\0', 300);
			fgets (input_string, 300, file_ptr);

			sscanf(input_string, "%s\t%d\t%[^\n]", this_keyname, &this_int_value, this_comment);

			if (!fits_update_key_lng(f_ptr, this_keyname, this_int_value, this_comment, status)) {
		
			} else {

				return 1;

			}

		}

	} else {

		return 1;

	}

	return 0;

}


/************************************************************************

 Function:		write_key_to_file
 Last Modified Date:    11/09/14
 Purpose:		writes an error code to file
 Required By:		all
 Additional Notes:	None

************************************************************************/

int write_key_to_file(char ERROR_CODES_FILE [], char REF_ERROR_CODES_FILE [], char fits_key [], int fits_key_value, char fits_key_comment [], char ERROR_CODES_FILE_WRITE_ACCESS []) {

	FILE *ERROR_CODES_FILE_ptr;
	ERROR_CODES_FILE_ptr = fopen(ERROR_CODES_FILE, ERROR_CODES_FILE_WRITE_ACCESS);

	FILE *REF_ERROR_CODES_FILE_ptr;
	REF_ERROR_CODES_FILE_ptr = fopen(REF_ERROR_CODES_FILE, "r");

	char verbose_comment [200];
	memset(verbose_comment, '\0', sizeof(char)*200);

	char severity [50];
	memset(severity, '\0', sizeof(char)*50);

	if (ERROR_CODES_FILE_ptr && REF_ERROR_CODES_FILE_ptr) {

		if(!get_error_description(REF_ERROR_CODES_FILE_ptr, fits_key_value, fits_key, verbose_comment, severity)) {

			printf("\nRoutine outcome\n");
			printf("---------------\n");

			fprintf(ERROR_CODES_FILE_ptr, "%s\t%d\t%s\n", fits_key, fits_key_value, fits_key_comment);
			printf("\n%s:\t%s.\n", severity, verbose_comment);
			printf("\t\tWrote error code value %d associated with key %s to file %s with comment \"%s\".\n\n", fits_key_value, fits_key, ERROR_CODES_FILE, fits_key_comment);

			fclose(ERROR_CODES_FILE_ptr);
			return 0;

		} else {

			return 1;

		}

	} else {	// failed to find output error codes file or reference error codes file

		if (!ERROR_CODES_FILE_ptr && REF_ERROR_CODES_FILE_ptr) {

			printf("\nWARNING:\tUnable to write error code to file. File %s doesn't exist.\n\n", ERROR_CODES_FILE);
			fclose(REF_ERROR_CODES_FILE_ptr);
			return 1;

		} else if (ERROR_CODES_FILE_ptr && !REF_ERROR_CODES_FILE_ptr) {

			printf("\nWARNING:\tUnable to write error code to file. File %s doesn't exist.\n\n", REF_ERROR_CODES_FILE);
			fclose(ERROR_CODES_FILE_ptr);
			return 1;

		} else {

			printf("\nWARNING:\tUnable to write error code to file. Files %s and %s don't exist.\n\n", REF_ERROR_CODES_FILE, ERROR_CODES_FILE);
			return 1;

		}


	}

}

