/************************************************************************

 File:				lotus_red_!TODO.c
 Last Modified Date:     	02/09/15

************************************************************************/

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "lotus_error_handling.h"
#include "lotus_functions.h"
#include "lotus_config.h"
#include "lotus_red_"/*!TODO*/".h"

// *********************************************************************

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != /*!TODO*/) {

		if(populate_env_variable(LO/*!TODO*/_BLURB_FILE, "L2_LO"/*!TODO*/"_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(SP/*!TODO*/_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STAT"/*!TODO*/"", -1, "Status flag for L2 "/*!TODO*/" routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {
		// ***********************************************************************
		// Redefine routine input parameters
		
		/*!TODO*/
		/*char str			= strdup(argv[1]);	
		int int				= strtol(argv[2], NULL, 0);
		double double			= strtod(argv[3], NULL);*/
		
		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, ""/*!TODO*/"", RETURN_FLAG, "Status flag for L2 "/*!TODO*/" routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}

