#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_sort_double.h>
#include "lotus_error_handling.h"
#include "lotus_functions.h"
#include "lotus_config.h"
#include "lotus_red_reformat.h"
#include "lotus_red_rebin.h"

int main (int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 5) {

		if(populate_env_variable(LORF_BLURB_FILE, "L2_LORF_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(LORF_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -1, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		char *input_f		= strdup(argv[1]);
		char *headers_f		= strdup(argv[2]);
		char *operation		= strdup(argv[3]);
		char *output_f		= strdup(argv[4]);

		// ***********************************************************************
		// Check operation choice is recognised

		if (strcmp(operation, "L1_IMAGE") && strcmp(operation, "LSS_NONSS") && strcmp(operation, "SPEC_NONSS") && strcmp(operation, "SPEC_SS")) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -2, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(headers_f);
			free(operation);
			free(output_f);

			return 1;

		}

		// ***********************************************************************
		// Check input file (ARG 1) exists and get parameters if so

		bool missing_input_file = FALSE;

		fitsfile *input_f_ptr;

		int input_f_maxdim = 2;
		int input_f_status = 0, input_f_bitpix, input_f_naxis;
		long input_f_naxes[2] = {1,1};

		if(check_file_exists(input_f)) {	// file doesn't exist	
	
			missing_input_file = TRUE;

		} else {

			if(!fits_open_file(&input_f_ptr, input_f, READONLY, &input_f_status)) {

				if(!populate_img_parameters(input_f, input_f_ptr, input_f_maxdim, &input_f_bitpix, &input_f_naxis, input_f_naxes, &input_f_status, "INPUT FRAME")) {} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -3, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, input_f_status); 

					free(input_f);
					free(headers_f);
					free(operation);
					free(output_f);

					if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1; 

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -4, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				free(input_f);
				free(headers_f);
				free(operation);
				free(output_f);

				return 1; 

			}

		}

		// ***********************************************************************
		// Open the headers file (ARG 2)

		fitsfile *headers_f_ptr;

		int headers_f_status = 0;

		if(!fits_open_file(&headers_f_ptr, headers_f, READONLY, &headers_f_status)) {} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -5, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, headers_f_status); 

			free(input_f);
			free(headers_f);
			free(operation);
			free(output_f);

			if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Create the output file (ARG 4) if it doesn't already exist, otherwise
		// open

		fitsfile *output_f_ptr;

		int output_f_status = 0;

		if(check_file_exists(output_f)) {	// file doesn't exist

			if (!fits_create_file(&output_f_ptr, output_f, &output_f_status)) {

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -6, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(headers_f);
				free(operation);
				free(output_f);

				if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
				if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 

				return 1; 
			
			}

		} else {

			if(!fits_open_file(&output_f_ptr, output_f, READWRITE, &output_f_status)) {

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -7, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(headers_f);
				free(operation);
				free(output_f);

				if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
				if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Print operation specifics to screen

		int hdunum;

		if(!fits_get_num_hdus(output_f_ptr, &hdunum, &output_f_status)) {} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -8, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			free(input_f);
			free(headers_f);
			free(operation);
			free(output_f);

			if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
			if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
			if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

			return 1; 

		}

		printf("\nOperation:\t%s\n", operation);
		printf("Current HDU:\t%d\n", hdunum+1);

		// SPECIFIC FRAME DATA HANDLING
		// ***********************************************************************

		int ii;

		// ***********************************************************************
		// Write blank extension to output file (ARG 4) if input file (ARG 1) is
		// missing

		if (missing_input_file == TRUE) {

			// CREATE BLANK EXTENSION IN OUTPUT FILE 			   (BLANK)
			// ***********************************************************************
			// 1.	Set output file (ARG 4) parameters

			long output_f_naxes[1] = {1};
			int output_f_bitpix = 16, output_f_naxis = 1;

			// 2.	Create a new IMG extension in output file (ARG 4)

			if(!fits_create_img(output_f_ptr, output_f_bitpix, output_f_naxis, output_f_naxes, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -9, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(headers_f);
				free(operation);
				free(output_f);

				if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
				if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
				if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

				return 1; 

			}

			RETURN_FLAG = 7;
			
		} else if (!strcmp(operation, "L1_IMAGE") || !strcmp(operation, "LSS_NONSS") || !strcmp(operation, "SPEC_NONSS") || !strcmp(operation, "SPEC_SS")) { 

			// COPY DATA FROM INPUT FILE TO OUTPUT FILE 			    (COPY) 
			// ***********************************************************************
			// 1.	Set output file (ARG 4) parameters

			long output_f_naxes[2] = {input_f_naxes[0], input_f_naxes[1]};
			int output_f_bitpix = input_f_bitpix, output_f_naxis = input_f_naxis;

			// 2.	Create a new IMG extension in output file (ARG 4)

			if(!fits_create_img(output_f_ptr, output_f_bitpix, output_f_naxis, output_f_naxes, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -10, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(headers_f);
				free(operation);
				free(output_f);

				if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
				if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
				if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

				return 1; 

			}

			// 3.	Copy values across from primary extension of input file (ARG 1)

			if(!fits_copy_data(input_f_ptr, output_f_ptr, &output_f_status)) {} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -11, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(headers_f);
				free(operation);
				free(output_f);

				if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
				if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
				if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

				return 1; 

			}

		}

		// ***********************************************************************
		// Add EXTNAME keyword

		if (!fits_write_key(output_f_ptr, TSTRING, "EXTNAME", operation, NULL, &output_f_status)) {} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -12, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			free(input_f);
			free(headers_f);
			free(operation);
			free(output_f);

			if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
			if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
			if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

			return 1; 

		}

		// GENERIC HEADERS HANDLING
		// ***********************************************************************

		if (missing_input_file == FALSE) {

			if(populate_env_variable(FITS_KEYS_TO_OMIT, "L2_FITS_KEYS_TO_OMIT_FILE")) {	// this file contains a list of keywords to be omitted from each extension

				RETURN_FLAG = 2;

			} else {

				// ***********************************************************************
				// Copy SCALE keywords from input file (ARG 1) to output file (ARG 4) if
				// they exist

				char card [FLEN_CARD];
				bool found_key;

				int nkeys;

				if(!fits_get_hdrspace(input_f_ptr, &nkeys, NULL, &input_f_status)) {} else {

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -13, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, headers_f_status); 

					free(input_f);
					free(headers_f);
					free(operation);
					free(output_f);

					if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
					if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
					if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

					return 1; 

				} 

				for (ii=0; ii<nkeys; ii++) {             

					found_key = FALSE;

					if(!fits_read_record(input_f_ptr, ii, card, &input_f_status)) {} else {

						write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -14, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
						fits_report_error(stdout, input_f_status); 

						free(input_f);
						free(headers_f);
						free(operation);
						free(output_f);

						if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
						if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
						if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

						return 1; 

					}

					if (fits_get_keyclass(card) == TYP_SCAL_KEY) {				// copy only scale keywords (BZERO/BSCALE)

						check_key_to_omit(FITS_KEYS_TO_OMIT, card, operation, &found_key);

						if (found_key == TRUE) {

							// printf("%d\n", TRUE);	// DEBUG

						} else {

							// printf("%d\t%s\n", fits_get_keyclass(card), card);
							fits_write_record(output_f_ptr, card, &output_f_status);

						}

					}

				}

				// ***********************************************************************
				// Copy headers from headers file (ARG 2) to output file (ARG 4)

				if(!fits_get_hdrspace(headers_f_ptr, &nkeys, NULL, &headers_f_status)) {} else {

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -15, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, headers_f_status); 

					free(input_f);
					free(headers_f);
					free(operation);
					free(output_f);

					if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
					if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
					if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

					return 1; 

				} 

				for (ii=0; ii<nkeys; ii++) {

				// printf("%d\t%d\n", ii, nkeys);	// DEBUG
		
					found_key = FALSE;

					if(!fits_read_record(headers_f_ptr, ii, card, &headers_f_status)) {} else {

						write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -16, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
						fits_report_error(stdout, headers_f_status); 

						free(input_f);
						free(headers_f);
						free(operation);
						free(output_f);

						if (fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
						if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
						if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

						return 1; 

					}

					if (fits_get_keyclass(card) >= TYP_WCS_KEY && fits_get_keyclass(card) != TYP_COMM_KEY) {

						check_key_to_omit(FITS_KEYS_TO_OMIT, card, operation, &found_key);

						if (found_key == TRUE) {

							// printf("%d\n", TRUE);	// DEBUG

						} else {

							// printf("%d\t%s\n", fits_get_keyclass(card), card);	// DEBUG

							fits_write_record(output_f_ptr, card, &output_f_status);

						}

					} // else { printf("%d\t%s\n", fits_get_keyclass(card), card); }	// DEBUG

				// printf("%d\n", &nkeys);	// DEBUG
	
				}

			}

		}

		// SPECIFIC HEADERS HANDLING (ADDITIONAL KEYS, ERROR CODES etc.)
		// ***********************************************************************

		if (!strcmp(operation, "L1_IMAGE")) {	

		} else if (!strcmp(operation, "LSS_NONSS")) {

			if (write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "LSS_CALIBRATION", 2, &output_f_status) == 1) {

				RETURN_FLAG = 3;

			}

			if (write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "STARTDATE", 2, &output_f_status) == 1) {

				RETURN_FLAG = 4;

			}

			if (write_error_codes_file_to_header(ERROR_CODES_FILE, output_f_ptr, &output_f_status) == 1) {

				RETURN_FLAG = 5;

			}

		} else if (!strcmp(operation, "SPEC_NONSS")) {

			if (write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "SPEC_CALIBRATION", 2, &output_f_status) == 1) {

				RETURN_FLAG = 3;
			
			}

			if (write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "STARTDATE", 2, &output_f_status) == 1) {

				RETURN_FLAG = 4;

			}

			if (write_error_codes_file_to_header(ERROR_CODES_FILE, output_f_ptr, &output_f_status) == 1) {

				RETURN_FLAG = 5;

			}	

		} else if (!strcmp(operation, "SPEC_SS")) {

                        if (write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "SPEC_CALIBRATION", 2, &output_f_status) == 1) {

                                RETURN_FLAG = 3;
                        
                        }

                        if (write_additional_keys_file_to_header(ADDITIONAL_KEYS_FILE, output_f_ptr, "STARTDATE", 2, &output_f_status) == 1) {

                                RETURN_FLAG = 4;

                        }

                        if (write_error_codes_file_to_header(ERROR_CODES_FILE, output_f_ptr, &output_f_status) == 1) {

                                RETURN_FLAG = 5;

                        }       

                } 

		// ***********************************************************************
		// Write checksums

	  	if (fits_write_chksum(output_f_ptr, &output_f_status)) {

			RETURN_FLAG = 6;	
	
		}

		// ***********************************************************************
		// Clean up heap memory

		free(input_f);
		free(headers_f);
		free(operation);
		free(output_f);

		// ***********************************************************************
		// Close input file (ARG 1), headers file (ARG 2) and output file (ARG 4)

		if (missing_input_file == FALSE) {

			if (fits_close_file(input_f_ptr, &input_f_status)) { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -17, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error (stdout, input_f_status); 

				if (fits_close_file(headers_f_ptr, &headers_f_status)) fits_report_error (stdout, headers_f_status); 
				if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

				return 1; 

			}

		}

		if (fits_close_file(headers_f_ptr, &headers_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -18, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, headers_f_status); 

			if (fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status); 

			return 1; 

		    }

		if (fits_close_file(output_f_ptr, &output_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", -19, "Status flag for L2 loreformat routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, output_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Return success

		//write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRF", RETURN_FLAG, "Status flag for L2 frreformat routine", ERROR_CODES_FILE_WRITE_ACCESS); // don't want to use this, otherwise consecutive calls to frreformat will result in the headers for the current file being populated with the resulting flag from the previous file

		return 0;

	}

}

