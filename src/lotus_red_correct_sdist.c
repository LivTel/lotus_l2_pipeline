/************************************************************************

 File:				lotus_red_correct_sdist.c
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
#include "lotus_red_correct_sdist.h"
#include "lotus_red_trace_sdist.h"
#include "gsl_poly.h"

// *********************************************************************

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 5) {

		if(populate_env_variable(LOCS_BLURB_FILE, "L2_LOCS_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(LOCS_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -1, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {
		// ***********************************************************************
		// Redefine routine input parameters
		
		char *input_f			= strdup(argv[1]);
		char *interpolation_type	= strdup(argv[2]);
		int conserve_flux		= strtol(argv[3], NULL, 0);		
		char *output_f			= strdup(argv[4]);		
		
		// ***********************************************************************
		// Open input file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *input_f_ptr;

		int input_f_maxdim = 2, input_f_status = 0, input_f_bitpix, input_f_naxis;
		long input_f_naxes [2] = {1,1};

		if(!fits_open_file(&input_f_ptr, input_f, IMG_READ_ACCURACY, &input_f_status)) {

			if(!populate_img_parameters(input_f, input_f_ptr, input_f_maxdim, &input_f_bitpix, &input_f_naxis, input_f_naxes, &input_f_status, "INPUT FRAME")) {

				if (input_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -2, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(input_f);
					free(output_f);					
					free(interpolation_type);
					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -3, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				free(input_f);
					free(output_f);					
				free(interpolation_type);
				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -4, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, input_f_status); 

			free(input_f);
			free(output_f);				
			free(interpolation_type);			

			return 1; 

		}
		
		// ***********************************************************************
		// Set the range limits

		int cut_x [2] = {1, input_f_naxes[0]};
		int cut_y [2] = {1, input_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from input file (ARG 1)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from input fits file (ARG 1)

		double input_f_pixels [nxelements];
		
		// ***********************************************************************
		// Get input fits file (ARG 1) values and store in 2D array

		int ii;

		double input_frame_values [nyelements][nxelements];
		memset(input_frame_values, 0, sizeof(double)*nxelements*nyelements);
		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(input_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(input_f_ptr, TDOUBLE, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					input_frame_values[fpixel[1]-1][ii] = input_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -5, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				free(input_f);
				free(output_f);					
				free(interpolation_type);				
				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		}	
		
		// ***********************************************************************
		// Open [LOTRACE_OUTPUTF_TRACES_FILE] input file
	
		FILE *inputfile;
	
		if (!check_file_exists(LOTRACE_OUTPUTF_TRACES_FILE)) { 

			inputfile = fopen(LOTRACE_OUTPUTF_TRACES_FILE , "r");

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -6, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			return 1;

		}	
		
		// ***********************************************************************
		// Find some [LOTRACE_OUTPUTF_TRACES_FILE] file details

		char input_string [500];

		bool find_polynomialorder_comment = FALSE;

		int polynomial_order;	

		char search_string_1 [20] = "# Polynomial Order:\0";	// this is the comment to be found from the [LOTRACE_OUTPUTF_TRACES_FILE] file

		while(!feof(inputfile)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, inputfile);	

			if (strncmp(input_string, search_string_1, strlen(search_string_1)) == 0) { 

				sscanf(input_string, "%*[^\t]%d", &polynomial_order);		// read all data up to tab as string ([^\t]), but do not store (*)
				find_polynomialorder_comment = TRUE;
				break;

			} 

		}

		if (find_polynomialorder_comment == FALSE) {	// error check - didn't find the comment in the [LOTRACE_OUTPUTF_TRACES_FILE] file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -7, "Status flag for L2 frcorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(output_f);				
			free(interpolation_type);				
			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1;

		}
		
		// ***********************************************************************
		// Rewind and extract coefficients from [LOTRACE_OUTPUTF_TRACES_FILE] file 

		rewind(inputfile);

		int token_index;	// this variable will hold which token we're dealing with
		int coeff_index;	// this variable will hold which coefficient we're dealing with
		double this_coeff;
		double this_chisquared;
	
		char *token;

		double coeffs [polynomial_order+1];
		memset(coeffs, 0, sizeof(double)*(polynomial_order+1));

		while(!feof(inputfile)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, inputfile);

			token_index = 0;
			coeff_index = 0;

			if (strtol(&input_string[0], NULL, 0) > 0) { 		// check the line begins with a positive number
				
				// ***********************************************************************
				// String tokenisation loop: 
				//
				// 1. init calls strtok() loading the function with input_string
				// 2. terminate when token is null
				// 3. we keep assigning tokens of input_string to token until termination by calling strtok with a NULL first argument
				// 
				// n.b. searching for tab or newline separators ('\t' and '\n')
				for (token=strtok(input_string, "\t\n"); token !=NULL; token = strtok(NULL, "\t\n")) {	
					if ((token_index >= 0) && (token_index <= polynomial_order)) { 			// coeff token
						this_coeff = strtod(token, NULL);
						//printf("%d\t%e\n", coeff_index, this_coeff);				// DEBUG
						coeffs[coeff_index] = this_coeff;
						coeff_index++;
					} else if (token_index == polynomial_order+1) {					// chisquared token
						this_chisquared = strtod(token, NULL);
					}
					token_index++;
				}
			}
		}	
		
		// ***********************************************************************
		// Determine the maximum offset from c0 (this is needed to avoid 
		// trying to interpolate past the limits). Use maximum offset as clip
		// for both sides.
		
		double c0 = coeffs[0];
		int max_offset = 0;
		for (ii=0; ii<nxelements; ii++) {
			int this_offset_ceil = abs(ceil(c0 - gsl_poly_eval(coeffs, polynomial_order+1, ii)));
			if (this_offset_ceil > max_offset)
				max_offset = this_offset_ceil;
		}	
		
		int nyelements_reb = nyelements-(2*max_offset);
		
		// ***********************************************************************
		// Do the rebinning (conserving flux where applicable)

		double reb_values[nyelements_reb][nxelements];
		memset(reb_values, 0, sizeof(double)*nyelements_reb*nxelements);		
		
		double this_pre_rebin_flux, this_post_rebin_flux;
		double this_column_values[nyelements];
		double this_column_values_reb[nyelements_reb];		
		double x_offsetted[nyelements];
		
		int jj;
		for (ii=0; ii<nxelements; ii++) {
			this_pre_rebin_flux = 0.;
			double this_offset = c0 - gsl_poly_eval(coeffs, polynomial_order+1, ii);
			memset(this_column_values, 0, sizeof(double)*nyelements);
			memset(this_column_values_reb, 0, sizeof(double)*nyelements);			
			for (jj=0; jj<nyelements; jj++) {
				this_column_values[jj] = input_frame_values[jj][ii];
				x_offsetted[jj] = jj + this_offset;
				this_pre_rebin_flux += input_frame_values[jj][ii];
			}

			interpolate(interpolation_type, x_offsetted, this_column_values, nyelements, max_offset, nyelements-max_offset, 1, this_column_values_reb);

			// get post rebin flux
			this_post_rebin_flux = 0.;
			for (jj=0; jj<nyelements_reb; jj++) {
				this_post_rebin_flux += this_column_values_reb[jj];
			}
			
			// apply conservation factor
			double conservation_factor = this_pre_rebin_flux/this_post_rebin_flux;
			//printf("%f\t%f\t%f\n", this_pre_rebin_flux, this_post_rebin_flux, conservation_factor);	// DEBUG
			if (conserve_flux == TRUE) {
				for (jj=0; jj<nyelements_reb; jj++) {
					reb_values[jj][ii] = this_column_values_reb[jj] * conservation_factor;
				} 
			} else {
				for (jj=0; jj<nyelements_reb; jj++) {				
					reb_values[jj][ii] = this_column_values_reb[jj];					
				}
			}	
		}
		
		// ***********************************************************************
		// Set output frame parameters	

		fitsfile *output_f_ptr;
	
		int output_f_status = 0;
		long output_f_naxes [2] = {nxelements,nyelements_reb};
	
		long output_f_fpixel = 1;

		// ***********************************************************************
		// Create [output_frame_values] array to hold the output data in the 
		// correct format
                
		int kk;
		double output_frame_values [nxelements*nyelements_reb];
		memset(output_frame_values, 0, sizeof(double)*nxelements*nyelements_reb);
		for (ii=0; ii<nyelements_reb; ii++) {	
			jj = ii * nxelements;			
			for (kk=0; kk<nxelements; kk++) {
				output_frame_values[jj] = reb_values[ii][kk];
				jj++;
			}
		}	
		
		// ***********************************************************************
		// Create and write [output_frame_values] to output file (ARG 4)
	
		if (!fits_create_file(&output_f_ptr, output_f, &output_f_status)) {
	
			if (!fits_create_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, output_f_naxes, &output_f_status)) {

				if (!fits_write_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], output_f_fpixel, nxelements*nyelements_reb, output_frame_values, &output_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -8, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, output_f_status); 

					free(input_f);
					free(output_f);				
					free(interpolation_type);
					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
					if(fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status);

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -9, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(output_f);				
				free(interpolation_type);
				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
				if(fits_close_file(output_f_ptr, &output_f_status)) fits_report_error (stdout, output_f_status);

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -10, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			free(input_f);
			free(output_f);				
			free(interpolation_type);
			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		}

		// ***********************************************************************
		// Free arrays on heap

		free(input_f);
		free(interpolation_type);
		free(output_f);
		
		// ***********************************************************************
		// Close [LOTRACE_OUTPUTF_TRACES_FILE] output file, input file (ARG 1) and
		// output file (ARG 4)
		
		if (fclose(inputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -11, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		}		

		if(fits_close_file(input_f_ptr, &input_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -12, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, input_f_status); 

			return 1; 

	    	}
	    	
		if(fits_close_file(output_f_ptr, &output_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", -13, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, output_f_status); 

			return 1; 

	    	}	    	
		
		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCO", RETURN_FLAG, "Status flag for L2 locorrect routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}

