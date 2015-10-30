#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "lotus_error_handling.h"
#include "lotus_functions.h"
#include "lotus_config.h"
#include "lotus_red_rebin.h"
#include "lotus_red_arcfit.h"

#include <gsl/gsl_math.h>

int main (int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 8) {

		if(populate_env_variable(LOR_BLURB_FILE, "L2_LOR_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(LOR_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -1, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters
	
		char *input_f	                = strdup(argv[1]);
		double start_wav		= strtod(argv[2], NULL);
		double end_wav			= strtod(argv[3], NULL);
		char *interpolation_type	= strdup(argv[4]);
		double dispersion		= strtod(argv[5], NULL);
		int conserve_flux		= strtol(argv[6], NULL, 0);
		char *output_f	                = strdup(argv[7]);

		// ***********************************************************************
		// Open input file (ARG 1), get parameters and perform any data format 
                // checks 

		fitsfile *input_f_ptr;

		int input_f_maxdim = 2;
		int input_f_status = 0, input_f_bitpix, input_f_naxis;
		long input_f_naxes [2] = {1,1};

		if(!fits_open_file(&input_f_ptr, input_f, READONLY, &input_f_status)) {

			if(!populate_img_parameters(input_f, input_f_ptr, input_f_maxdim, &input_f_bitpix, &input_f_naxis, input_f_naxes, &input_f_status, "INPUT FRAME")) {

				if (input_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -2, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(input_f);
					free(interpolation_type);
					free(output_f);

					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -3, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status); 

				free(input_f);
				free(interpolation_type);
				free(output_f);

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -4, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, input_f_status); 

			free(input_f);
			free(interpolation_type);
			free(output_f);

			return 1; 

		}

		// ***********************************************************************
		// Set the range limits using input fits file (ARG 1)

		int cut_x [2] = {1, input_f_naxes[0]};
		int cut_y [2] = {1, input_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from input fits file (ARG 1)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from input fits file (ARG 1)

		double input_f_pixels [nxelements];

		// ***********************************************************************
		// Open [LOARCFIT_OUTPUTF_WAVFITS_FILE] dispersion solutions file

		FILE *dispersion_solutions_f;
	
		if (!check_file_exists(LOARCFIT_OUTPUTF_WAVFITS_FILE)) { 

			dispersion_solutions_f = fopen(LOARCFIT_OUTPUTF_WAVFITS_FILE , "r");

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -5, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(interpolation_type);
			free(output_f);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1;

		}	

		// ***********************************************************************
		// Find some [LOARCFIT_OUTPUTF_WAVFITS_FILE] file details

		char input_string [500];

		bool find_polynomialorder_comment = FALSE;

		int polynomial_order;	

		char search_string_1 [20] = "# Polynomial Order:\0";	// this is the comment to be found from the [LOARCFIT_OUTPUTF_WAVFITS_FILE] file

		while(!feof(dispersion_solutions_f)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, dispersion_solutions_f);	

			if (strncmp(input_string, search_string_1, strlen(search_string_1)) == 0) { 

				sscanf(input_string, "%*[^\t]%d", &polynomial_order);		// read all data up to tab as string ([^\t]), but do not store (*)
				find_polynomialorder_comment = TRUE;
				break;


			} 

		}

		if (find_polynomialorder_comment == FALSE) {	// error check - didn't find the comment in the [LOARCFIT_OUTPUTF_WAVFITS_FILE] file

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -6, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(interpolation_type);
			free(output_f);

			fclose(dispersion_solutions_f);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1;

		}

		// ***********************************************************************
		// Rewind and extract coefficients from [LOARCFIT_OUTPUTF_WAVFITS_FILE]
		// file 

		rewind(dispersion_solutions_f);

		int token_index;	// this variable will hold which token we're dealing with
		int coeff_index;	// this variable will hold which coefficient we're dealing with
		double this_coeff;
		double this_chisquared;
	
		char *token;

		double coeffs [polynomial_order+1];
		memset(coeffs, 0, sizeof(double)*(polynomial_order+1));

		while(!feof(dispersion_solutions_f)) {

			memset(input_string, '\0', sizeof(char)*500);
	
			fgets(input_string, 500, dispersion_solutions_f);

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
                                  
                                        if (token_index == 0) {    
                                        } else if ((token_index >= 1) && (token_index <= polynomial_order+1)) { 	// coeff token

						this_coeff = strtod(token, NULL);
						// printf("%d\t%e\n", coeff_index, this_coeff);		// DEBUG
						coeffs[coeff_index] = this_coeff;
						coeff_index++;

					} else if (token_index == polynomial_order+2) {					// chisquared token

						this_chisquared = strtod(token, NULL);
                                                //printf("%f\n", this_chisquared);                      // DEBUG

					}

					token_index++;

				}

			}

		}

		// ***********************************************************************
		// Find wavelength extremities from [LOARCFIT_OUTPUTF_WAVFITS_FILE] file
		// and ensure the input constraints [start_wav] (ARG 2) and [end_wav]
		// (ARG 3) don't lie outside these boundaries

		double smallest_wav, largest_wav;

		int ii;

		for (ii=0; ii<=polynomial_order; ii++) {
		    
			smallest_wav += coeffs[ii]*pow(0+INDEXING_CORRECTION, ii);
			largest_wav += coeffs[ii]*pow((cut_x[1]-1)+INDEXING_CORRECTION, ii);

		}

		// ***********************************************************************
	        // Need to find pixel indexes for starting/ending wavelength positions

	        double this_element_wav;
	        int first_element_index, last_element_index;
                int jj;
		for (ii=0; ii<nxelements; ii++) {
                        this_element_wav = 0.0;
                        for (jj=0; jj<=polynomial_order; jj++) {
                                this_element_wav += coeffs[jj]*pow(ii,jj);
                        }
		      	if (this_element_wav >= start_wav) {	// the current index, ii, represents the first pixel with a wavelength >= start_wav. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function
		       		break;
		     	}
		     	
                        first_element_index = ii;		     	
                }

		// printf("%d\t%f\n", ii, this_element_wav);	// DEBUG

		for (ii=nxelements; ii>=0; ii--) {
                        this_element_wav = 0.0;
                        for (jj=0; jj<=polynomial_order; jj++) {
                                this_element_wav += coeffs[jj]*pow(ii,jj);
                        }
                        if (this_element_wav <= end_wav) {    // the current index, ii, represents the first pixel with a wavelength <= end_wav. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function
                                break;
                        }

			last_element_index = ii;

	        }
	        
	        // printf("%d\t%f\n", ii, this_element_wav);     // DEBUG
	
		printf("\nWavelength boundaries");
		printf("\n---------------------\n");

		printf("\nInherent minimum wavelength:\t%.2f Å", smallest_wav);
		printf("\nInherent maximum wavelength:\t%.2f Å\n", largest_wav);

		if (start_wav < smallest_wav) {         // Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function
		  
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -7, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(interpolation_type);
			free(output_f);

			fclose(dispersion_solutions_f);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		} else if (end_wav > largest_wav) {     // Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -8, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(interpolation_type);
			free(output_f);

			fclose(dispersion_solutions_f);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1; 

		}

		// ***********************************************************************
	        // Set the bin wavelengths 

		int num_bins = 0;

		if (!gsl_fcmp((end_wav-start_wav)/dispersion, rint((end_wav-start_wav)/dispersion), 1e-5)) {	// check to see if nearest integer is within tolerance value	

			num_bins = rint((end_wav-start_wav)/dispersion) + 1;					// if TRUE, round

		} else {

			num_bins = floor((end_wav-start_wav)/dispersion) + 1;					// if FALSE, floor

		}

		// printf("%d\n", num_bins);						// DEBUG

		double bin_wavelengths [num_bins];
		memset(bin_wavelengths, 0, sizeof(double)*num_bins);

		for (ii=0; ii<num_bins; ii++) {

			bin_wavelengths[ii] = start_wav + dispersion*ii;		
			// printf("%f\n", bin_wavelengths[ii]);				// DEBUG
		
		}	

		// printf("%f\t%f\n", bin_wavelengths[0], bin_wavelengths[num_bins-1]);	// DEBUG


		// REBIN INPUT FRAME (ARG 1) AND CONSERVE FLUX IF APPLICABLE
		// ***********************************************************************
		// 1.	Open input frame

		int this_row_index;

		double x_wav [nxelements];

		double output_frame_values [nyelements][num_bins];
		memset(output_frame_values, 0, sizeof(double)*nyelements*num_bins);

		double output_f_pixels [num_bins];
		memset(output_f_pixels, 0, sizeof(double)*(num_bins));

		double this_pre_rebin_row_flux, this_post_rebin_row_flux;

		double conservation_factor;	

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			this_row_index = fpixel[1] - 1;

			memset(input_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(input_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, input_f_pixels, NULL, &input_f_status)) {

				// 2.	Calculate pre-rebin total fluxes

				this_pre_rebin_row_flux = 0.0;

				for (ii=first_element_index; ii<=last_element_index; ii++) {

					this_pre_rebin_row_flux += input_f_pixels[ii];

				}

				// 3.	Create pixel-wavelength translation array and perform interpolation

				memset(x_wav, 0, sizeof(double)*nxelements);

				for (ii=0; ii<nxelements; ii++) {

					for (jj=0; jj<=polynomial_order; jj++) {
					    
						x_wav[ii] += coeffs[jj]*pow(ii+INDEXING_CORRECTION,jj);

					}

					// printf("%d\t%f\n", ii, x_wav[ii]); // DEBUG

				}

				// for (ii=0; ii< nxelements; ii++) printf("\n%f\t%f", x_wav[ii], input_f_pixels[ii]);	// DEBUG

				if (interpolate(interpolation_type, x_wav, input_f_pixels, nxelements, bin_wavelengths[0], bin_wavelengths[num_bins-1], dispersion, output_f_pixels)) {
		
					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -9, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(input_f);
					free(interpolation_type);
					free(output_f);

					fclose(dispersion_solutions_f);

					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

					return 1; 

				}

				// 4.	Calculate post-rebin total fluxes

				this_post_rebin_row_flux = 0.0;

				for (ii=0; ii<num_bins; ii++) {

					this_post_rebin_row_flux += output_f_pixels[ii];

				}

				// 5.	Conserve flux if applicable

				conservation_factor = this_pre_rebin_row_flux/this_post_rebin_row_flux;

				// printf("%f\t%f\t%f\n", this_pre_rebin_row_flux, this_post_rebin_row_flux, conservation_factor);	// DEBUG

				for (ii=0; ii<num_bins; ii++) {

					if (conserve_flux == TRUE) {

						output_frame_values[this_row_index][ii] = output_f_pixels[ii]*conservation_factor;

					} else {

						output_frame_values[this_row_index][ii] = output_f_pixels[ii];

					}

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -10, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, input_f_status);

				free(input_f);
				free(interpolation_type);
				free(output_f);

				fclose(dispersion_solutions_f);

				if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

				return 1; 

			}

		}

		// 6.	Create [LOREBIN_OUTPUTF_REBIN_WAVFITS_FILE] output file and print
		// 	a few parameters

		FILE *outputfile;
		outputfile = fopen(LOREBIN_OUTPUTF_REBIN_WAVFITS_FILE, FILE_WRITE_ACCESS);

		if (!outputfile) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -11, "Status flag for L2 frrebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(input_f);
			free(interpolation_type);
			free(output_f);

			fclose(dispersion_solutions_f);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 

			return 1;

		}

		char timestr [80];
		memset(timestr, '\0', sizeof(char)*80);

		find_time(timestr);

		fprintf(outputfile, "#### %s ####\n\n", LOREBIN_OUTPUTF_REBIN_WAVFITS_FILE);
	        fprintf(outputfile, "# Rebinning wavelength fit parameters.\n\n");
                fprintf(outputfile, "# Run Datetime:\t\t%s\n\n", timestr);
	        fprintf(outputfile, "# Target Filename:\t%s\n\n", input_f);
	        fprintf(outputfile, "# Starting Wavelength:\t%.2f\n", bin_wavelengths[0]);
	        fprintf(outputfile, "# Dispersion:\t\t%.2f\n", dispersion);
		fprintf(outputfile, "%d", EOF);

		// 7.	Write these values to the [ADDITIONAL_KEYS_FILE] file

		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CTYPE1", "Wavelength", "Type of co-ordinate on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CUNIT1", "Angstroms", "Units for axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CRVAL1", bin_wavelengths[0], "[pixel] Value at ref. pixel on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CDELT1", dispersion, "[pixel] Pixel scale on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CRPIX1", 1.0, "[pixel] Reference pixel on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CTYPE2", "a2", "Type of co-ordinate on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CUNIT2", "Pixels", "Units for axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CRVAL2", 1, "[pixel] Value at ref. pixel on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CDELT2", 1, "[pixel] Pixel scale on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "LSS_CALIBRATION", "CRPIX2", 1, "[pixel] Reference pixel on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);

		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CTYPE1", "Wavelength", "Type of co-ordinate on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CUNIT1", "Angstroms", "Units for axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CRVAL1", bin_wavelengths[0], "[pixel] Value at ref. pixel on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CDELT1", dispersion, "[pixel] Pixel scale on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CRPIX1", 1.0, "[pixel] Reference pixel on axis 1", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CTYPE2", "a2", "Type of co-ordinate on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_str(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CUNIT2", "Pixels", "Units for axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CRVAL2", 1, "[pixel] Value at ref. pixel on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CDELT2", 1, "[pixel] Pixel scale on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);
		write_additional_key_to_file_dbl(ADDITIONAL_KEYS_FILE, "SPEC_CALIBRATION", "CRPIX2", 1, "[pixel] Reference pixel on axis 2", ADDITIONAL_KEYS_FILE_WRITE_ACCESS);

		// ***********************************************************************
		// Set output frame parameters

		fitsfile *output_f_ptr;
	
		int output_f_status = 0;
		long output_f_naxes [2] = {num_bins,nyelements};
	
		long output_f_fpixel = 1;

		// ***********************************************************************
		// Create [output_frame_values_1D] array to hold the output data in the 
                // correct format

		double output_frame_values_1D [num_bins*nyelements];
		memset(output_frame_values_1D, 0, sizeof(double)*num_bins*nyelements);
                int kk;
		for (ii=0; ii<nyelements; ii++) {
	
			jj = ii * num_bins;
	
			for (kk=0; kk<num_bins; kk++) {
	
				output_frame_values_1D[jj] = output_frame_values[ii][kk];
				jj++;

			}
		
		}

		// ***********************************************************************
		// Create and write [output_frame_values_1D] to output file (ARG 5)	
	
		if (!fits_create_file(&output_f_ptr, output_f, &output_f_status)) {
	
			if (!fits_create_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, output_f_naxes, &output_f_status)) {

				if (!fits_write_img(output_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], output_f_fpixel, num_bins * nyelements, output_frame_values_1D, &output_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -12, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
					fits_report_error(stdout, output_f_status); 

					free(input_f);
					free(interpolation_type);
					free(output_f);

					fclose(dispersion_solutions_f);
					fclose(outputfile);

					if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
					if(fits_close_file(output_f_ptr, &output_f_status)); 

					return 1; 

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -13, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, output_f_status); 

				free(input_f);
				free(interpolation_type);
				free(output_f);

				fclose(dispersion_solutions_f);
				fclose(outputfile);

                                if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
                                if(fits_close_file(output_f_ptr, &output_f_status)); 

				return 1; 

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -14, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, output_f_status); 

			free(input_f);
			free(interpolation_type);
			free(output_f);

			fclose(dispersion_solutions_f);
			fclose(outputfile);

			if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status);  

			return 1; 

		}

		// ***********************************************************************
		// Clean up heap memory

		free(input_f);
		free(interpolation_type);
		free(output_f);

		// ***********************************************************************
		// Close input file (ARG 1), output file (ARG 7) and 
		// [FRARCFIT_OUTPUTF_WAVFITS_FILE] file

		if (fclose(dispersion_solutions_f)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -15, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

			fclose(outputfile);

                        if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
                        if(fits_close_file(output_f_ptr, &output_f_status)); 

			return 1; 

		}

		if (fclose(outputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -16, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

                        if(fits_close_file(input_f_ptr, &input_f_status)) fits_report_error (stdout, input_f_status); 
                        if(fits_close_file(output_f_ptr, &output_f_status)); 

			return 1; 

		}

		if(fits_close_file(input_f_ptr, &input_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -17, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, input_f_status); 

                        if(fits_close_file(output_f_ptr, &output_f_status)); 

			return 1; 

	    	}

		if(fits_close_file(output_f_ptr, &output_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", -18, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, output_f_status); 

			return 1; 

	    	}

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATRE", RETURN_FLAG, "Status flag for L2 lorebin routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}

