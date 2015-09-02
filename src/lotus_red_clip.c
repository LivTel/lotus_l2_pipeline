/************************************************************************

 File:				lotus_red_clip.c
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
#include "lotus_red_clip.h"

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>

// *********************************************************************

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 13) {

		if(populate_env_variable(LOC_BLURB_FILE, "L2_LOC_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(LOC_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -1, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

		return 1;

	} else {

		// ***********************************************************************
		// Redefine routine input parameters

		char *cont_f			= strdup(argv[1]);
		char *in_f			= strdup(argv[2]);		
		int bin_size_px			= strtol(argv[3], NULL, 0);
		double bg_percentile		= strtod(argv[4], NULL);		
		double clip_sigma		= strtod(argv[5], NULL);
		double thresh_sigma		= strtod(argv[6], NULL);		
		int scan_window_size_px		= strtol(argv[7], NULL, 0);
		int scan_window_nsigma		= strtol(argv[8], NULL, 0);
		int min_spectrum_width_px	= strtol(argv[9], NULL, 0);
                int force_lo_px                 = strtol(argv[10], NULL, 0);
                int force_hi_px                 = strtol(argv[11], NULL, 0);
		char *out_f			= strdup(argv[12]);
		
		// ***********************************************************************
		// Open cont file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *cont_f_ptr;

		int cont_f_maxdim = 2, cont_f_status = 0, cont_f_bitpix, cont_f_naxis;
		long cont_f_naxes [2] = {1,1};

		if(!fits_open_file(&cont_f_ptr, cont_f, IMG_READ_ACCURACY, &cont_f_status)) {

			if(!populate_img_parameters(cont_f, cont_f_ptr, cont_f_maxdim, &cont_f_bitpix, &cont_f_naxis, cont_f_naxes, &cont_f_status, "CONTINUUM FRAME")) {

				if (cont_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -2, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

					free(cont_f);
					free(in_f);					
					free(out_f);
					if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -3, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cont_f_status); 

				free(cont_f);
				free(in_f);					
				free(out_f);
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -4, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
			fits_report_error(stdout, cont_f_status); 

			free(cont_f);
			free(in_f);					
			free(out_f);
			
			return 1; 

		}
		
		// ***********************************************************************
		// Open input file (ARG 2), get parameters and perform any data format 
		// checks

		fitsfile *in_f_ptr;

		int in_f_maxdim = 2, in_f_status = 0, in_f_bitpix, in_f_naxis;
		long in_f_naxes [2] = {1,1};

		if(!fits_open_file(&in_f_ptr, in_f, IMG_READ_ACCURACY, &in_f_status)) {

			if(!populate_img_parameters(in_f, in_f_ptr, in_f_maxdim, &in_f_bitpix, &in_f_naxis, in_f_naxes, &in_f_status, "INPUT FRAME")) {

				if (in_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -5, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

					free(cont_f);
					free(in_f);					
					free(out_f);
					if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 
					if(fits_close_file(in_f_ptr, &in_f_status)) fits_report_error (stdout, in_f_status);					

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -6, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
				fits_report_error(stdout, in_f_status); 

				free(cont_f);
				free(in_f);					
				free(out_f);
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 
				if(fits_close_file(in_f_ptr, &in_f_status)) fits_report_error (stdout, in_f_status);

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -7, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
			fits_report_error(stdout, in_f_status); 

			free(cont_f);
			free(in_f);					
			free(out_f);
			
			return 1; 

		}	
		
		// ***********************************************************************
		// Check consistency of arc/target/continuum fits files (ARGS 1, 2 and 3)
	
		printf("\nConsistency check");
		printf("\n-----------------\n");

		printf("\nBits per pixel:\t\t");

		if (cont_f_bitpix != in_f_bitpix) { 	// if a = b and b = c then a must = c

			printf("FAIL\n"); 
                        RETURN_FLAG = 3;

		} else { 

			printf("OK\n"); 

		} 

		printf("Number of axes:\t\t");

		if (cont_f_naxis != in_f_naxis) {	
	
			printf("FAIL\n"); 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -9, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			free(cont_f);
			free(in_f);					
			free(out_f);
			
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 
			if(fits_close_file(in_f_ptr, &in_f_status)) fits_report_error (stdout, in_f_status);
			
			return 1; 

		} else { 

			printf("OK\n"); 

		} 
	
		printf("First axis dimension:\t");

		if (cont_f_naxes[0] != in_f_naxes[0]) {	
	
			printf("FAIL\n"); 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -10, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			free(cont_f);
			free(in_f);					
			free(out_f);
			
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 
			if(fits_close_file(in_f_ptr, &in_f_status)) fits_report_error (stdout, in_f_status);

			return 1; 

		} else { 

			printf("OK\n"); 

		} 

		printf("Second axis dimension:\t");

		if (cont_f_naxes[1] != in_f_naxes[1]) {

			printf("FAIL\n"); 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -11, "Status flag for L2 loclip routinee", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			free(cont_f);
			free(in_f);					
			free(out_f);
			
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 
			if(fits_close_file(in_f_ptr, &in_f_status)) fits_report_error (stdout, in_f_status);
			
			return 1; 

		} else { 

			printf("OK\n");
			

		}				

		// ***********************************************************************
		// Set the range limits n.b. this should be an arbitrary choice if all 
		// files have identical parameters

		int cut_x [2] = {1, cont_f_naxes[0]};
		int cut_y [2] = {1, cont_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from continuum and input fits 
		// file (ARGS 1 and 2)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from continuum and input fits file 
		// (ARGS 1 and 2)

		double cont_f_pixels [nxelements];
		double in_f_pixels [nxelements];

		// ***********************************************************************
		// Get continuum fits file (ARG 1) values and store in 2D array

		int ii;

		double cont_frame_values [nyelements][nxelements];
		memset(cont_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(cont_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(cont_f_ptr, TDOUBLE, fpixel, nxelements, NULL, cont_f_pixels, NULL, &cont_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					cont_frame_values[fpixel[1]-1][ii] = cont_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -12, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cont_f_status); 

				free(cont_f);
				free(in_f);					
				free(out_f);				
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

				return 1; 

			}

		}
		
		// ***********************************************************************
		// Get input fits file (ARG 2) values and store in 2D array

		double in_frame_values [nyelements][nxelements];
		memset(in_frame_values, 0, sizeof(double)*nxelements*nyelements);

		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(in_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(in_f_ptr, TDOUBLE, fpixel, nxelements, NULL, in_f_pixels, NULL, &in_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					in_frame_values[fpixel[1]-1][ii] = in_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -13, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
				fits_report_error(stdout, cont_f_status); 

				free(cont_f);
				free(in_f);					
				free(out_f);				
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

				return 1; 

			}

		}
				
		
		// BIN ARRAY AND FIND BACKGROUND
		// ***********************************************************************		
		// 1.	Bin array according to bin width given by [bin_size_px]
			
		int disp_nelements = nxelements, spat_nelements = nyelements;
		
		int disp_nelements_binned = (int)floor(disp_nelements/bin_size_px);			
		double this_frame_values_binned[spat_nelements][disp_nelements_binned];
		memset(this_frame_values_binned, 0, sizeof(double)*spat_nelements*disp_nelements_binned);
		
		double this_bin_value;
		int bin_number = 0;
		int jj;
		for (jj=0; jj<spat_nelements; jj++) {
			this_bin_value = 0;
			bin_number = 0;
			for (ii=0; ii<disp_nelements; ii++) {
				if (ii % bin_size_px == 0 && ii != 0) {
					this_frame_values_binned[jj][bin_number] = this_bin_value;
					bin_number++;
					this_bin_value = 0;
				}
				this_bin_value += cont_frame_values[jj][ii];
			}
		}	
		
		// 2.	Transform to 1D array
		double this_frame_values_binned_1D[spat_nelements*disp_nelements_binned];
		memset(this_frame_values_binned_1D, 0, sizeof(double)*spat_nelements*disp_nelements_binned);
		
		int idx = 0;
		for (jj=0; jj<spat_nelements; jj++) {
			for (ii=0; ii<disp_nelements_binned; ii++) {
				this_frame_values_binned_1D[idx] = this_frame_values_binned[jj][ii];
				idx++;
			}
		}	

		double this_frame_values_binned_1D_sorted[spat_nelements*disp_nelements_binned];
		memcpy(this_frame_values_binned_1D_sorted, this_frame_values_binned_1D, sizeof(double)*spat_nelements*disp_nelements_binned);	
		gsl_sort(this_frame_values_binned_1D_sorted, 1, spat_nelements*disp_nelements_binned);
			
		int bg_nelements = (int)floor(spat_nelements*disp_nelements_binned*bg_percentile);
		double bg_values [bg_nelements];
		idx = 0;
		for (jj=0; jj<spat_nelements*disp_nelements_binned; jj++) {
			bg_values[idx] = this_frame_values_binned_1D_sorted[jj];
			idx++;
			if (idx == bg_nelements)
					break;
		}
		
		/*if (idx != bg_nelements) {
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -14, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			free(cont_f);
			free(in_f);					
			free(out_f);				
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

			return 1;
		}*/		
		
		double start_mean = gsl_stats_mean(bg_values, 1, bg_nelements);
		double start_sd	  = gsl_stats_sd(bg_values, 1, bg_nelements);		
	
		// 3.	Perform iterative sigma clip of frame to ascertain approximate background level of frame
		double this_frame_values_binned_mean = gsl_stats_mean(this_frame_values_binned_1D, 1, disp_nelements_binned*spat_nelements);
		double this_bg_values_mean, this_bg_values_sd;
		int final_num_retained_indexes;
		int retain_indexes[spat_nelements*disp_nelements_binned]; 
		
		printf("\nBinned background level determination");
		printf("\n-------------------------------------\n");
		
		printf("\nStart mean:\t\t\t%f", start_mean);
		printf("\nStart SD:\t\t\t%f", start_sd);
		iterative_sigma_clip(this_frame_values_binned_1D, spat_nelements*disp_nelements_binned, clip_sigma, retain_indexes, start_mean, start_sd, &this_bg_values_mean, &this_bg_values_sd, &final_num_retained_indexes, TRUE);
		printf("\nFinal mean:\t\t\t%f", this_bg_values_mean);
		printf("\nFinal SD:\t\t\t%f\n", this_bg_values_sd);		
		
		printf("\nFrame Mean (counts):\t\t%.2f", this_frame_values_binned_mean);		
		printf("\nBG Mean (counts):\t\t%.2f", this_bg_values_mean);
		printf("\nBG SD  (counts):\t\t%.2f", this_bg_values_sd);
		printf("\nWindow thresh (counts):\t\t> %.2f\n", this_bg_values_mean + scan_window_nsigma*this_bg_values_sd);
		
		// 4.	Check that there's a significant difference between BG and frame means.
		if (this_bg_values_mean + this_bg_values_sd*thresh_sigma >= this_frame_values_binned_mean) {
			RETURN_FLAG = 2;	
		}
		

		// FIND EDGES OF SPECTRUM
		// ***********************************************************************
		// 1.	Cycle through pixel value array
		
		double this_spat_values[spat_nelements];
		double this_spat_values_der[spat_nelements-1];
		int next_el_spat, this_el_disp, this_el_spat;
		
		double edges_t[disp_nelements_binned];
		double edges_b[disp_nelements_binned];
		int nedges_t = 0, nedges_b = 0;
		int kk;
		for (jj=0; jj<disp_nelements_binned; jj++) {
		
			memset(this_spat_values, 0, sizeof(double)*spat_nelements);
			memset(this_spat_values_der, 0, sizeof(double)*spat_nelements-1);						
			
			// 2.	Accumulate values and derivatives for this bin
			for (ii=0; ii<spat_nelements; ii++) {
				this_el_disp = jj;
				this_el_spat = ii;
				next_el_spat = ii+1;
				this_spat_values[ii] = this_frame_values_binned[this_el_spat][this_el_disp];
				if (ii != spat_nelements-1)
					this_spat_values_der[ii] = this_frame_values_binned[next_el_spat][this_el_disp] - this_frame_values_binned[this_el_spat][this_el_disp];
			}
			
			// 3.	Get index of sorted derivatives (in ascending order)
			size_t this_spat_values_der_idx [spat_nelements-1];	
			gsl_sort_index(this_spat_values_der_idx, this_spat_values_der, 1, spat_nelements-1);	
			
			// 4.	FORWARD DIRECTION (find start of spectrum)
			//	Starting with the highest valued derivative, check that all succesive pixels values within 
			// 	[scan_window_size_px] are higher than [this_spat_values_sd] * [scan_window_nsigma]. 
			//	If true, flag as edge and break. If false, proceed to next highest derivative.
			for (ii=spat_nelements-2; ii>=0; ii--) {
				int this_pk_idx = this_spat_values_der_idx[ii] + 1;		// +1 as we want the pixel at the higher end of the derivative calculation	
				if (this_pk_idx + scan_window_size_px >= spat_nelements) {	// check to make sure window doesn't fall off edge of CCD
					continue;
				} else {
					int is_false_trigger = false;
					for (kk=this_pk_idx; kk<this_pk_idx+scan_window_size_px; kk++) {		// check pixels in window are greater than background
						if (fabs(this_spat_values[kk] - this_bg_values_mean) <= this_bg_values_sd*scan_window_nsigma) {
							is_false_trigger = true;
							break;
						}
					}	
					if (!is_false_trigger) {
						edges_b[nedges_b] = (double)this_pk_idx;
						nedges_b++;
						break;
					}
				}			
			}	
			
			// 5.	REVERSE DIRECTION (find end of spectrum)
			//	Starting with the lowest valued derivative, check that all succesive pixels values within 
			// 	[scan_window_size_px] are higher than [this_spat_values_sd] * [scan_window_nsigma]. 
			//	If true, flag as edge and break. If false, proceed to next lowest derivative.
			for (ii=0; ii<spat_nelements-1; ii++) {
				int this_pk_idx = this_spat_values_der_idx[ii] - 1;		// -1 as we want the pixel at the higher end of the derivative calculation	
				if (this_pk_idx - scan_window_size_px < 0) {			// check to make sure window doesn't fall off edge of CCD
					continue;
				} else {
					int is_false_trigger = false;
					for (kk=this_pk_idx; kk>this_pk_idx+scan_window_size_px; kk--) {	// check pixels in window are greater than background
						if (fabs(this_spat_values[kk] - this_bg_values_mean) <= this_bg_values_sd*scan_window_nsigma) {
							is_false_trigger = true;
							break;
						}
					}	
					if (!is_false_trigger) {
						edges_t[nedges_t] = (double)this_pk_idx;
						nedges_t++;
						break;
					}
				}			
			}				
			
		}
		
		// 6.	Take median of both edge arrays
		double edges_b_sorted [nedges_b];
		double edges_t_sorted [nedges_t];
		memcpy(edges_b_sorted, edges_b, sizeof(double)*nedges_b);		
		memcpy(edges_t_sorted, edges_t, sizeof(double)*nedges_t);	
		gsl_sort(edges_b_sorted, 1, nedges_b);
		gsl_sort(edges_t_sorted, 1, nedges_t);
				
		int median_edges_b = floor(gsl_stats_median_from_sorted_data(edges_b_sorted, 1, nedges_b));
		int median_edges_t = ceil(gsl_stats_median_from_sorted_data(edges_t_sorted, 1, nedges_t));
                
                if (force_lo_px != -1)
                    median_edges_b = force_lo_px;
                if (force_hi_px != -1)
                    median_edges_t = force_hi_px;
                
		int spectrum_width = (int)fabs(median_edges_t - median_edges_b);
		
		printf("\nEdges detected");
		printf("\n--------------\n");
		printf("\nMedian bottom position (px):\t%d", median_edges_b);
		printf("\nMedian top position (px):\t%d", median_edges_t);	
		printf("\nWidth (px):\t\t\t%d\n", spectrum_width);
		
		if (spectrum_width < min_spectrum_width_px) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -14, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			free(cont_f);
			free(in_f);					
			free(out_f);			
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

			return 1; 		

		}
		
		// ***********************************************************************
		// Create [out_frame_values] array to hold the output data in the correct 
		// format
		double out_frame_values [(disp_nelements * (spectrum_width + median_edges_b)) + 1];
		memset(out_frame_values, 0, sizeof(double)*(disp_nelements *  (spectrum_width + median_edges_b)) + 1);
		for (jj=median_edges_b; jj<=spectrum_width+median_edges_b; jj++) {
	
			ii = (jj-median_edges_b) * disp_nelements;
	
			for (kk=0; kk<disp_nelements; kk++) {
	
				out_frame_values[ii] = in_frame_values[jj][kk];
				ii++;
	
			}
	
		}	
		
		// ***********************************************************************
		// Set output frame (ARG 9) parameters

		fitsfile *out_f_ptr;

		int out_f_status = 0;
		long out_f_naxes [2] = {disp_nelements, spectrum_width};
		long out_f_fpixel = 1;		
		
		// ***********************************************************************
		// Create and write trimmed file to output file (ARG 9)
	
		if (!fits_create_file(&out_f_ptr, out_f, &out_f_status)) {
	
			if (!fits_create_img(out_f_ptr, INTERMEDIATE_IMG_ACCURACY[0], 2, out_f_naxes, &out_f_status)) {

				if (!fits_write_img(out_f_ptr, INTERMEDIATE_IMG_ACCURACY[1], out_f_fpixel, (disp_nelements * spectrum_width) + 1, out_frame_values, &out_f_status)) {

				} else { 

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -15, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

					free(cont_f);
					free(in_f);					
					free(out_f);			
					if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 
					if(fits_close_file(out_f_ptr, &out_f_status)) fits_report_error (stdout, out_f_status); 					

					return 1; 	

				}

			} else {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -16, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

				free(cont_f);
				free(in_f);					
				free(out_f);			
				if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 
				if(fits_close_file(out_f_ptr, &out_f_status)) fits_report_error (stdout, out_f_status); 
				
				return 1; 	

			}

		} else {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -17, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

			free(cont_f);
			free(in_f);					
			free(out_f);			
			if(fits_close_file(cont_f_ptr, &cont_f_status)) fits_report_error (stdout, cont_f_status); 

			return 1; 	

		}	

		// ***********************************************************************
		// Clean up heap memory

		free(cont_f);
		free(in_f);					
		free(out_f);

		// ***********************************************************************
		// Close continuum file (ARG 1), input file (ARG 2) and output file 
		// (ARG 9)

		if(fits_close_file(cont_f_ptr, &cont_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -18, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
			fits_report_error (stdout, cont_f_status); 

			return 1; 

	    	}
	    	
		if(fits_close_file(in_f_ptr, &in_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -19, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
			fits_report_error (stdout, in_f_status); 

			return 1; 

	    	}
	    	
		if(fits_close_file(out_f_ptr, &out_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", -20, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);
			fits_report_error (stdout, out_f_status); 

			return 1; 

	    	}	    	

		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATCL", RETURN_FLAG, "Status flag for L2 loclip routine", ERROR_CODES_INITIAL_FILE_WRITE_ACCESS);

		return 0;

	}

}


