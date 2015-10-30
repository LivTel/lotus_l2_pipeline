#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "lotus_error_handling.h"
#include "lotus_functions.h"
#include "lotus_config.h"
#include "lotus_red_find.h"

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>

int main(int argc, char *argv []) {

	if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

		printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

	}

	if (argc != 15) {

		if(populate_env_variable(LOF_BLURB_FILE, "L2_LOF_BLURB_FILE")) {

			RETURN_FLAG = 1;

		} else {

			print_file(LOF_BLURB_FILE);

		}

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -1, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 1;

	} else {
		// ***********************************************************************
		// Redefine routine input parameters
		
		char *target_f				= strdup(argv[1]);	
		int bin_size_px				= strtol(argv[2], NULL, 0);
		double bg_percentile			= strtod(argv[3], NULL);
		double clip_sigma			= strtod(argv[4], NULL);
		int median_filter_width_px		= strtol(argv[5], NULL, 0);	
		double min_SNR				= strtod(argv[6], NULL);
		int min_spatial_width_px		= strtol(argv[7], NULL, 0);
                int finding_window_lo_px                = strtol(argv[8], NULL, 0);
                int finding_window_hi_px                = strtol(argv[9], NULL, 0);                
		int max_centering_num_px		= strtol(argv[10], NULL, 0);		
		int centroid_half_window_size_px	= strtol(argv[11], NULL, 0);
		int min_used_bins			= strtol(argv[12], NULL, 0);
		int window_x_lo				= strtol(argv[13], NULL, 0);
		int window_x_hi				= strtol(argv[14], NULL, 0);
		
		// ***********************************************************************
		// Open target file (ARG 1), get parameters and perform any data format 
		// checks

		fitsfile *target_f_ptr;

		int target_f_maxdim = 2, target_f_status = 0, target_f_bitpix, target_f_naxis;
		long target_f_naxes [2] = {1,1};

		if(!fits_open_file(&target_f_ptr, target_f, IMG_READ_ACCURACY, &target_f_status)) {

			if(!populate_img_parameters(target_f, target_f_ptr, target_f_maxdim, &target_f_bitpix, &target_f_naxis, target_f_naxes, &target_f_status, "TARGET FRAME")) {

				if (target_f_naxis != 2) {	// any data format checks here

					write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -2, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

					free(target_f);
					if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

					return 1;
	
				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -3, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, target_f_status); 

				free(target_f);
				if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

				return 1; 

			}

		} else { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -4, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error(stdout, target_f_status); 

			free(target_f);

			return 1; 

		}
		
		// ***********************************************************************
		// Set the range limits

		int cut_x [2] = {window_x_lo, window_x_hi};
		int cut_y [2] = {1, target_f_naxes[1]};

		// ***********************************************************************
		// Set parameters used when reading data from target file (ARG 1)

		long fpixel [2] = {cut_x[0], cut_y[0]};
		long nxelements = (cut_x[1] - cut_x[0]) + 1;
		long nyelements = (cut_y[1] - cut_y[0]) + 1;

		// ***********************************************************************
		// Create arrays to store pixel values from target fits file (ARG 1)

		double target_f_pixels [nxelements];
		
		// ***********************************************************************
		// Get target fits file (ARG 1) values and store in 2D array

		int ii;

		double target_frame_values [nyelements][nxelements];
		memset(target_frame_values, 0, sizeof(double)*nxelements*nyelements);
		for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

			memset(target_f_pixels, 0, sizeof(double)*nxelements);

			if(!fits_read_pix(target_f_ptr, TDOUBLE, fpixel, nxelements, NULL, target_f_pixels, NULL, &target_f_status)) {

				for (ii=0; ii<nxelements; ii++) {

					target_frame_values[fpixel[1]-1][ii] = target_f_pixels[ii];

				}

			} else { 

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -5, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);
				fits_report_error(stdout, target_f_status); 

				free(target_f);									
				if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

				return 1; 

			}

		}
		
		// FIND VALUES OF PEAK CENTROID ALONG DISPERSION AXIS
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
				this_bin_value += target_frame_values[jj][ii];
			}
		}

		printf("\nFinding peaks");
		printf("\n-------------------------------------\n");	
		double peaks[disp_nelements_binned];
		int num_bins_used = 0;
		for (ii=0; ii<disp_nelements_binned; ii++) {

			// 1a.	Establish if any target flux is in this bin
			// 	First find the mean/sd of the [bg_percentile]th lowest valued pixels as an initial parameters for sigma clip			
			double this_spat_values[spat_nelements];
			double this_spat_values_sorted[spat_nelements];
			for (jj=0; jj<spat_nelements; jj++) {
				this_spat_values[jj] = this_frame_values_binned[jj][ii];
			}			
			memcpy(this_spat_values_sorted, this_spat_values, sizeof(double)*spat_nelements);	
			gsl_sort(this_spat_values_sorted, 1, spat_nelements);
			
			int bg_nelements = (int)floor(spat_nelements*bg_percentile);
			double bg_values [bg_nelements];
			int idx = 0;
			for (jj=0; jj<spat_nelements; jj++) {
				if (this_spat_values_sorted[jj] != 0) {			// avoid 0s set from median filter edges
					bg_values[idx] = this_spat_values_sorted[jj];
					idx++;
					if (idx == bg_nelements)
						break;
				}
			}

			if (idx != bg_nelements) {
				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -6, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

				free(target_f);
				if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

				return 1;
			}
			
			double start_mean = gsl_stats_mean(bg_values, 1, bg_nelements);
			double start_sd	  = gsl_stats_sd(bg_values, 1, bg_nelements);

			// 1b.	Iterative sigma clip around dataset with the initial guess
			int retain_indexes[spat_nelements];
			double final_mean, final_sd;
			int final_num_retained_indexes;
			
			printf("\nBin:\t\t\t\t%d", ii);	
			printf("\nStart mean:\t\t\t%f", start_mean);
			printf("\nStart SD:\t\t\t%f", start_sd);
			iterative_sigma_clip(this_spat_values, spat_nelements, clip_sigma, retain_indexes, start_mean, start_sd, &final_mean, &final_sd, &final_num_retained_indexes, FALSE);
			printf("\nFinal mean:\t\t\t%f", final_mean);
			printf("\nFinal SD:\t\t\t%f", final_sd);
			
			// 2.	Smooth array with median filter
			double this_spat_values_smoothed[spat_nelements];			
			memset(this_spat_values_smoothed, 0, sizeof(double)*spat_nelements);
			median_filter(this_spat_values, this_spat_values_smoothed, spat_nelements, median_filter_width_px);
			
			// 3.	Ascertain if this bin contains target flux
			int num_pixels_contain_target_flux = 0;
			for (jj=0; jj<spat_nelements-1; jj++) {
				if (this_spat_values_smoothed[jj] > final_mean + final_sd*min_SNR) {
					num_pixels_contain_target_flux++;
				}
			}
			printf("\nNum pixels (target):\t\t%d", num_pixels_contain_target_flux);
			
			printf("\nDoes bin contain target flux?\t");
			if (num_pixels_contain_target_flux >= min_spatial_width_px) {
				printf("Yes\n");
			} else {
				printf("No\n");
				peaks[ii] = -1;
				continue;
			}
			
			// 3.	Take derivatives
			double this_spat_values_der[spat_nelements-1];
			memset(this_spat_values_der, 0, sizeof(double)*spat_nelements-1);
			for (jj=1; jj<spat_nelements; jj++) {
				this_spat_values_der[jj-1] = this_frame_values_binned[jj][ii] - this_frame_values_binned[jj-1][ii];
			}
			
			// 4.	Smooth derivatives
			double this_spat_values_der_smoothed[spat_nelements-1];	
			memcpy(this_spat_values_der_smoothed, this_spat_values_der, sizeof(double)*spat_nelements-1);
			median_filter(this_spat_values_der, this_spat_values_der_smoothed, spat_nelements-1, median_filter_width_px);				
			
			// 5.	Pick most positive gradient from window, retaining proper index   
                        double this_spat_values_der_smoothed_windowed[spat_nelements-1]; 
                        memcpy(this_spat_values_der_smoothed_windowed, this_spat_values_der_smoothed, sizeof(double)*spat_nelements-1);
                        for (jj=0; jj<spat_nelements; jj++) {
                            if (jj >= finding_window_lo_px && jj <= finding_window_hi_px) {
                                this_spat_values_der_smoothed_windowed[jj] = this_spat_values_der_smoothed_windowed[jj];
                            } else {
                                this_spat_values_der_smoothed_windowed[jj] = -1;
                            }  
                        }
			size_t this_pk_idx = gsl_stats_max_index(this_spat_values_der_smoothed_windowed, 1, spat_nelements-1);
			printf("Start peak index:\t\t%d\n", this_pk_idx);	
			
			// 6.	Using this index, walk through centering window [max_centering_num_px] and find derivative turnover point
			printf("Found turnover:\t\t\t");			
			bool found_turnover = FALSE;
			for (jj=this_pk_idx; jj<this_pk_idx+max_centering_num_px; jj++) {
				if (this_spat_values_der_smoothed[jj] < 0.) {
					this_pk_idx = jj;
					found_turnover = TRUE;
					break;
				}
			}
			
			if (found_turnover) {
				printf("Yes\n");
				printf("End peak index:\t\t\t%d\n", this_pk_idx);	
                                num_bins_used++;
			} else {
				printf("No\n");
				peaks[ii] = -1;
				continue;
			}

			// 7.	Get parabolic centroid using the full centroid window [centroid_half_window_size_px]
			double this_pk_window_idxs[1 + (2*centroid_half_window_size_px)];
			double this_pk_window_vals[1 + (2*centroid_half_window_size_px)];
			
			memset(this_pk_window_idxs, 0, sizeof(double)*(1 + (2*centroid_half_window_size_px)));
			memset(this_pk_window_vals, 0, sizeof(double)*(1 + (2*centroid_half_window_size_px)));
			idx = 0;
			for (jj=this_pk_idx-centroid_half_window_size_px; jj<=this_pk_idx+centroid_half_window_size_px; jj++) {
				this_pk_window_idxs[idx] = jj;
				this_pk_window_vals[idx] = this_frame_values_binned[jj][ii];
				idx++;
			}	
			
			int order = 2;
			double coeffs [order+1];  
			memset(coeffs, 0, sizeof(double)*order+1);
			double chi_squared;
			if (calc_least_sq_fit(2, 1 + (2*centroid_half_window_size_px), this_pk_window_idxs, this_pk_window_vals, coeffs, &chi_squared)) {

				write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -7, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

				free(target_f);
				if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

				return 1; 		

			}
			
			double fitted_peak_idx = -coeffs[1]/(2*coeffs[2]);
			printf("Fitted peak index:\t\t%f\n", fitted_peak_idx);

			peaks[ii] = fitted_peak_idx;	
			
			/*for (jj=0; jj<spat_nelements-1; jj++) {
				printf("%d\t%f\t%f\n", jj, this_spat_values_der_smoothed[jj], this_spat_values_der[jj]);
			}*/ // DEBUG	
			
		}	
		
		printf("\nNum bins used:\t%d\n", num_bins_used);		
		
		if (num_bins_used < min_used_bins) {
			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -8, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(target_f);
			if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

			return 1;
		}		
		
		// ***********************************************************************
		// Create [FRFIND_OUTPUTF_PEAKS_FILE] output file and print a few 
		// parameters

		FILE *outputfile;
		outputfile = fopen(LOFIND_OUTPUTF_PEAKS_FILE, FILE_WRITE_ACCESS);

		if (!outputfile) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -9, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

			free(target_f);
			if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

			return 1;


		}

		char timestr [80];
		memset(timestr, '\0', sizeof(char)*80);

		find_time(timestr);

		fprintf(outputfile, "#### %s ####\n\n", LOFIND_OUTPUTF_PEAKS_FILE);
		fprintf(outputfile, "# Lists the coordinates of the peaks found using the lofind routine.\n\n");
		fprintf(outputfile, "# Run filename:\t%s\n", target_f);
		fprintf(outputfile, "# Run datetime:\t%s\n\n", timestr);
		
		for (ii=0; ii<disp_nelements_binned; ii++) {
			if (peaks[ii] == -1)
				continue;
			
			fprintf(outputfile, "%f\t%f\n", cut_x[0]+(ii*bin_size_px) + (double)bin_size_px/2., peaks[ii]);
		}
		
		fprintf(outputfile, "%d", EOF);		
		
		// ***********************************************************************
		// Clean up heap memory

		free(target_f);

		// ***********************************************************************
		// Close [FRFIND_OUTPUTF_PEAKS_FILE] output file and target file (ARG 1)
		
		if (fclose(outputfile)) {

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -10, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

			if(fits_close_file(target_f_ptr, &target_f_status)) fits_report_error (stdout, target_f_status); 

			return 1; 

		}

		if(fits_close_file(target_f_ptr, &target_f_status)) { 

			write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", -11, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);
			fits_report_error (stdout, target_f_status); 

			return 1; 

	    	}		
		
		// ***********************************************************************
		// Write success to [ERROR_CODES_FILE]

		write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATFI", RETURN_FLAG, "Status flag for L2 lofind routine", ERROR_CODES_FILE_WRITE_ACCESS);

		return 0;

	}

}


