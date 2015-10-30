#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "lotus_error_handling.h"
#include "lotus_functions.h"
#include "lotus_config.h"
#include "lotus_red_arcfit.h"

#include <gsl/gsl_statistics_int.h>

int main (int argc, char *argv []) {

        if(populate_env_variable(REF_ERROR_CODES_FILE, "L2_ERROR_CODES_FILE")) {

                printf("\nUnable to populate [REF_ERROR_CODES_FILE] variable with corresponding environment variable. Routine will proceed without error handling\n");

        }

        if (argc != 11) {

                if(populate_env_variable(LOA_BLURB_FILE, "L2_LOA_BLURB_FILE")) {

                        RETURN_FLAG = 1;

                } else {

                        print_file(LOA_BLURB_FILE);

                }

                write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -1, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                return 1;

        } else {
                        
                // ***********************************************************************
                // Redefine routine input parameters
        
                char *ext_arc_f                 = strdup(argv[1]);
                int min_dist                    = strtol(argv[2], NULL, 0);
                int half_aperture_num_pix       = strtol(argv[3], NULL, 0);
                int derivative_tol              = strtol(argv[4], NULL, 0);
                int derivative_tol_ref_px       = strtol(argv[5], NULL, 0);
                char *arc_line_list_filename    = strdup(argv[6]);
                int max_pix_diff                = strtol(argv[7], NULL, 0);
                int min_matched_lines           = strtol(argv[8], NULL, 0);
                double max_av_wavelength_diff   = strtol(argv[9], NULL, 0);
                int fit_order                   = strtol(argv[10], NULL, 0);

                // ***********************************************************************
                // Open ext arc file (ARG 1), get parameters and perform any data format 
                // checks

                fitsfile *ext_arc_f_ptr;

                int ext_arc_f_maxdim = 2, ext_arc_f_status = 0, ext_arc_f_bitpix, ext_arc_f_naxis;
                long ext_arc_f_naxes [2] = {1,1};

                if(!fits_open_file(&ext_arc_f_ptr, ext_arc_f, READONLY, &ext_arc_f_status)) {

                        if(!populate_img_parameters(ext_arc_f, ext_arc_f_ptr, ext_arc_f_maxdim, &ext_arc_f_bitpix, &ext_arc_f_naxis, ext_arc_f_naxes, &ext_arc_f_status, "ARC FRAME")) {

                                if (ext_arc_f_naxis != 2) {     // any data format checks here

                                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -2, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                                        free(ext_arc_f);
                                        free(arc_line_list_filename);

                                        if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                                        return 1;
        
                                }

                        } else { 

                                write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -3, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
                                fits_report_error(stdout, ext_arc_f_status); 

                                free(ext_arc_f);
                                free(arc_line_list_filename);

                                if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                                return 1; 

                        }

                } else { 

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -4, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
                        fits_report_error(stdout, ext_arc_f_status); 

                        free(ext_arc_f);
                        free(arc_line_list_filename);

                        return 1; 

                }

                // ***********************************************************************
                // Set the range limits using arc fits file

                int cut_x [2] = {1, ext_arc_f_naxes[0]};
                int cut_y [2] = {1, ext_arc_f_naxes[1]};

                // ***********************************************************************
                // Set parameters used when reading data from the arc file (ARG 1)

                long fpixel [2] = {cut_x[0], cut_y[0]};
                long nxelements = (cut_x[1] - cut_x[0]) + 1;
                long nyelements = (cut_y[1] - cut_y[0]) + 1;

                // ***********************************************************************
                // Create arrays to store pixel values from the arc file (ARG 1)

                double ext_arc_f_pixels [nxelements];

                // ***********************************************************************
                // Get arc fits file (ARG 1) values and store in 2D array

                int ii, jj;

                double ext_arc_frame_values [nyelements][nxelements];
                memset(ext_arc_frame_values, 0, sizeof(double)*nxelements*nyelements);

                for (fpixel[1] = cut_y[0]; fpixel[1] <= cut_y[1]; fpixel[1]++) {

                        memset(ext_arc_f_pixels, 0, sizeof(double)*nxelements);

                        if(!fits_read_pix(ext_arc_f_ptr, IMG_READ_ACCURACY, fpixel, nxelements, NULL, ext_arc_f_pixels, NULL, &ext_arc_f_status)) {

                                for (ii=0; ii<nxelements; ii++) {

                                        ext_arc_frame_values[fpixel[1]-1][ii] = ext_arc_f_pixels[ii];

                                }

                        } else { 

                                write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -5, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
                                fits_report_error(stdout, ext_arc_f_status); 

                                free(ext_arc_f);
                                free(arc_line_list_filename);

                                if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                                return 1; 

                        }

                }      
                
                // ***********************************************************************
                // Collapse array along spatial axis
                
                double col_ext_arc_frame_values [nxelements];
                memset(col_ext_arc_frame_values, 0, sizeof(double)*nxelements);
                for (ii=0; ii<nxelements; ii++) {
                        for (jj=0; jj<nyelements; jj++) {
                            col_ext_arc_frame_values[ii] += ext_arc_frame_values[jj][ii];
                        }
                }                            
                
                // ARC LINE MATCHING
                // ***********************************************************************
                // 1.   Find peaks in arc frame

                int num_peaks = 0;      
                int peaks[nxelements];                                            

                find_peaks(nxelements, col_ext_arc_frame_values, peaks, &num_peaks, min_dist, half_aperture_num_pix, derivative_tol, derivative_tol_ref_px, INDEXING_CORRECTION);
                
                printf("\nArc line matching");
                printf("\n-----------------\n\n");
                printf("Candidate lines found:\t\t\t\t\t%d\n", num_peaks);

                // 2.   Find parabolic centroids of contiguous peaks in arc 
                //      frame

                double peak_centroids [num_peaks];
                memset(peak_centroids, 0, sizeof(double)*num_peaks); 

                find_centroid_parabolic(col_ext_arc_frame_values, peaks, num_peaks, peak_centroids, INDEXING_CORRECTION);

                // 3.   Open reference arc line list file (ARG 6) and count
                //      number of reference lines
        
                FILE *arc_line_list_f;
        
                if (!check_file_exists(arc_line_list_filename)) { 

                        arc_line_list_f = fopen(arc_line_list_filename , "r");

                } else {

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -7, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                        free(ext_arc_f);
                        free(arc_line_list_filename);

                        if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status);  

                        return 1;

                }       
        
                char input_string[300];

                int arc_line_list_num_ref_lines = 0;
        
                while(!feof(arc_line_list_f)) {

                        memset(input_string, '\0', sizeof(char)*300);

                        fgets(input_string, 300, arc_line_list_f);      

                        if (strtol(&input_string[0], NULL, 0) > 0) {            // check the line begins with a positive number (usable)
        
                                arc_line_list_num_ref_lines++;                          

                        }
        
                }

                printf("Number of reference lines:\t\t\t\t%d\n", arc_line_list_num_ref_lines);

                // 4.   Read arc line list file (ARG 6) for reference
                //      wavelengths and their corresponding pixel positions

                rewind(arc_line_list_f);

                double arc_peak_centroids [arc_line_list_num_ref_lines];
                memset(arc_peak_centroids, 0, sizeof(double)*arc_line_list_num_ref_lines);

                double arc_peak_wavelengths [arc_line_list_num_ref_lines];
                memset(arc_peak_wavelengths, 0, sizeof(double)*arc_line_list_num_ref_lines);

                double pixel_channel, wavelength;

                int arc_line_list_line_index = 0;

                while(!feof(arc_line_list_f)) {

                        memset(input_string, '\0', sizeof(char)*300);

                        fgets(input_string, 300, arc_line_list_f);      

                        if (strtol(&input_string[0], NULL, 0) > 0) {            // check the line begins with a positive number (usable)
        
                                sscanf(input_string, "%lf\t%lf\t", &pixel_channel, &wavelength); 

                                arc_peak_centroids[arc_line_list_line_index] = pixel_channel;
                                arc_peak_wavelengths[arc_line_list_line_index] = wavelength;    

                                arc_line_list_line_index++;     

                        }
        
                }

                // for (ii=0; ii<arc_line_list_line_index; ii++) printf("\n%f\t%f", arc_peak_centroids[ii], arc_peak_wavelengths[ii]);  // DEBUG

                // 5.   Match identified peaks to lines in reference file

                double this_diff, least_diff, least_diff_index;
                bool matched_line;
                int matched_line_count = 0;

                int matched_line_indexes [num_peaks];

                for (ii=0; ii<num_peaks; ii++) {

                        matched_line_indexes[ii] = -1;                  // -1 indicates no match

                }

                double matched_line_diffs [num_peaks];
                memset(matched_line_diffs, 0, sizeof(double)*num_peaks);  

                int duplicate_index;
                
                for (ii=0; ii<num_peaks; ii++) {                                                        // for each identified candidate peak

                        least_diff = 0;
                        matched_line = FALSE;

                        for (jj=0; jj<arc_line_list_num_ref_lines; jj++) {                              // then take each line in the reference arc line list file

                                this_diff = fabs(peak_centroids[ii] - arc_peak_centroids[jj]);          // calculate the difference between the identified and reference arc line

                                if (this_diff <= max_pix_diff) {                                        // comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function
                        
                                        if (matched_line == FALSE || this_diff < least_diff) {          // is this the first line found or a line with a smaller difference between ref arc line and matched?

                                                matched_line = TRUE;

                                                least_diff = this_diff;
                                                least_diff_index = jj;

                                        }
                                
                                }

                        }

                        if (matched_line == TRUE) {

                                if (lsearch_int(matched_line_indexes, least_diff_index, num_peaks) != -1) {                     // has this line already been allocated to another peak?

                                        duplicate_index = lsearch_int(matched_line_indexes, least_diff_index, num_peaks);

                                        if(least_diff < matched_line_diffs[duplicate_index]) {                                  // does it have a smaller difference? - comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

                                                matched_line_indexes[duplicate_index] = -1;                                     // reset the duplicate values, no increment required if duplicate
                                                matched_line_diffs[duplicate_index] = 0;

                                                matched_line_indexes[ii] = least_diff_index;
                                                matched_line_diffs[ii] = least_diff;
                                                
                                        }

                                } else {

                                        matched_line_indexes[ii] = least_diff_index;
                                        matched_line_diffs[ii] = least_diff;
                                        matched_line_count++;

                                }

                        } else {

                                matched_line_indexes[ii] = -1;                                          // not matched

                        }

                }

                // 6.   Did we match enough lines?

                if (matched_line_count < min_matched_lines) {

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -8, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                        free(ext_arc_f);
                        free(arc_line_list_filename);

                        fclose(arc_line_list_f);

                        if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 
                        return 1;

                }

                printf("Number of matched lines:\t\t\t\t%d\n", matched_line_count);

                printf("\nIndex\tList wavelength\tList centroid\tAv. channel\tChannel difference\n");           // Channel difference is the difference between the reference arc list pixel location and the identified pixel location
                printf("\t(Å)\t\t(px)\t\t(px)\t\t(px)\n\n");    

                for (ii=0; ii<num_peaks; ii++) {

                        if (matched_line_indexes[ii] == -1) {   // this line was unmatched

                                // printf("%d\t%s\t\t%s\t\t%.2f\t\t%s\n", ii, "", "", average_pixel_channels[ii], "");  // DEBUG
                                continue;

                        } else {

                                printf("%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n", ii, arc_peak_wavelengths[matched_line_indexes[ii]], arc_peak_centroids[matched_line_indexes[ii]], peak_centroids[ii], matched_line_diffs[ii]);

                        }

                }


                // 7.   And do they sample the distribution well?

                // for (ii=0; ii<arc_line_list_line_index; ii++) printf("\n%f\t%f\t%d", arc_peak_wavelengths[ii], arc_peak_wavelengths[matched_line_indexes[ii]], matched_line_indexes[ii]);    // DEBUG

                double list_total_dist = 0.0;
                double list_av_dist;

                for (ii=0; ii<arc_line_list_num_ref_lines-1; ii++) {

                        for (jj=ii+1; jj<arc_line_list_num_ref_lines; jj++) {
                                
                                list_total_dist += arc_peak_wavelengths[jj] - arc_peak_wavelengths[ii];

                                // printf("\n%f", arc_peak_wavelengths[jj] - arc_peak_wavelengths[ii]);         // DEBUG

                        }

                }

                list_av_dist = list_total_dist / powf(arc_line_list_line_index, 2);

                double matched_line_wavelengths [matched_line_count];
                memset(matched_line_wavelengths, 0, sizeof(double)*matched_line_count);
                
                int this_matched_line_index = 0;

                for (ii=0; ii<num_peaks; ii++) {
                
                        if (matched_line_indexes[ii] == -1) {   // this line was unmatched

                                continue;

                        } else {

                                matched_line_wavelengths[this_matched_line_index] = arc_peak_wavelengths[matched_line_indexes[ii]];
                                this_matched_line_index++;

                        }

                }

                double sample_total_dist = 0.0;
                double sample_av_dist;

                for (ii=0; ii<matched_line_count-1; ii++) {

                        for (jj=ii+1; jj<matched_line_count; jj++) {
                                
                                sample_total_dist += matched_line_wavelengths[jj] - matched_line_wavelengths[ii];

                                // printf("\n%f", matched_line_wavelengths[jj] - matched_line_wavelengths[ii]);         // DEBUG

                        }

                }

                sample_av_dist = sample_total_dist / powf(matched_line_count, 2);

                double sample_list_diff = abs(list_av_dist-sample_av_dist);

                printf("\nAverage distance between lines in reference arc list:\t%.1f", list_av_dist);
                printf("\nAverage distance between lines in matched line list:\t%.1f\n", sample_av_dist);

                if (sample_list_diff > max_av_wavelength_diff) {        // comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

                        RETURN_FLAG = 2;

                }

                // FIND THE DISPERSION SOLUTION AND WRITE TO 
                // [LOARCFIT_OUTPUTF_WAVFITS_FILE] OUTPUT FILE
                // ***********************************************************************

                // 1.   Perform a few checks to ensure the input fitting parameters 
                //      are sensible

                if ((fit_order < LOARCFIT_VAR_POLYORDER_LO) || (fit_order > LOARCFIT_VAR_POLYORDER_HI)) {       // Check [fit_order] is within config limits

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -9, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                        free(ext_arc_f);
                        free(arc_line_list_filename);

                        fclose(arc_line_list_f);

                        if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                        return 1; 

                }

                // 2.   Create [LOARCFIT_OUTPUTF_WAVFITS_FILE] output file and print a few 
                //      parameters

                FILE *outputfile;
                outputfile = fopen(LOARCFIT_OUTPUTF_WAVFITS_FILE, FILE_WRITE_ACCESS);

                if (!outputfile) { 

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -10, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                        free(ext_arc_f);
                        free(arc_line_list_filename);

                        fclose(arc_line_list_f);

                        if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                        return 1;

                }

                char timestr [80];
                memset(timestr, '\0', sizeof(char)*80);

                find_time(timestr);

                fprintf(outputfile, "#### %s ####\n\n", LOARCFIT_OUTPUTF_WAVFITS_FILE);
                fprintf(outputfile, "# List of pixel channel to wavelength dispersion solution coefficients and corresponding chi-squared found using the loarcfit program.\n\n");
                fprintf(outputfile, "# Run Datetime:\t\t%s\n\n", timestr);
                fprintf(outputfile, "# Arc Filename:\t\t%s\n", ext_arc_f);
                fprintf(outputfile, "# Polynomial Order:\t%d\n\n", fit_order);

                // 3.   Populate arrays and perform fits

                double coord_x [matched_line_count];
                double coord_y [matched_line_count];

                double coeffs [fit_order+1];

                double this_chi_squared;

                int line_count;

                memset(coord_x, 0, sizeof(double)*matched_line_count);
                memset(coord_y, 0, sizeof(double)*matched_line_count);
                memset(coeffs, 0, sizeof(double)*fit_order+1);

                line_count = 0;

                for (ii=0; ii<num_peaks; ii++) {
                
                    if (matched_line_indexes[ii] == -1) {                           // this line was unmatched

                        continue;

                    } else {

                        coord_x[line_count] = peak_centroids[ii];
                        coord_y[line_count] = arc_peak_wavelengths[matched_line_indexes[ii]];
                        line_count++;

                    }

                }

                if (calc_least_sq_fit(fit_order, matched_line_count, coord_x, coord_y, coeffs, &this_chi_squared)) {    

                    write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -11, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                    free(ext_arc_f);
                    free(arc_line_list_filename);

                    fclose(arc_line_list_f);
                    fclose(outputfile);

                    if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                    return 1; 

                }

                // Output solutions to [LOARCFIT_OUTPUTF_WAVFITS_FILE] file 

                fprintf(outputfile, "%d\t", 1);

                for (ii=0; ii<=fit_order; ii++) {
                  
                    fprintf(outputfile, LOARCFIT_VAR_ACCURACY_COEFFS, coeffs[ii]);
                    fprintf(outputfile, "\t");

                }

                fprintf(outputfile, LOARCFIT_VAR_ACCURACY_CHISQ, this_chi_squared);
                fprintf(outputfile, "\n");

                // printf("%d\t%f\n", ii, this_chi_squared);    // DEBUG

                fprintf(outputfile, "%d", EOF);

                printf("\nFitting results");
                printf("\n---------------\n");
                printf("χ2:\t\t%.2f\n", this_chi_squared);

                // 4.   Perform a few checks to ensure the chi squareds are sensible 

                if ((this_chi_squared < LOARCFIT_VAR_CHISQUARED_MIN) || (this_chi_squared > LOARCFIT_VAR_CHISQUARED_MAX)) {

                        RETURN_FLAG = 3;

                }

                // ***********************************************************************
                // Clean up heap memory

                free(ext_arc_f);
                free(arc_line_list_filename);
                
                // ***********************************************************************
                // Close input files (ARGS 1,2 and 3), arc list file (ARG 6) and 
                // [LOARCFIT_OUTPUTF_WAVFITS_FILE] file

                if (fclose(arc_line_list_f)) {

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -12, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                        fclose(outputfile);

                        if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                        return 1; 

                }

                if (fclose(outputfile)) {

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -13, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                        if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) fits_report_error (stdout, ext_arc_f_status); 

                        return 1; 

                }

                if(fits_close_file(ext_arc_f_ptr, &ext_arc_f_status)) { 

                        write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", -14, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);
                        fits_report_error (stdout, ext_arc_f_status); 

                        return 1; 

                }
                
                // ***********************************************************************
                // Write success to [ERROR_CODES_FILE]

                write_key_to_file(ERROR_CODES_FILE, REF_ERROR_CODES_FILE, "L2STATAR", RETURN_FLAG, "Status flag for L2 loarcfit routine", ERROR_CODES_FILE_WRITE_ACCESS);

                return 0;

        }

}

