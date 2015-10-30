#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/stat.h>
#include <ctype.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl_randist.h>

/************************************************************************

 Function:		find_edges
 Purpose:		finds the edges in a dataset
 Additional Notes:	


 [INDEXING CORRECTION]:

 This variable defines how we offset the arrays.

 If the correction is set to TRUE, [peaks] will store the x value of the
 identified peak and not the ii value (x-1).

************************************************************************/
/*
int find_edges(int nxelements, double row_values [], int peaks [], int * num_peaks, int min_dist, int half_aperture_num_pix, int derivative_tol, int derivative_tol_ref_px, bool INDEXING_CORRECTION) {

	int ii, jj;
	int peak_count = 0;
	
	for (ii = half_aperture_num_pix; ii < nxelements - half_aperture_num_pix; ii++) {  		// this is the pixel we're considering (ii)
	
		for (jj = ii - half_aperture_num_pix; jj <= ii + half_aperture_num_pix; jj++) {		// this is the aperture we're considering

			if (jj == ii) { 

				continue; // this is the pixel in question, skip

			} else if (row_values[jj] > row_values[ii]) { 

				break;	  // this pixel is not a peak within the specified aperture

			} else if ((jj == ii + half_aperture_num_pix) && ((row_values[ii] - row_values[ii-derivative_tol_ref_px]) > derivative_tol && (row_values[ii] - row_values[ii+derivative_tol_ref_px]) > derivative_tol)) {	// reached last pixel in aperture => pixel ii is a peak within it. Check that it is sufficiently brighter than the derivative tolerance value / pixel. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

				if (peak_count == 0) {	// if this is the first peak, we can just add it without checking the minimum distance parameter (no other peaks exist)

					peaks[peak_count] = ii + INDEXING_CORRECTION;	// see INDEXING_CORRECTION
					peak_count++;

				} else if (ii - peaks[peak_count-1] > min_dist) { // else we need to check the minimum distance criteria is satisfied
	
					peaks[peak_count] = ii + INDEXING_CORRECTION;	// see INDEXING_CORRECTION
					peak_count++;
	
				}
	
			}
	
		}
	
	}
		
	*num_peaks = peak_count;

	return 0;

}
*/

/************************************************************************

 Function:		calc_least_sq_fit
 Purpose:		finds the coefficients of the [order] fit, for
			[equations], whose values are defined in [array_x]
			and [array_y]. stores in [coeff] where coeff[0] 
			is the 0th order of fit with a corresponding 
			[chi_squared]
 Additional Notes:


 Performs least-squares fitting to a general model y = Xc using the GSL
 library.

 			Array formats
		 	-------------

 where m is the order and n is the number of equations.

 y = [	y_0, y_1, y_2, ..., y_n;	]

 X = [	1, (x_0)^1, (x_0)^2, ..., (x_0)^m;
	1, (x_1)^1, (x_1)^2, ..., (x_1)^m;
	1, (x_2)^1, (x_2)^2, ..., (x_2)^m;
	1, (x_3)^1, (x_3)^2, ..., (x_3)^m;
	.	.	.	.	.
	.	.	.	.	.
	.	.	.	.	.
	1, (x_n)^1, (x_n)^2, ..., (x_n)^m; ]

 c = [	c_0, c_1, c_2, ..., c_n+1;	]

 n.b. matrices are are denoted in (rows x columns) format.

************************************************************************/

int calc_least_sq_fit(int order, int equations, double array_x [], double array_y [], double coeffs [], double *chi_squared) {

	if (order >= equations)	{	// then we don't have enough equations to solve for all parameters

		return 1;

	}

	int ii, jj;

	gsl_matrix *X, *cov;
	gsl_vector *y, *c;

	double chisq;

	X = gsl_matrix_alloc(equations, order+1);	// matrix of predictor variables
	y = gsl_vector_alloc(equations);		// vector of observed values

	c = gsl_vector_alloc(order+1);			// coefficient vector
	cov = gsl_matrix_alloc(order+1, order+1);	// covariance matrix

	for (ii=0; ii<equations; ii++) {					// for each equation
			
		for (jj=0; jj<order+1; jj++) {
			   
			gsl_matrix_set(X, ii, jj, powf(array_x[ii], jj)); 	// populate X matrix with the value of x^m for each order m.

		}
			   
		gsl_vector_set(y, ii, array_y[ii]);				// and populate y matrix

	}

	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(equations, order+1);	// define a workspace
	gsl_multifit_linear(X, y, c, cov, &chisq, work);					// perform the fit
	gsl_multifit_linear_free(work);								// free the workspace

	#define C(i) (gsl_vector_get(c,(i)))

	for (ii=0; ii<order+1; ii++) {

		coeffs[ii] = C(ii);	// populate coeff array
		
		// printf("%e\t", coeffs[ii]);	if (ii==order) { printf("\n"); } 	// DEBUG

	}

	*chi_squared = chisq;		// write chisq
	     
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	     
	return 0;

}

/************************************************************************

 Function:		calculate_cross_correlation_offset
 Purpose:		calculates the best offset for two given
			datasets, x and y
 Additional Notes:	

 The cross correlation coefficient, r, is determined by the equation:

 r =

		sum_i[(x(i) - mean_x) * y(i-offset) - mean_y)]
 ------------------------------------------------------------------------
 sqrt_(sum_i((x(i) - mean_x)^2)) * sqrt_(sum_i((y(i-offset) - mean_y)^2)) 

 a +ve offset indicates the series lags the reference series
 a -ve offset indicates the series leads the references series		

************************************************************************/

int calculate_cross_correlation_offset(double x [], double y [], int n, int max_offset, int *this_best_offset, double *this_best_r) {

	double mean_x = gsl_stats_mean(x, 1, n);
	double mean_y = gsl_stats_mean(y, 1, n);

	int offset, best_offset;
	double r, best_r;

	for (offset=-max_offset; offset<=max_offset; offset++) {
	
		double sx, sy, numerator = 0.0, den_x_term = 0.0, den_y_term = 0.0, denominator;

		int ii, jj;
	
		for (ii=0; ii<n; ii++) {
		
			jj = ii + offset;

			if (jj < 0 || jj >= n) {	// element is outside of array boundaries

				continue;

			} else {

				sx = (x[ii] - mean_x);
				sy = (y[jj] - mean_y);

				numerator += sx * sy;
				den_x_term += powf(sx, 2);
				den_y_term += powf(sy, 2);
			
			}

			denominator = sqrt(den_x_term) * sqrt(den_y_term);
	
		}

		r = numerator / denominator;

		if ((offset == -max_offset) || (r > best_r)) {	// Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

			best_r = r;
			best_offset = offset;

		} 

		// printf("%d\t%f\t%f\t%f\n", offset, numerator/denominator, numerator, denominator);	// DEBUG

	}

	*this_best_offset = best_offset;
	*this_best_r = best_r;

	// printf("%d\t%f\n", best_offset, best_r);		// DEBUG

	return 0;

}

/************************************************************************

 Function:		check_file_exists
 Purpose:		checks a file to see if it exists
 Additional Notes:	None

************************************************************************/

int check_file_exists(char filename []) {

	FILE *fp = fopen(filename,"r");
	
	if (fp) {

		fclose(fp);
		return 0;
	
	} else {
	
		return 1;
	
	}

}

/************************************************************************

 Function:		check_key_to_omit
 Purpose:		checks the [FITS_KEYS_TO_OMIT] file for the 
			existence of [card] in hdu number [hdunum]
			and assigns [found_key] to TRUE if found
 Additional Notes:	None

************************************************************************/

int check_key_to_omit(char * FITS_KEYS_TO_OMIT, char * card, char * operation, int * found_key) {

	FILE *FITS_KEYS_TO_OMIT_FILE_ptr;
	FITS_KEYS_TO_OMIT_FILE_ptr = fopen(FITS_KEYS_TO_OMIT, "r");

	char input_string [300];

	char file_fits_op [100];
	char file_fits_key [81];

	while(!feof(FITS_KEYS_TO_OMIT_FILE_ptr)) {

		memset(input_string, '\0', sizeof(char)*300);

		memset(file_fits_op, '\0', sizeof(char)*100);
		memset(file_fits_key, '\0', sizeof(char)*81);

		fgets (input_string, 300, FITS_KEYS_TO_OMIT_FILE_ptr);
	
		sscanf(input_string, "%[^\t]\t%[^\n]", file_fits_op, file_fits_key);	// last argument in sscanf means read characters as one string up to newline char

		if (isalpha(input_string[0])) {	// check the line begins with a letter

			// printf("%s\t%s\n", file_fits_op, file_fits_key);										// DEBUG
			// printf("%s\t%.5s\n", operation, card);											// DEBUG
			// printf("%d\t%d\n", (strncmp(card,file_fits_key,strlen(file_fits_key)) == 0), (strcmp(file_fits_op, operation) == 0));	// DEBUG

			if ((strncmp(card,file_fits_key,strlen(file_fits_key)) == 0) && (strcmp(file_fits_op, operation) == 0)) {

				// printf("%d\n", TRUE);	// DEBUG
				*found_key = TRUE;
				break;
		
			}

		}

	} 

	fclose(FITS_KEYS_TO_OMIT_FILE_ptr);

	return 0;
	
}

/************************************************************************

 Function:		find_centroid_parabolic	
 Purpose:		finds the parabolic centroid of each peak in a
			dataset given the position of the peaks
 Additional Notes:	


 [INDEXING CORRECTION]:

 This variable defines how if we offset the arrays.

 If the correction is set to TRUE, [peaks_centroids] will store the x 
 value of the identified peak and not the ii value (x-1).

************************************************************************/

int find_centroid_parabolic(double row_values [], int peaks [], int num_peaks, double peak_centroids [], bool INDEXING_CORRECTION) {
	     
	int ii, jj;

	int order = 2;		// this is a quadratic procedure
	int equations = 3;	// using 3 pixels (i.e. fully defined)

	double array_x [equations];
	memset(array_x, 0, sizeof(double)*equations);

	double array_y [equations];
	memset(array_y, 0, sizeof(double)*equations);

	double coeffs [order+1];  
	memset(coeffs, 0, sizeof(double)*order+1);
 
	double chi_squared;

	for (ii=0; ii<num_peaks; ii++) {		// for each peak

		for (jj=0; jj<=order; jj++) {		// populate array_x and array_y with the appropriate variables
			
			array_x[jj] = peaks[ii] + (jj-1);			
			array_y[jj] = row_values[peaks[ii]+(jj-1)-INDEXING_CORRECTION];	// see INDEXING_CORRECTION

		}

		if (calc_least_sq_fit(order, equations, array_x, array_y, coeffs, &chi_squared)) {

			return 1;

		} else {

			peak_centroids[ii] = -coeffs[1]/(2*coeffs[2]);	// store centroid (derivative of y = ax^2 + bx + c is -b/2a)

		}

	}

	return 0;

}

/************************************************************************

 Function:		find_peaks
 Purpose:		finds the peaks in a dataset
 Additional Notes:	


 [INDEXING CORRECTION]:

 This variable defines how if we offset the arrays.

 If the correction is set to TRUE, [peaks] will store the x value of the
 identified peak and not the ii value (x-1).

************************************************************************/

int find_peaks(int nxelements, double row_values [], int peaks [], int * num_peaks, int min_dist, int half_aperture_num_pix, int derivative_tol, int derivative_tol_ref_px, bool INDEXING_CORRECTION) {

	int ii, jj;
	int peak_count = 0;
	
	for (ii = half_aperture_num_pix; ii < nxelements - half_aperture_num_pix; ii++) {  		// this is the pixel we're considering (ii)
	
		for (jj = ii - half_aperture_num_pix; jj <= ii + half_aperture_num_pix; jj++) {		// this is the aperture we're considering

			if (jj == ii) { 

				continue; // this is the pixel in question, skip

			} else if (row_values[jj] > row_values[ii]) { 

				break;	  // this pixel is not a peak within the specified aperture

			} else if ((jj == ii + half_aperture_num_pix) && ((row_values[ii] - row_values[ii-derivative_tol_ref_px]) > derivative_tol && (row_values[ii] - row_values[ii+derivative_tol_ref_px]) > derivative_tol)) {	// reached last pixel in aperture => pixel ii is a peak within it. Check that it is sufficiently brighter than the derivative tolerance value / pixel. Comparing doubles but accuracy isn't a necessity so don't need gsl_fcmp function

				if (peak_count == 0) {	// if this is the first peak, we can just add it without checking the minimum distance parameter (no other peaks exist)

					peaks[peak_count] = ii + INDEXING_CORRECTION;	// see INDEXING_CORRECTION
					peak_count++;

				} else if (ii - peaks[peak_count-1] > min_dist) { // else we need to check the minimum distance criteria is satisfied
	
					peaks[peak_count] = ii + INDEXING_CORRECTION;	// see INDEXING_CORRECTION
					peak_count++;
	
				}
	
			}
	
		}
	
	}
		
	*num_peaks = peak_count;

	return 0;

}

/************************************************************************

 Function:		find_peaks_contiguous
 Purpose:		finds contiguous peaks in a dataset
 Additional Notes:	None

************************************************************************/

int find_peaks_contiguous(int nxelements, int nyelements, double ** frame_values, int ** peaks, int * num_peaks, int min_dist, int half_aperture_num_pix, int derivative_tol, int derivative_tol_ref_px, int pix_tolerance, bool INDEXING_CORRECTION) {

	int ii, jj;

	double row_values [nxelements];

	int this_row_peaks [nxelements];					// array to store peak locations for a single iteration

	int this_row_num_peaks;							// variable to store number of peaks for a single iteration

	int all_rows_peaks [nyelements][nxelements];				// array to store peak locations for all iterations
	memset(all_rows_peaks, 0, sizeof(int)*nyelements*nxelements);	

	int all_rows_num_peaks [nyelements];					// array to store number of peaks for all iterations
	memset(all_rows_num_peaks, 0, sizeof(int)*nyelements*nxelements);

	for (jj=0; jj<nyelements; jj++) {

		memset(row_values, 0, sizeof(double)*nxelements);
		memset(this_row_peaks, 0, sizeof(int)*nxelements);	

		for (ii=0; ii<nxelements; ii++) {

			row_values[ii] = frame_values[jj][ii];

		}

		find_peaks(nxelements, row_values, this_row_peaks, &this_row_num_peaks, min_dist, half_aperture_num_pix, derivative_tol, derivative_tol_ref_px, INDEXING_CORRECTION);

		for (ii=0; ii<this_row_num_peaks; ii++) {

			all_rows_peaks[jj][ii] = this_row_peaks[ii];		// copy peaks locations across

		}

		all_rows_num_peaks[jj] = this_row_num_peaks;			// copy number of peaks across

	}

	int kk, index_to_insert = 0;

	bool insert;

	int considered_peaks [nyelements];	// array to hold locations of the considered peak

	double this_pix_difference = 0;
	double least_pix_difference;

	for (ii=0; ii<all_rows_num_peaks[0]; ii++) {									// for each peak on the first row

		memset(considered_peaks, 0, sizeof(int)*nyelements);

		considered_peaks[0] = all_rows_peaks[0][ii];

		for (jj=1; jj<nyelements; jj++) {									// cycle each futher row

			insert = FALSE;

			for (kk=0; kk<all_rows_num_peaks[jj]; kk++) {							// and each peak on the row

				this_pix_difference = abs(all_rows_peaks[jj][kk] - all_rows_peaks[0][ii]);		// calculate the pixel distance between this peak and the peak from the first row

				if (this_pix_difference <= pix_tolerance) {						// found a suitable candidate peak within the right tolerance

					if ((insert == FALSE) || (this_pix_difference < least_pix_difference)) {	// if it's the first cycle OR there's a closer peak

						insert = TRUE;
						considered_peaks[jj] = all_rows_peaks[jj][kk];
						least_pix_difference = this_pix_difference;				// this is always set on each iteration of jj as insert is set to FALSE

					}
					
				} 

			}

			if (insert == FALSE) {

				break;

			}

		}

		if (insert == TRUE) {	// considered_peaks can be inserted

			for (jj=0; jj<nyelements; jj++) {

				peaks[jj][index_to_insert] = considered_peaks[jj];

			}

			*num_peaks = *num_peaks + 1;
			index_to_insert++;

		}

	}

	return 0;

}

/************************************************************************

 Function:		find_time
 Purpose:		finds the time
 Additional Notes:	None

************************************************************************/

int find_time (char timestr []) {

	struct tm *time_now;
	time_t ltime;

	ltime = time(NULL);

	time_now = localtime(&ltime);

	strftime(timestr, 80, "%d-%m-%Y %H:%M:%S", time_now);

	return 0;

}

/************************************************************************

 Function:		flip_array_dbl
 Purpose:		Flips an array
 Additional Notes:	None

************************************************************************/

int flip_array_dbl(double array [], int size) {

	int ii;

	double flip_array [size];
	memset(flip_array, 0, sizeof(double)*size);

	for (ii=0; ii<size; ii++) {

		flip_array[size-1-ii] = array[ii];	// store reversed order of array elements to temporary array [flip_array]
	
	}

	memset(array, 0, sizeof(double)*size);

	for (ii=0; ii<size; ii++) {

		array[ii] = flip_array[ii];		// rewrite array [array] with the reversed array [flip_array]
	
	}

	return 0;

}

/************************************************************************

 Function:		interpolate
 Purpose:		interpolates a dataset [x] to find all values
			between [interpolation_start] and 
			[interpolation_end] with a spacing of [spacing]
 Additional Notes:	


 Interpolation types are specified by the GSL library.


************************************************************************/

int interpolate(char interpolation_type [], double x_wav [], double x_val [], int nxelements, double interpolation_start, double interpolation_end, double spacing, double x_val_out []) {
 
	gsl_spline *spline;

	if (strcmp(interpolation_type, "linear") == 0) {

		spline = gsl_spline_alloc(gsl_interp_linear, nxelements);

	} else if (strcmp(interpolation_type, "polynomial") == 0) {

		spline = gsl_spline_alloc(gsl_interp_polynomial, nxelements);

	} else if (strcmp(interpolation_type, "cspline") == 0) {

		spline = gsl_spline_alloc(gsl_interp_cspline, nxelements);

	} else if (strcmp(interpolation_type, "cspline_periodic") == 0) {

		spline = gsl_spline_alloc(gsl_interp_cspline_periodic, nxelements);

	} else if (strcmp(interpolation_type, "akima") == 0) {

		spline = gsl_spline_alloc(gsl_interp_akima, nxelements);

	} else if (strcmp(interpolation_type, "akima_periodic") == 0) {

		spline = gsl_spline_alloc(gsl_interp_akima_periodic, nxelements);

	} else {

		return 1;

	}



	gsl_interp_accel *acc = gsl_interp_accel_alloc();
		     
	gsl_spline_init(spline, x_wav, x_val, nxelements);

	int this_interpolation_index = 0;

	double xi;

	for (xi = interpolation_start; gsl_fcmp(xi, interpolation_end+spacing, 1e-5); xi += spacing) {	// checking to see if xi is equal to interpolation_end+spacing (i.e. no more iterations)	
		//printf("\n%f\t%f\t%g\t%d", xi, interpolation_end+spacing, x_val_out[this_interpolation_index], gsl_fcmp(xi, interpolation_end+spacing, 1e-5));	// DEBUG
		x_val_out[this_interpolation_index] = gsl_spline_eval(spline, xi, acc);

		this_interpolation_index++;
	}
		
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);

	return 0;

}

/************************************************************************

 Function:		iterative_sigma_clip
 Purpose:		Perform iterative sigma clipping on a dataset of 
			n values
 Additional Notes:	

 Before each iteration, a ceiling value is worked out which is found by:

  ceiling = [median] + [clip_sigma] * [sd]

 The elements in the dataset [values] are then cycled and any value that
 is greater than the ceiling value for this iteration are flagged FALSE
 in the [retain_indexes] array. This procedure is recursively done until
 the number of indexes retained for an iteration equals that of the last
 iteration.

 Along with the [retain_indexes] array, the final mean [final_mean],
 final standard deviation [final_sd] and final number of retained 
 indexes [final_num_retained_indexes] variables are populated.

************************************************************************/

int iterative_sigma_clip(double values [], int n, float clip_sigma, int retain_indexes [], double start_av, double start_sd, double * final_mean, double * final_sd, int * final_num_retained_indexes, int use_median) {
	
	int ii;
		
	int this_iteration_num_retained_indexes, last_iteration_num_retained_indexes = n;

	int this_iteration_values_index;

	int iteration_count = 0;

	bool continue_iterations = TRUE;

	for (ii=0; ii<n; ii++) {	// initialise to TRUE under assumption that all are to be retained

		retain_indexes[ii] = TRUE;

	}
	
	double mean, median, sd;
	double values_sorted [n];
	
	if (start_av == 0 || start_sd == 0) {
		memcpy(values_sorted, values, sizeof(double)*n);		
		gsl_sort(values_sorted, 1, n);

		mean = gsl_stats_mean(values_sorted, 1, n);
		median = gsl_stats_median_from_sorted_data(values_sorted, 1, n);
		sd = gsl_stats_sd(values_sorted, 1, n);		
	} else {
		mean	= start_av;
		median	= start_av;
		sd	= start_sd;
	}
	
	double ceiling;
	if (!use_median) {
		ceiling = mean + (double) clip_sigma*sd;
	} else {
		ceiling = median + (double) clip_sigma*sd;	
	}
	
	while (continue_iterations == TRUE) {

		this_iteration_num_retained_indexes = n;

		// ***********************************************************************
		// Determine which indexes to retain given current ceiling value

		for (ii=0; ii<n; ii++) {

			if (values[ii] > ceiling) {	// this index, ii, is not to be retained
	
				retain_indexes[ii] = FALSE;
				this_iteration_num_retained_indexes--;

			}

		}

		//	printf("%f\t%f\t%f\n", median, sd, ceiling);	// DEBUG

		// ***********************************************************************
		// Check to see if the number of retained indexes is the same as the
		// previous iteration. If this is the first iteration, check to see if the
		// number of retained indexes is equal to the original number of values.
		// If so, stop the iterations.

		if ((this_iteration_num_retained_indexes == last_iteration_num_retained_indexes) || this_iteration_num_retained_indexes <= 1) {

			break;

		}

		iteration_count++;	

		// ***********************************************************************
		// Create new arrays of retained values and work out new ceiling value

		double *this_iteration_values;
		this_iteration_values = (double *) malloc((this_iteration_num_retained_indexes)*sizeof(double));

		this_iteration_values_index = 0;

		for (ii=0; ii<n; ii++) {

			if (retain_indexes[ii] == TRUE) {

				this_iteration_values[this_iteration_values_index] = values[ii];		

				this_iteration_values_index++;

			}

		}

		printf("\nIteration:\t\t\t%d\n", iteration_count);
		printf("Ceiling:\t\t\t%.3e\n", ceiling);
		printf("Number of indexes retained:\t%d", this_iteration_num_retained_indexes);
		
		double this_iteration_values_sorted [this_iteration_num_retained_indexes];
		memcpy(this_iteration_values_sorted, this_iteration_values, sizeof(double)*this_iteration_num_retained_indexes);		
		gsl_sort(this_iteration_values_sorted, 1, this_iteration_num_retained_indexes);
	
		mean = gsl_stats_mean(this_iteration_values_sorted, 1, this_iteration_num_retained_indexes);
		median = gsl_stats_median_from_sorted_data(this_iteration_values_sorted, 1, this_iteration_num_retained_indexes);
		sd = gsl_stats_sd(this_iteration_values_sorted, 1, this_iteration_num_retained_indexes);

		double ceiling;
		if (!use_median) {
			ceiling = mean + (double) clip_sigma*sd;
		} else {
			ceiling = median + (double) clip_sigma*sd;	
		}

		last_iteration_num_retained_indexes = this_iteration_num_retained_indexes;
	
		free(this_iteration_values);
	}

	*final_mean = mean;
	*final_sd = sd;

	*final_num_retained_indexes = this_iteration_num_retained_indexes;

	return 0;

}

/************************************************************************

 Function:		lsearch_int
 Purpose:		Performs a search for value [key] in int [array]
			of size [size]
 Additional Notes:	None

************************************************************************/

int lsearch_int(int array [], int key, int size) {

	int n;

	for (n=0; n<size; n++) {

		if (array[n] == key) { 

			return n;	// returns array index where value was found

		} 

	} 

	return -1;			// otherwise exited normally but didn't find value

}

/************************************************************************

 Function:		median_filter
 Purpose:		Applies a median filter to a dataset
 Additional Notes:	

 This algorithm applies no padding to the start/end of the array.

************************************************************************/

int median_filter(double row_values [], double smoothed_row_values [], int nxelements, int median_half_filter_size) {

	int ii, jj;

	double this_iteration_filter_values [median_half_filter_size*2];

	for (ii=median_half_filter_size; ii<nxelements-median_half_filter_size; ii++) {
		
		memset(this_iteration_filter_values, 0, sizeof(double)*median_half_filter_size*2);

		for (jj=-median_half_filter_size; jj<median_half_filter_size; jj++) {

			this_iteration_filter_values[jj+median_half_filter_size] = row_values[ii+jj];
			
		}	
		double this_iteration_filter_values_sorted [median_half_filter_size*2];
		memcpy(this_iteration_filter_values_sorted, this_iteration_filter_values, sizeof(double)*median_half_filter_size*2);

		gsl_sort(this_iteration_filter_values_sorted, 1, median_half_filter_size*2);
		double this_iteration_filter_median = gsl_stats_median_from_sorted_data(this_iteration_filter_values_sorted, 1, median_half_filter_size*2);	

		// printf("%f\n", this_iteration_filter_median);	// DEBUG

		smoothed_row_values[ii] = this_iteration_filter_median;
	}

	return 0;

}

/************************************************************************

 Function:		populate_env_variable
 Purpose:		populates global variables with corresponding
			environment variables from the shell
 Additional Notes:	None

************************************************************************/

int populate_env_variable(char var_to_populate [], char env_var_name []) {

	if (!getenv(env_var_name)) {

		return 1;

	} else {

		strcpy(var_to_populate, getenv(env_var_name));

		return 0;

	}

}

/************************************************************************

 Function:		populate_img_parameters
 Purpose:		populates variables with corresponding parameters
			and prints to stdout
 Additional Notes:	None

************************************************************************/

int populate_img_parameters(char f [], fitsfile *f_ptr, int maxdim, int *bitpix, int *naxis, long naxes [], int *status, char title_of_img []) {

	if(!fits_get_img_param(f_ptr, maxdim, bitpix, naxis, naxes, status)) {

		printf("\nFile parameters");
		printf("\n---------------\n");
		printf("\n");
		printf("Frame name:\t\t%s\n", title_of_img);
		printf("Relative path:\t\t%s\n", f);
		printf("Bits per pixel:\t\t%d\n", *bitpix);
		printf("Number of axes:\t\t%d\n", *naxis);
		printf("First axis dimension:\t%ld\n", naxes[0]);
		printf("Second axis dimension:\t%ld\n", naxes[1]);

		return 0;
		
	} else {

		return 1;

	}

}

/************************************************************************

 Function:		print_file
 Purpose:		prints content of file to screen
 Additional Notes:	None

************************************************************************/

int print_file(char text_file [200]) {

	FILE *text_file_ptr;
	text_file_ptr = fopen(text_file, "r");

	char input_string [300];

	if (text_file_ptr) {
	
		while(!feof(text_file_ptr)) {

			memset(input_string, '\0', 300);

			fgets (input_string, 300, text_file_ptr);
	
			printf("%s", input_string);
	
		}

	} else {

		printf("\nWARNING:\tUnable to print file. File %s doesn't exist.\n", text_file);
		return 1;

	}

	return 0;

}

/************************************************************************

 Function:              strdup
 Purpose:               duplicates a string
 Additional Notes:      None

************************************************************************/

char *strdup(const char *str) {

	int n = strlen(str) + 1;
	char *dup = malloc(n);

	if(dup) {

		strcpy(dup, str);

	}

	return dup;

}

/************************************************************************

 Function:		write_additional_keys_file_to_header
 Purpose:		writes an additional keys file to a header
 Additional Notes:	None

************************************************************************/

int write_additional_keys_file_to_header(char ADDITIONAL_KEYS_FILE [], fitsfile *f_ptr, char operation [], int decimals, int *status) {

	FILE *file_ptr;
	file_ptr = fopen(ADDITIONAL_KEYS_FILE, "r");

	char input_string [300];

	char this_operation [300];
	char this_keyname [300];
	char this_string_value [300];
	double this_double_value;
	char this_comment [300];

	if (file_ptr) {
	
		while(!feof(file_ptr)) {

			memset(input_string, '\0', 300);
			fgets (input_string, 300, file_ptr);

			if (strncmp(input_string, "str", 3) == 0) { 		// we're dealing with a string

				sscanf(input_string, "%*s\t%s\t%s\t%[^\t]\t%[^\n]", this_operation, this_keyname, this_string_value, this_comment);

				if (strncmp(operation, this_operation, strlen(operation)) == 0) { 	// we need to insert this key

					if (!fits_update_key_str(f_ptr, this_keyname, this_string_value, this_comment, status)) {
		
					} else {

						return 1;

					}

				}

			} else if (strncmp(input_string, "dbl", 3) == 0) { 	// we're dealing with a double

				sscanf(input_string, "%*s\t%s\t%s\t%lf\t%[^\n]", this_operation, this_keyname, &this_double_value, this_comment);

				if (strncmp(operation, this_operation, strlen(operation)) == 0) { 	// we need to insert this key

					if (!fits_update_key_dbl(f_ptr, this_keyname, this_double_value, decimals, this_comment, status)) {
		
					} else {

						return 1;

					}

				}

			}
	
		}

	} else {

		return 1;

	}

	return 0;

}

/************************************************************************

 Function:		write_additional_key_to_file_dbl
 Purpose:		writes an additional key to file (double)
 Additional Notes:	None

************************************************************************/

int write_additional_key_to_file_dbl(char ADDITIONAL_KEYS_FILE [], char id [], char fits_key [], double fits_key_value, char fits_key_comment [], char ADDITIONAL_KEYS_FILE_WRITE_ACCESS []) {

	FILE *ADDITIONAL_KEYS_FILE_ptr;
	ADDITIONAL_KEYS_FILE_ptr = fopen(ADDITIONAL_KEYS_FILE, ADDITIONAL_KEYS_FILE_WRITE_ACCESS);

	if (ADDITIONAL_KEYS_FILE_ptr) {

		fprintf(ADDITIONAL_KEYS_FILE_ptr, "dbl\t%s\t%s\t%f\t%s\n", id, fits_key, fits_key_value, fits_key_comment);
		fclose(ADDITIONAL_KEYS_FILE_ptr);
		return 0;

	} else {

		printf("\nWARNING:\tUnable to write additional key to file. File %s doesn't exist.\n\n", ADDITIONAL_KEYS_FILE);
		fclose(ADDITIONAL_KEYS_FILE_ptr);
		return 1;

	}

}

/************************************************************************

 Function:		write_additional_key_to_file_str
 Purpose:		writes an additional key to file (string)
 Additional Notes:	None

************************************************************************/

int write_additional_key_to_file_str(char ADDITIONAL_KEYS_FILE [], char id [], char fits_key [], char fits_key_value [], char fits_key_comment [], char ADDITIONAL_KEYS_FILE_WRITE_ACCESS []) {

	FILE *ADDITIONAL_KEYS_FILE_ptr;
	ADDITIONAL_KEYS_FILE_ptr = fopen(ADDITIONAL_KEYS_FILE, ADDITIONAL_KEYS_FILE_WRITE_ACCESS);

	if (ADDITIONAL_KEYS_FILE_ptr) {

		fprintf(ADDITIONAL_KEYS_FILE_ptr, "str\t%s\t%s\t%s\t%s\n", id, fits_key, fits_key_value, fits_key_comment);
		fclose(ADDITIONAL_KEYS_FILE_ptr);
		return 0;

	} else {

		printf("\nWARNING:\tUnable to write additional key to file. File %s doesn't exist.\n\n", ADDITIONAL_KEYS_FILE);
		fclose(ADDITIONAL_KEYS_FILE_ptr);
		return 1;

	}

}

