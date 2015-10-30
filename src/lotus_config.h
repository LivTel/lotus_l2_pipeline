char HEADER_FILE [100];
int RETURN_FLAG					= 0;

bool INDEXING_CORRECTION			= TRUE;

char FILE_WRITE_ACCESS [2] 			= "w";				// r readonly; w overwrite; a+ append

int IMG_READ_ACCURACY				= 82;
int INTERMEDIATE_IMG_ACCURACY[2]		= {-32, 82};			// Definition follows types from http://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node20.html {create_img, write_img}

char FITS_KEYS_TO_OMIT [100];							// variable to hold location of FITS_KEYS_TO_OMIT file

