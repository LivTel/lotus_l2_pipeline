##################################
# Makefile 			 #
# Last Modified Date:   02/09/15 #
##################################

#############################
# compilers, archivers etc. #
#############################

CC		= /usr/bin/gcc
AR		= /usr/bin/ar

##################
# compiler flags #
##################

CFLAGS		= -Wall -g -std=c99

#########################
# csh terminal commands #
#########################

CAT		= /bin/cat
GREP		= /bin/grep
MV		= /bin/mv
RM		= /bin/rm -i

################
# $(LIBS) flag #
################

LIBS		= -L$(L2_LIB_DIR)#						rmb/desktop,laptop
#LIBS		= -L/space/home/dev/src/cfitsio3310 -L$(L2_LIB_DIR)#		lt/ltdevsrv

####################
# $(INCLUDES) flag #
####################

INCLUDES	= -I/usr/include/cfitsio -I/usr/include/gsl#			rmb/desktop,laptop
#INCLUDES	= -I/space/home/dev/src/cfitsio3310 -I/usr/include/gsl#		lt/ltdevsrv

################
# make default #
################

all: lotus_red_clip lotus_red_find lotus_red_trace_sdist lotus_red_correct_sdist lotus_red_rebin lotus_red_extract lotus_red_reformat

#################
# make routines #
#################

lotus_red_clip: lotus_red_clip.c liblotus
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) lotus_red_clip.c -llotus -lgsl -lgslcblas -lcfitsio -lm -o loclip
	$(MV) loclip $(L2_BIN_DIR)
	
lotus_red_find: lotus_red_find.c liblotus
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) lotus_red_find.c -llotus -lgsl -lgslcblas -lcfitsio -lm -o lofind
	$(MV) lofind $(L2_BIN_DIR)
	
lotus_red_trace_sdist: lotus_red_trace_sdist.c liblotus
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) lotus_red_trace_sdist.c -llotus -lgsl -lgslcblas -lcfitsio -lm -o lotrace
	$(MV) lotrace $(L2_BIN_DIR)	
	
lotus_red_correct_sdist: lotus_red_correct_sdist.c liblotus
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) lotus_red_correct_sdist.c -llotus -lgsl -lgslcblas -lcfitsio -lm -o locorrect
	$(MV) locorrect $(L2_BIN_DIR)	
	
lotus_red_rebin: lotus_red_rebin.c liblotus
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) lotus_red_rebin.c -llotus -lgsl -lgslcblas -lcfitsio -lm -o lorebin
	$(MV) lorebin $(L2_BIN_DIR)	
	
lotus_red_extract: lotus_red_extract.c liblotus
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) lotus_red_extract.c -llotus -lgsl -lgslcblas -lcfitsio -lm -o loextract
	$(MV) loextract $(L2_BIN_DIR)	
	
lotus_red_reformat: lotus_red_reformat.c liblotus
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) lotus_red_reformat.c -llotus -lgsl -lgslcblas -lcfitsio -lm -o loreformat
	$(MV) loreformat $(L2_BIN_DIR)		

lotus_functions: lotus_functions.c
	$(CC) $(CFLAGS) $(INCLUDES) -c lotus_functions.c

lotus_error_handling: lotus_error_handling.c
	$(CC) $(CFLAGS) $(INCLUDES) -c lotus_error_handling.c

############
# liblotus #
############

liblotus: lotus_functions lotus_error_handling
	$(AR) rcs liblotus.a lotus_functions.o lotus_error_handling.o
	$(MV) liblotus.a $(L2_LIB_DIR)
	$(RM) -f *.o

##############
# make clean #
##############

clean:
	$(RM) $(L2_LIB_DIR)/*
	$(RM) $(L2_BIN_DIR)/*

