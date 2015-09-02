from matplotlib import pyplot as plt
import pyfits
import sys
import numpy as np
from optparse import OptionParser
import os
from subprocess import Popen, PIPE
from shutil import copyfile, move
import time
from datetime import date
import ConfigParser

from L2_exec_err import errors
from L2_analyse_plotpeaks import execute as a_pp_execute
from L2_analyse_plotimage import execute as a_pi_execute
from L2_analyse_plotspec import execute as a_ps_execute

L2_BIN_DIR 	        = os.environ['L2_BIN_DIR']
L2_TEST_DIR 	        = os.environ['L2_TEST_DIR']
L2_SCRIPT_DIR	        = os.environ['L2_SCRIPT_DIR']
L2_MAN_DIR	        = os.environ['L2_MAN_DIR']
L2_CONFIG_DIR           = os.environ['L2_CONFIG_DIR']
L2_INI_DIR              = os.environ['L2_INI_DIR']
L2_LOOKUP_TABLES_DIR    = os.environ['L2_LOOKUP_TABLES_DIR']

clip            = L2_BIN_DIR + "/loclip"
find            = L2_BIN_DIR + "/lofind"
trace           = L2_BIN_DIR + "/lotrace"
correct         = L2_BIN_DIR + "/locorrect"
extract         = L2_BIN_DIR + "/loextract"    
rebin           = L2_BIN_DIR + "/lorebin"
reformat        = L2_BIN_DIR + "/loreformat"
    
plot_peaks      = L2_SCRIPT_DIR + "/L2_analyse_plotpeaks.py"
plot_image      = L2_SCRIPT_DIR + "/L2_analyse_plotimage.py"
plot_spec       = L2_SCRIPT_DIR + "/L2_analyse_plotspec.py"    

def print_header():
    with open(L2_MAN_DIR + "/HEADER") as f:
        for line in f:
	    print line.strip('\n')
    
def print_routine(routine):
    bar = []
    for i in range(len(routine)):
        bar.append("*")
    print ''.join(bar) + '****'
    print '* ' + routine + ' *'
    print ''.join(bar) + '****' 
    
def print_notification(message):
    print "* " + message
    print
  
def read_ini(err, path):
    if not os.path.exists(path):
        err.set_code(28)
    ini = ConfigParser.ConfigParser()
    ini.read(path)
    cfg = {}
    for section in ini.sections():
        cfg[section] = {}
        for option in ini.options(section):
            cfg[section][option] = str(ini.get(section, option))  
    return cfg
  
def rename_dat_files(suffix):
    for i in os.listdir("."):
        if i.endswith(".dat"):
            if not i.startswith("p_"):
                filename = os.path.splitext(os.path.basename(i))[0]
                ext = os.path.splitext(os.path.basename(i))[1]
                move(i, "p_" + filename + "_" + suffix + ext)    
  
def rewrite_error_codes_file(error_codes_file, old_header_key="", new_header_key="", add_to_desc="", omit=False):
    with open("new_error_codes", "w") as f_new:
        with open(error_codes_file) as f_old:
            for line in f_old:
                if omit:
                    if not line.startswith(old_header_key):
                        f_new.write(line)
                else:
                    if not line.startswith(old_header_key):
                        f_new.write(line)
                    else:
                        key         = line.split('\t')[0]
                        code        = line.split('\t')[1]   
                        desc        = line.split('\t')[2].strip('\n')       
                        f_new.write(new_header_key + "\t" + code + "\t" + desc + " " + add_to_desc + "\n")
    os.remove(error_codes_file)
    move("new_error_codes", error_codes_file) 
                
def search_lookup_table(lookup_table_path, this_DATEOBS, this_CCDXBIN, this_CCDYBIN):
  
    a_path = []
    a_from_date = []
    a_from_time = []
    a_to_date = []
    a_to_time = []
    a_binning = []
    with open(lookup_table_path) as f:
        for line in f:
            if line != "\n":
                this_path = line.split('\t')[0].strip('\n').strip()
                this_binning = line.split('\t')[1].strip('\n').strip()
                this_from_date = line.split('\t')[2].strip('\n').strip()
                this_from_time = line.split('\t')[3].strip('\n').strip()
                a_path.append(this_path)
                a_from_date.append(this_from_date)
                a_from_time.append(this_from_time)
                a_binning.append(this_binning)
                
                this_to_date = line.split('\t')[4].strip('\n').strip()
                if ("now" in this_to_date):
                    today = date.today()
                    this_to_date = today.strftime("%d/%m/%y")
                    this_to_time = time.strftime("%H:%M:%S")
                else:
                    this_to_date = line.split('\t')[4].strip('\n').strip()
                    this_to_time = line.split('\t')[5].strip('\n').strip()
                
                a_to_date.append(this_to_date)
                a_to_time.append(this_to_time)
                 
    this_file_datetime   = time.strptime(this_DATEOBS, "%Y-%m-%dT%H:%M:%S")
    this_file_binning    = this_CCDXBIN + "x" + this_CCDYBIN
    
    chosen_entry = None
    for i in range(len(a_path)):
        this_entry_from_time = time.strptime(a_from_date[i] + " " + a_from_time[i], "%d/%m/%y %H:%M:%S")
        this_entry_to_time = time.strptime(a_to_date[i] + " " + a_to_time[i], "%d/%m/%y %H:%M:%S")
        this_entry_binning = a_binning[i]
        
        if this_file_datetime >= this_entry_from_time and this_file_datetime <= this_entry_to_time and this_file_binning == this_entry_binning:
            chosen_entry = a_path[i]
            break     
            
    return chosen_entry
    
def chk_ref_run(f_ref, f_cont):
  
     # get basename of files
    ref                 = os.path.splitext(os.path.basename(f_ref))[0]
    cont                = os.path.splitext(os.path.basename(f_cont))[0]
  
    out_ref_tr_filename  = "tmp.fits"  
    out_ref_cor_filename = "tmp2.fits"      
    ref_pre_sdist_plot   = "ref_pre_sdist_plot.png"   
    ref_post_sdist_plot  = "ref_post_sdist_plot.png"
  
    err = errors()

    # input sanity checks
    if not all([f_ref, f_cont]):
        err.set_code(1)
    elif not os.path.exists(f_ref):
        err.set_code(3)   
    elif not os.path.exists(f_cont):
        err.set_code(4) 
        
    # move files to working directory, redefine paths and change to working directory
    try:
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
            copyfile(f_ref, work_dir + "/" + ref + ref_suffix + ".fits")
            copyfile(f_cont, work_dir + "/" + cont + cont_suffix + ".fits")
        else:
            if clobber:
                for i in os.listdir(work_dir):
                    os.remove(work_dir + "/" + i)
                copyfile(f_ref, work_dir + "/" + ref + ref_suffix + ".fits")
                copyfile(f_cont, work_dir + "/" + cont + cont_suffix + ".fits")
            else:
                err.set_code(14)
    except OSError:
        err.set_code(15)
    
    f_ref = ref + ref_suffix + ".fits"
    f_cont = cont + cont_suffix + ".fits"

    os.chdir(work_dir)        
        
    # -----------------------------------
    # - DETERMINE SUITABLE CONFIG FILES -
    # -----------------------------------
    
    print_routine("Finding suitable config file.")
    print 
    
    f_ref_fits = pyfits.open(f_ref)
    try:
        f_ref_fits_hdr_DATEOBS = f_ref_fits[0].header['DATE-OBS'].strip()
        f_ref_fits_hdr_CCDXBIN = str(f_ref_fits[0].header['CCDXBIN']).strip()
        f_ref_fits_hdr_CCDYBIN = str(f_ref_fits[0].header['CCDYBIN']).strip()
    except KeyError:
        f_ref_fits.close()
        err.set_code(16)
    f_ref_fits.close()        
        
    config_tab_path      = L2_LOOKUP_TABLES_DIR + "/" + "/config.tab"
   
    chosen_config_entry = search_lookup_table(config_tab_path, f_ref_fits_hdr_DATEOBS, f_ref_fits_hdr_CCDXBIN, f_ref_fits_hdr_CCDYBIN)
    if chosen_config_entry is None:
        print_notification("Failed.") 
        err.set_code(29)
        
    chosen_config_file_path = L2_INI_DIR + "/" + chosen_config_entry

    print_notification("Success. Using file " + chosen_config_file_path)    
    cfg = read_ini(err, chosen_config_file_path)            
        
    # -------------------------
    # - CLIP SPECTRA (LOCLIP) -
    # -------------------------
    
    print_routine("Trim spectra (lotrim)")
    
    in_ref_filename = f_ref
    in_cont_filename = f_cont

    output = Popen([clip, in_cont_filename, in_ref_filename, cfg['loclip']['bin_size_px'], cfg['loclip']['bg_percentile'], cfg['loclip']['clip_sigma'], cfg['loclip']['thresh_sigma'], \
      cfg['loclip']['scan_window_size_px'], cfg['loclip']['scan_window_nsigma'], cfg['loclip']['min_spectrum_width_px'], cfg['loclip']['force_bottom_px'], cfg['loclip']['force_top_px'], \
        out_ref_tr_filename], stdout=PIPE)   
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(20, is_fatal=False)
    
    # ---------------------------------------------
    # - FIND PEAKS OF REFERENCE SPECTRUM (LOFIND) -
    # ---------------------------------------------
    
    print_routine("Find peaks of reference spectrum (lofind)")   
    
    # use previous file as input
    in_ref_filename = out_ref_tr_filename

    output = Popen([find, in_ref_filename, cfg['lofind_ref']['bin_size_px'], cfg['lofind_ref']['bg_percentile'], cfg['lofind_ref']['clip_sigma'], \
      cfg['lofind_ref']['median_filter_width_px'], cfg['lofind_ref']['min_snr'], cfg['lofind_ref']['min_spatial_width_px'], cfg['lofind_ref']['finding_window_lo_px'], \
        cfg['lofind_ref']['finding_window_hi_px'], cfg['lofind_ref']['max_centering_num_px'], cfg['lofind_ref']['centroid_half_window_size_px'], \
          cfg['lofind_ref']['min_used_bins']], stdout=PIPE)
    print output.stdout.read()  
    output.wait()
    if output.returncode != 0:
        err.set_code(21, is_fatal=False)   
        
    # ----------------------------------------------
    # - FIND SDIST OF REFERENCE SPECTRUM (LOTRACE) -
    # ----------------------------------------------
    
    print_routine("Find sdist of reference spectrum (lotrace)") 
    
    output = Popen([trace, cfg['lotrace']['polynomial_order']], stdout=PIPE)
    print output.stdout.read()  
    output.wait()
    if output.returncode != 0:
        err.set_code(22, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATTR", "L2STATT1", add_to_desc="(ref uncorrected)")    
    
    # ---------------------------------------------------------
    # - PLOT TRACE OF REFERENCE SPECTRUM PRE SDIST CORRECTION -
    # ---------------------------------------------------------

    print_routine("Plot trace of reference spectrum pre sdist correction") 
    print
    
    try:
        a_pp_execute(in_ref_filename, cfg['general']['lofind_output_file'], cfg['general']['lotrace_output_file'], ref_pre_sdist_plot, "Reference pre SDIST correction", \
        cfg['general']['max_curvature_post_cor']) 
    except IOError:
        pass
    
    if os.path.exists(ref_pre_sdist_plot):
        print_notification("Success.")
    else:
        print_notification("Failed.") 
        err.set_code(6, is_fatal=False)
        
    # -----------------------------------------
    # - CORRECT SPECTRA FOR SDIST (LOCORRECT) -
    # -----------------------------------------
    
    print_routine("Correct spectra for sdist (locorrect)")      

    output = Popen([correct, in_ref_filename, cfg['locorrect']['interpolation_type'], cfg['locorrect']['conserve_flux'], out_ref_cor_filename], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(23, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATCO", "L2STATOT", add_to_desc="(ref)")
    
    # ---------------------------------------------
    # - FIND PEAKS OF REFERENCE SPECTRUM (LOFIND) -
    # ---------------------------------------------
    
    print_routine("Find peaks of reference spectrum (lofind)")       
    
    in_ref_filename = out_ref_cor_filename

    output = Popen([find, in_ref_filename, cfg['lofind_ref']['bin_size_px'], cfg['lofind_ref']['bg_percentile'], cfg['lofind_ref']['clip_sigma'], \
      cfg['lofind_ref']['median_filter_width_px'], cfg['lofind_ref']['min_snr'], cfg['lofind_ref']['min_spatial_width_px'], cfg['lofind_ref']['finding_window_lo_px'], \
        cfg['lofind_ref']['finding_window_hi_px'], cfg['lofind_ref']['max_centering_num_px'], cfg['lofind_ref']['centroid_half_window_size_px'], \
          cfg['lofind_ref']['min_used_bins']], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(21, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATFI", omit=True)    
    
    # ----------------------------------------------
    # - FIND SDIST OF REFERENCE SPECTRUM (LOTRACE) -
    # ----------------------------------------------
    
    print_routine("Find sdist of reference spectrum (lotrace)")       
    
    output = Popen([trace, cfg['lotrace']['polynomial_order']], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(22, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATTR", omit=True)

    # ----------------------------------------------------------------------------------------------------------
    # - PLOT TRACE OF SPECTRUM POST SDIST CORRECTION AND CHECK RETURN CODE FOR SIGNIFICANT REMAINING CURVATURE -
    # ----------------------------------------------------------------------------------------------------------
    
    print_routine("Plot trace of spectrum post sdist correction and check for curvature (l2pp)") 
    print
    
    try:
        rtn = a_pp_execute(in_ref_filename, cfg['general']['lofind_output_file'], cfg['general']['lotrace_output_file'], ref_post_sdist_plot, "Reference post SDIST correction", \
        cfg['general']['max_curvature_post_cor']) 
    except IOError:
        rtn = -1
    
    if os.path.exists(ref_post_sdist_plot) and rtn == 0:
        print_notification("Success.")
    elif not os.path.exists(ref_post_sdist_plot):
        print_notification("Failed.") 
        err.set_code(7, is_fatal=False)   
    elif rtn != 0:
        print_notification("Failed.") 
        err.set_code(8)     
    
    # make sure no negative routine error codes for essential routines
    rtn_keys    = []    
    rtn_codes   = []
    with open(cfg['general']['error_codes_file']) as f:
        for line in f:
            if line.startswith("L2"):
                rtn_keys.append(str(line.split('\t')[0]))
                rtn_codes.append(int(line.split('\t')[1]))

    for idx, i in enumerate(rtn_keys):
        if i == "L2STATF2" or i == "L2STATX1" or i=="L2STATX2":
            continue                                                    # skip 1D extraction related error codes.
        if rtn_codes[idx] < 0: 
            err.set_code(13)             
             
    err.set_code(0)   	# this is a bit of a bodge, it disregards the current error code!  
    

def full_run(f_target, f_ref, f_cont, work_dir, clobber):

    err = errors()

    # input sanity checks
    if not all([f_target, f_ref, f_cont]):
        err.set_code(1)
    elif not os.path.exists(f_target):
        err.set_code(2)   
    elif not os.path.exists(f_ref):
        err.set_code(3)   
    elif not os.path.exists(f_cont):	
        err.set_code(4)    
      
    # get basename of files
    target              = os.path.splitext(os.path.basename(f_target))[0]
    ref                 = os.path.splitext(os.path.basename(f_ref))[0]
    cont                = os.path.splitext(os.path.basename(f_cont))[0]          

    # OUTPUT
    ## L2
    output_target       = target[:-1] + "2.fits"       
    
    ## plots
    ref_pre_sdist_plot  	= "ref_pre_sdist_plot.png"
    ref_post_sdist_plot 	= "ref_post_sdist_plot.png"
    L1_IMAGE_plot	        = "L1_IMAGE.png"
    SPEC_NONSS_plot	        = os.path.splitext(os.path.basename(output_target))[0] + "_SPEC_NONSS.png"
    SPEC_SS_plot                = os.path.splitext(os.path.basename(output_target))[0] + "_SPEC_SS.png"
    montage_plot                = os.path.splitext(os.path.basename(output_target))[0] + "_output.png"  

    # move files to working directory, redefine paths and change to working directory
    try:
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
            copyfile(f_target, work_dir + "/" + target + target_suffix + ".fits")
            copyfile(f_ref, work_dir + "/" + ref + ref_suffix + ".fits")
            copyfile(f_cont, work_dir + "/" + cont + cont_suffix + ".fits") 
        else:
            if clobber:
                for i in os.listdir(work_dir):
	            os.remove(work_dir + "/" + i)
                copyfile(f_target, work_dir + "/" + target + target_suffix + ".fits")
                copyfile(f_ref, work_dir + "/" + ref + ref_suffix + ".fits")
                copyfile(f_cont, work_dir + "/" + cont + cont_suffix + ".fits")
            else:
	        err.set_code(14)
    except OSError:
        err.set_code(15)
    
    f_target = target + target_suffix + ".fits"
    f_ref = ref + ref_suffix + ".fits"
    f_cont = cont + cont_suffix + ".fits"

    os.chdir(work_dir)
    
    # ---------
    # - START -
    # ---------  
    
    today = date.today()
    now_date = today.strftime("%d-%m-%y")
    now_time = time.strftime("%H:%M:%S")
    # add L2DATE key to additional_keys
    with open("additional_keys", 'w') as f:
        f.write("str\tSTARTDATE\tL2DATE\t" + now_date + " " + now_time + "\twhen this reduction was performed\n")     
        
    st_unix = time.time()    
    
    # ----------------------------------
    # - DETERMINE SUITABLE CONFIG FILE -
    # ----------------------------------
    
    print_routine("Finding suitable config files.")
    print 
    
    f_ref_fits = pyfits.open(f_ref)
    try:
        f_ref_fits_hdr_DATEOBS = f_ref_fits[0].header['DATE-OBS'].strip()
        f_ref_fits_hdr_CCDXBIN = str(f_ref_fits[0].header['CCDXBIN']).strip()
        f_ref_fits_hdr_CCDYBIN = str(f_ref_fits[0].header['CCDYBIN']).strip()       
    except KeyError:
        f_ref_fits.close()
        err.set_code(16)
    f_ref_fits.close()          
        
    config_tab_path      = L2_LOOKUP_TABLES_DIR + "/" + "/config.tab"
    
    chosen_config_entry = search_lookup_table(config_tab_path, f_ref_fits_hdr_DATEOBS, f_ref_fits_hdr_CCDXBIN, f_ref_fits_hdr_CCDYBIN)
    if chosen_config_entry is None:
        print_notification("Failed.") 
        err.set_code(29)
        
    chosen_config_file_path = L2_INI_DIR + "/" + "/" + chosen_config_entry

    print_notification("Success. Using file " + chosen_config_file_path)    
    cfg = read_ini(err, chosen_config_file_path)    
         
    # -------------------------
    # - CLIP SPECTRA (LOCLIP) -
    # -------------------------
    
    print_routine("Trim spectra (lotrim)")
    
    in_target_filename = f_target
    in_ref_filename = f_ref
    in_cont_filename = f_cont
    out_target_filename = target + target_suffix + trim_suffix + ".fits"
    out_ref_filename = ref + ref_suffix + trim_suffix + ".fits"
    out_cont_filename = cont + cont_suffix + trim_suffix + ".fits"

    ## target
    output = Popen([clip, in_cont_filename, in_target_filename, cfg['loclip']['bin_size_px'], cfg['loclip']['bg_percentile'], cfg['loclip']['clip_sigma'], cfg['loclip']['thresh_sigma'], \
      cfg['loclip']['scan_window_size_px'], cfg['loclip']['scan_window_nsigma'], cfg['loclip']['min_spectrum_width_px'], cfg['loclip']['force_bottom_px'], cfg['loclip']['force_top_px'], \
        out_target_filename], stdout=PIPE)  
    print output.stdout.read()
    output.wait()
    if output.returncode != 0:
        err.set_code(20, is_fatal=False)    
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATCL", "L2STATCT", add_to_desc="(target)")

    ## reference
    output = Popen([clip, in_cont_filename, in_ref_filename, cfg['loclip']['bin_size_px'], cfg['loclip']['bg_percentile'], cfg['loclip']['clip_sigma'], cfg['loclip']['thresh_sigma'], \
      cfg['loclip']['scan_window_size_px'], cfg['loclip']['scan_window_nsigma'], cfg['loclip']['min_spectrum_width_px'], cfg['loclip']['force_bottom_px'], cfg['loclip']['force_top_px'], \
        out_ref_filename], stdout=PIPE)   
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(20, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATCL", "L2STATCR", add_to_desc="(ref)")
    
    ## continuum
    output = Popen([clip, in_cont_filename, in_cont_filename, cfg['loclip']['bin_size_px'], cfg['loclip']['bg_percentile'], cfg['loclip']['clip_sigma'], cfg['loclip']['thresh_sigma'], \
      cfg['loclip']['scan_window_size_px'], cfg['loclip']['scan_window_nsigma'], cfg['loclip']['min_spectrum_width_px'], cfg['loclip']['force_bottom_px'], cfg['loclip']['force_top_px'], \
        out_cont_filename], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(20, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATCL", "L2STATCC", add_to_desc="(continuum)")

    # ---------------------------------------------
    # - FIND PEAKS OF REFERENCE SPECTRUM (LOFIND) -
    # ---------------------------------------------
    
    print_routine("Find peaks of reference spectrum (lofind)")   
    
    in_ref_filename = ref + ref_suffix + trim_suffix + ".fits"

    output = Popen([find, in_ref_filename, cfg['lofind_ref']['bin_size_px'], cfg['lofind_ref']['bg_percentile'], cfg['lofind_ref']['clip_sigma'], \
      cfg['lofind_ref']['median_filter_width_px'], cfg['lofind_ref']['min_snr'], cfg['lofind_ref']['min_spatial_width_px'], cfg['lofind_ref']['finding_window_lo_px'], \
        cfg['lofind_ref']['finding_window_hi_px'], cfg['lofind_ref']['max_centering_num_px'], cfg['lofind_ref']['centroid_half_window_size_px'], \
          cfg['lofind_ref']['min_used_bins']], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(21, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATFI", "L2STATF1", add_to_desc="(ref uncorrected)")

    # ----------------------------------------------
    # - FIND SDIST OF REFERENCE SPECTRUM (LOTRACE) -
    # ----------------------------------------------
    
    print_routine("Find sdist of reference spectrum (lotrace)") 
    
    output = Popen([trace, cfg['lotrace']['polynomial_order']], stdout=PIPE)
    print output.stdout.read()  
    output.wait()
    if output.returncode != 0:
        err.set_code(22, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATTR", "L2STATT1", add_to_desc="(ref uncorrected)")
  
    # ---------------------------------------------------------
    # - PLOT TRACE OF REFERENCE SPECTRUM PRE SDIST CORRECTION -
    # ---------------------------------------------------------
    
    print_routine("Plot trace of reference spectrum pre sdist correction") 
    print
    
    in_ref_filename = ref + ref_suffix + trim_suffix + ".fits"
    
    try:
        a_pp_execute(in_ref_filename, cfg['general']['lofind_output_file'], cfg['general']['lotrace_output_file'], ref_pre_sdist_plot, "Reference pre SDIST correction", \
      cfg['general']['max_curvature_post_cor']) 
    except IOError:
        pass

    if os.path.exists(ref_pre_sdist_plot):
        print_notification("Success.")
    else:
        print_notification("Failed.") 
        err.set_code(6, is_fatal=False)

    # -----------------------------------------
    # - CORRECT SPECTRA FOR SDIST (LOCORRECT) -
    # -----------------------------------------
    
    print_routine("Correct spectra for sdist (locorrect)")      

    in_target_filename = target + target_suffix + trim_suffix + ".fits"
    in_ref_filename = ref + ref_suffix + trim_suffix + ".fits"
    in_cont_filename = cont + cont_suffix + trim_suffix + ".fits"

    out_target_filename = target + target_suffix + trim_suffix + cor_suffix + ".fits"
    out_ref_filename = ref + ref_suffix + trim_suffix + cor_suffix + ".fits"
    out_cont_filename = cont + cont_suffix + trim_suffix + cor_suffix + ".fits"

    ## target
    output = Popen([correct, in_target_filename, cfg['locorrect']['interpolation_type'], cfg['locorrect']['conserve_flux'], out_target_filename], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(23, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATCO", "L2STATOT", add_to_desc="(target)")
    
    ## reference
    output = Popen([correct, in_ref_filename, cfg['locorrect']['interpolation_type'], cfg['locorrect']['conserve_flux'], out_ref_filename], stdout=PIPE)
    print output.stdout.read()  
    output.wait()
    if output.returncode != 0:
        err.set_code(23, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATCO", "L2STATOR", add_to_desc="(ref)")
    
    ## continuum
    output = Popen([correct, in_cont_filename, "linear", "1", out_cont_filename], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(23, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATCO", "L2STATOC", add_to_desc="(continuum)")
   
    # rename dat files to avoid conflict
    rename_dat_files("ref_uncorrected")

    # ---------------------------------------------
    # - FIND PEAKS OF REFERENCE SPECTRUM (LOFIND) -
    # ---------------------------------------------
    
    print_routine("Find peaks of reference spectrum (lofind)")       
    
    in_ref_filename = ref + ref_suffix + trim_suffix + cor_suffix + ".fits"

    output = Popen([find, in_ref_filename, cfg['lofind_ref']['bin_size_px'], cfg['lofind_ref']['bg_percentile'], cfg['lofind_ref']['clip_sigma'], \
      cfg['lofind_ref']['median_filter_width_px'], cfg['lofind_ref']['min_snr'], cfg['lofind_ref']['min_spatial_width_px'], cfg['lofind_ref']['finding_window_lo_px'], \
        cfg['lofind_ref']['finding_window_hi_px'], cfg['lofind_ref']['max_centering_num_px'], cfg['lofind_ref']['centroid_half_window_size_px'], \
          cfg['lofind_ref']['min_used_bins']], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(21, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATFI", omit=True)

    # ----------------------------------------------
    # - FIND SDIST OF REFERENCE SPECTRUM (LOTRACE) -
    # ----------------------------------------------
    
    print_routine("Find sdist of reference spectrum (lotrace)")       
    
    output = Popen([trace, cfg['lotrace']['polynomial_order']], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(22, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATTR", omit=True)

    # ----------------------------------------------------------------------------------------------------------
    # - PLOT TRACE OF SPECTRUM POST SDIST CORRECTION AND CHECK RETURN CODE FOR SIGNIFICANT REMAINING CURVATURE -
    # ----------------------------------------------------------------------------------------------------------
    
    print_routine("Plot trace of spectrum post sdist correction and check for curvature (l2pp)") 
    print
    
    in_ref_filename = ref + ref_suffix + trim_suffix + cor_suffix + ".fits"
    
    try:
        rtn = a_pp_execute(in_ref_filename, cfg['general']['lofind_output_file'], cfg['general']['lotrace_output_file'], ref_post_sdist_plot, "Reference post SDIST correction", \
        cfg['general']['max_curvature_post_cor']) 
    except IOError:
        rtn = -1
    
    if os.path.exists(ref_post_sdist_plot) and rtn == 0:
        print_notification("Success.")
    elif not os.path.exists(ref_post_sdist_plot):
        print_notification("Failed.") 
        err.set_code(7, is_fatal=False)   
    elif rtn != 0:
        print_notification("Failed.") 
        err.set_code(8) 
        
    # rename dat files to avoid conflict
    rename_dat_files("ref_corrected")       
                
    # -------------------------------------------------
    # - FIND PIXEL TO WAVELENGTH SOLUTIONS (loarcfit) -
    # -------------------------------------------------
    
    print_routine("Copy over hardcoded dispersion solution from solutions file")   

    copyfile(L2_CONFIG_DIR + "/ARC_SOLUTION.dat", "./loarcfit_wavfits.dat")

    print_notification("\nSuccess.")

    # -------------------
    # - REBIN (lorebin) -
    # -------------------
    
    print_routine("Rebin data spectrally (lorebin)") 
    
    in_target_filename = target + target_suffix + trim_suffix + cor_suffix + ".fits"
    out_target_filename = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ".fits"

    output = Popen([rebin, in_target_filename, cfg['lorebin']['start_wav'], cfg['lorebin']['end_wav'], cfg['lorebin']['interpolation_type'], cfg['lorebin']['dispersion'], \
       cfg['lorebin']['conserve_flux'], out_target_filename], stdout=PIPE)
    print output.stdout.read()  
    output.wait()
    if output.returncode != 0:
        err.set_code(25, is_fatal=False)     
    
    # rename dat files to avoid conflict
    rename_dat_files("arc_corrected")         
                
    # ---------------------------------------------------------------
    # - FIND POSITION OF TARGET SPECTRUM WITH A SINGLE BIN (LOFIND) -
    # ---------------------------------------------------------------
    print_routine("Find peaks of target spectrum (lofind)")        
    in_target_filename = target + target_suffix + trim_suffix + cor_suffix + ".fits"

    output = Popen([find, in_target_filename, cfg['lofind_target']['bin_size_px'], cfg['lofind_target']['bg_percentile'], cfg['lofind_target']['clip_sigma'], \
      cfg['lofind_target']['median_filter_width_px'], cfg['lofind_target']['min_snr'], cfg['lofind_target']['min_spatial_width_px'], cfg['lofind_target']['finding_window_lo_px'], \
        cfg['lofind_target']['finding_window_hi_px'], cfg['lofind_target']['max_centering_num_px'], cfg['lofind_target']['centroid_half_window_size_px'], \
          cfg['lofind_target']['min_used_bins']], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(21, is_fatal=False) 
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATFI", "L2STATF2", "(target corrected)")      
    
    # -------------------------------
    # - EXTRACT SPECTRA (LOEXTRACT) -
    # -------------------------------

    print_routine("Extract NONSS spectra (loextract)")
    
    in_target_filename = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ".fits"
    out_target_filename = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ext_suffix + ".fits"

    output = Popen([extract, in_target_filename, cfg['loextract_nonss']['method'], cfg['loextract_nonss']['ss_method'], cfg['loextract_nonss']['target_half_aperture_px'], \
      cfg['loextract_nonss']['sky_window_half_aperture_px'], out_target_filename], stdout=PIPE)
    print output.stdout.read() 
    output.wait()
    if output.returncode != 0:
        err.set_code(26, is_fatal=False)     
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATEX", "L2STATX1", "(target NONSS)") 
    
    print_routine("Extract SS spectra (loextract)")
    
    in_target_filename = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ".fits"
    out_target_filename = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ext_suffix + ss_suffix + ".fits"
    
    output = Popen([extract, in_target_filename, cfg['loextract_ss']['method'], cfg['loextract_ss']['ss_method'], cfg['loextract_ss']['target_half_aperture_px'], \
      cfg['loextract_ss']['sky_window_half_aperture_px'], out_target_filename], stdout=PIPE)
    print output.stdout.read()
    output.wait()
    if output.returncode != 0:
        err.set_code(26, is_fatal=False)      
    rewrite_error_codes_file(cfg['general']['error_codes_file'], "L2STATEX", "L2STATX2", "(target SS)")     
 
    # rename dat files to avoid conflict
    rename_dat_files("target_corrected")
            
    # ------------------------------
    # - REFORMAT FILE (LOREFORMAT) -
    # ------------------------------
    
    print_routine("Reformat spectra (loreformat)")
    
    in_target_headers_filename = target + target_suffix + ".fits"    
    in_target_filename_L1_IMAGE = target + target_suffix + ".fits"    
    in_target_filename_LSS_NONSS = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ".fits"
    in_target_filename_SPEC_NONSS = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ext_suffix + ".fits"
    in_target_filename_SPEC_SS = target + target_suffix + trim_suffix + cor_suffix + reb_suffix + ext_suffix + ss_suffix + ".fits"

    out_target_filename = target[:-1] + "2.fits"
    
    ## make L1 extension
    output = Popen([reformat, in_target_filename_L1_IMAGE, in_target_headers_filename, "L1_IMAGE", out_target_filename], stdout=PIPE)
    print output.stdout.read()   
    
    ## make additional extensions
    for op in cfg['loreformat']['operations'].split(','):
        if op.strip() == "LSS_NONSS":
          in_f = in_target_filename_LSS_NONSS
        elif op.strip() == "SPEC_NONSS":
          in_f = in_target_filename_SPEC_NONSS
        elif op.strip() == "SPEC_SS":
          in_f = in_target_filename_SPEC_SS
        else:
          continue
        
        output = Popen([reformat, in_f, in_target_headers_filename, op.strip(), out_target_filename], stdout=PIPE)
        print output.stdout.read() 
        output.wait()
        if output.returncode != 0:
            err.set_code(27, is_fatal=False)  
    
    # ----------------------------------------------
    # - GENERATE RASTER PLOT OF L1_IMAGE extension -
    # ----------------------------------------------   
    
    print_routine("Plot extensions of output file (l2pi)")  
    print
    
    # use previously defined output filename from loreformat
    in_target_filename = out_target_filename
    
    a_pi_execute(in_target_filename, "L1_IMAGE", L1_IMAGE_plot, "L1_IMAGE")
    
    if os.path.exists(L1_IMAGE_plot):
        print_notification("Success.")
    else:
        print_notification("Failed.") 
        err.set_code(9, is_fatal=False)    
        
    # ------------------------------------------------
    # - GENERATE RASTER PLOT OF SPEC_NONSS extension -
    # ------------------------------------------------    
    
    print_routine("Plot extensions of output file (l2ps)")   
    print
    
    # use previously defined output filename from loreformat
    in_target_filename = out_target_filename  
    
    try:
        a_ps_execute(in_target_filename, "SPEC_NONSS", SPEC_NONSS_plot, "SPEC_NONSS", "green")
    except KeyError:
        pass
        
    if os.path.exists(SPEC_NONSS_plot):
        print_notification("Success.")
    else:
        print_notification("Failed.") 
        err.set_code(10, is_fatal=False)  
        
    # ------------------------------------------------
    # - GENERATE RASTER PLOT OF SPEC_NONSS extension -
    # ------------------------------------------------  
    
    print_routine("Plot extensions of output file (l2ps)")   
    print
    
    # use previously defined output filename from loreformat
    in_target_filename = out_target_filename   
        
    try:
        a_ps_execute(in_target_filename, "SPEC_SS", SPEC_SS_plot, "SPEC_SS", "blue")
    except KeyError:
        pass
        
    if os.path.exists(SPEC_SS_plot):
        print_notification("Success.")
    else:
        print_notification("Failed.") 
        err.set_code(11, is_fatal=False)   
        
    # -------------------------
    # - GENERATE MONTAGE PLOT -
    # -------------------------  
    
    print_routine("Montage extensions of output file (l2pi/l2ps)")       
    print
    
    # use previously defined output filename from loreformat
    in_target_filename = out_target_filename       
    
    fig = plt.figure()
    fig.suptitle("Raster image of L1_IMAGE and SPEC_* extensions for file " + in_target_filename, fontsize=12)
    fig.add_subplot(211)
    a_pi_execute(in_target_filename, "L1_IMAGE", "", "", save=False, hold=True)

    fig.add_subplot(212)
    try:
        a_ps_execute(in_target_filename, "SPEC_NONSS", "",  "", "green", leg_title="SPEC_NONSS", save=False, hold=True)
        a_ps_execute(in_target_filename, "SPEC_SS", "", "", "blue", legend=True, leg_title="SPEC_SS", save=False, hold=True) 
    except KeyError:
        pass
    plt.savefig(montage_plot, bbox_inches="tight")
    
    if os.path.exists(montage_plot):
        print_notification("Success.")
    else:
        print_notification("Failed.") 
        err.set_code(12, is_fatal=False)   

    # -------
    # - END -
    # -------       
        
    print_routine("Results")       
    print        
    fi_unix = time.time()
    exec_time = fi_unix - st_unix
    print_notification("Execution time: " + str(exec_time) + "s.")
    
    # make sure no negative routine error codes for essential routines
    rtn_keys    = []    
    rtn_codes   = []
    with open(cfg['general']['error_codes_file']) as f:
        for line in f:
            if line.startswith("L2"):
                rtn_keys.append(str(line.split('\t')[0]))
                rtn_codes.append(int(line.split('\t')[1]))

    for idx, i in enumerate(rtn_keys):
        if i == "L2STATF2" or i == "L2STATX1" or i=="L2STATX2":
            continue                                                    # skip 1D extraction related error codes.
        if rtn_codes[idx] < 0: 
            err.set_code(13)    

    err.set_code(0) 	# this is a bit of a bodge, it disregards the current error code!
	        
if __name__ == "__main__":
  
    print_header()
    
    parser = OptionParser()

    parser.add_option('--t', dest='f_target', action='store', default=L2_TEST_DIR + "/1H0323/v_e_20141115_14_1_0_1.fits", help="path to target file")
    parser.add_option('--r', dest='f_ref', action='store', default=L2_TEST_DIR + "/1H0323/v_e_20141115_14_1_0_1.fits", help="path to reference file")
    parser.add_option('--c', dest='f_cont', action='store', default=L2_TEST_DIR + "/1H0323/v_w_20141121_2_1_0_1.fits", help="path to continuum file")
    parser.add_option('--dir', dest='work_dir', action='store', default="test", help="path to working dir")
    parser.add_option('--rc', dest='ref_chk', action='store_true', help="perform reference frame check only")
    parser.add_option('--o', dest='clobber', action='store_true')
    (options, args) = parser.parse_args()

    f_target = options.f_target
    f_ref = options.f_ref
    f_cont = options.f_cont
    work_dir = options.work_dir
    ref_chk = options.ref_chk
    clobber = options.clobber
    
    # DEFINE EXTENSIONS    
    target_suffix       = "_target"
    ref_suffix          = "_ref"
    cont_suffix         = "_cont"
    arc_suffix          = "_arc"
    trim_suffix         = "_tr"
    cor_suffix          = "_cor"
    reb_suffix          = "_reb"
    ext_suffix          = "_ex"
    ss_suffix           = "_ss"    

    if ref_chk:
      chk_ref_run(f_ref, f_cont)
    else:
      full_run(f_target, f_ref, f_cont, work_dir, clobber)	        
    
