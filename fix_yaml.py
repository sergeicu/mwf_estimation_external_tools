#!/home/ch215616/miniconda2/envs/tch2/bin/python

"""This function corrects the echo ordering for .nrrd files which have been generated with DICOM to NRRD conversion tool. In certain cases, the aforementioned DICOM to NRRD tool incorrectly orders the echoes of the signal. This function corrects that mistake. 

Usage: 

    python fix_yaml.py --input $inputdir <path_to_nrrd_dir>
    python fix_yaml.py --input $inputdir --outdir <existing directory to save to>
    

Function was tested for Myelin Water Function acquisition files (CPMG sequence).
Works with single slice (2D) and volumetric (3D) data. 

Specify directory or file with .nrrd extension. Each .nrrd file must have a matching .yaml file.

The script will check each .yaml file for presence of multiple echoes. If the order is incorrect, it will fix the .nrrd file ordering. 

If no explicit output directory is specified in the arguments, a new file name will be created in the same directory as the original input file with a suffix to its name - "_echoes_reordered.nrrd"


"""

import os
import glob
import sys
import argparse
import shutil

import numpy as np
import nrrd     
       
def check_echo_ordering(file): 

    # check the the yaml file exists 
    assert os.path.exists(file), 'file does not exist'

    # fetch the yaml file 
    with open(file) as f:    
        lines = f.readlines()

    # find a starting line in yaml file where the list of echoes starts  
    for i,line in enumerate(lines): 
        if 'VOLUME_LIST' in line: 
            startingline = i+2       

    #get list of echo times as list of strings (read to the end of the file) 
    echoes_str = lines[startingline:]

    # extract echo times as integers from the list of strings
    #echo_list = [[int(s) for s in echo_str.split() if s.isdigit()] for echo_str in echoes_str]
    echo_list = [[float(s) for s in echo_str.split() if s.replace(".","").isdigit()] for echo_str in echoes_str]

    # proceed only if list is not empty (i.e. echo list exists) 
    if echo_list and echo_list[0]:
        multi_echo_file = True 
        original_echo_order = echo_list 

        # extract the lists from list (sometimes the yaml file output can cause the list = [[9],[18],...] insted of list = [9,18,...]
        if isinstance(original_echo_order[0],list): # if is list, correct  
            original_echo_order = [i[0] for i in original_echo_order]

        # check if the orders match 
        if original_echo_order != sorted(original_echo_order):     
            print(f"Incorrect ordering of echoes detected in YAML file for: {os.path.basename(file)}.")
            print(f"{original_echo_order}")
            print(f"Reordering echoes")
            proceed = True
        else:
            print(f"Echoes are in tact for {os.path.basename(file)}")

            if args.debug: 
                print(f"Echoes are in tact {original_echo_order}")
            proceed = False

    else: 
        multi_echo_file = False
        # do not proceed if no echo lists found 
        proceed = False
        original_echo_order = []

    return proceed, original_echo_order, multi_echo_file

def get_correcting_index(original_echo_order): 
    # sort the list of echoes into correct order, but save the index of the old (incorrect) order 
    correcting_index = [b[0] for b in sorted(enumerate(original_echo_order),key=lambda i:i[1])]

    # debug (check that the original echo order sorted by ascending order manually is the corrected echo order with sort function)
    corrected_echo_order_manual = [original_echo_order[i] for i in correcting_index]
    assert corrected_echo_order_manual == sorted(original_echo_order)

    return correcting_index 

def load_nrrd(file,original_echo_order):
    # load .nrrd file 
    nrrd_file = file.replace('.yaml','.nrrd')
    im, header = nrrd.read(nrrd_file)

    # check that the list of echoes matches the correct dimension of the image 
    if im.shape[0] != len(original_echo_order): 
        print('The first dimension of the NRRD file does not match the list of echoes')
        print(f"The nrrd dimensions are {im.shape}")
        print(f"Original echo list from YAML file: {original_echo_order}")
        sys.exit("Exiting. ")

    return nrrd_file, im, header

def sort_nrrd(im, correcting_index):

    if im.ndim == 4: 
        imnew = [im[i,:,:,:] for i in correcting_index]
    elif im.ndim == 3: 
        imnew = [im[i,:,:] for i in correcting_index]
    else: 
        sys.exit("NRRD file must be 3D or 4D. Please check the file dimensions.")

    imnew = np.array(imnew)    

    return imnew 

def load_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', type=str, required=True,help='directory to check for incorrect echoes. Must match nrrd and yaml filenames') 
    parser.add_argument('-o','--outdir', type=str, default=None, help='directory to write corrected filenames. If not specified, files written to the same directory with a suffix _echoes_reordered.nii.gz')
    parser.add_argument('--debug', action='store_true', help='print out the name of every filename being checked and list its echoes and if they are in tact') 
    
    args = parser.parse_args()

    return args 

    
def process(inputs, outdir):
    # check if file or folder 
    if os.path.isfile(inputs): 
        # check if extension is .nrrd 
        assert inputs.endswith('.nrrd') or inputs.endswith('.yaml'), "Please specify path to .nrrd or .yaml file, or to a directory"
        
        # must specify input as .yaml file
        if inputs.endswith('.nrrd'):             
            
            path, ext = os.path.splitext(inputs)
            yaml_file = path+'.yaml'
            # check if corresponding .yaml file exists 
            assert os.path.exists(yaml_file), f"Matching yaml file does not exist: {yaml_file}"
        elif inputs.endswith('.yaml'):
            yaml_file = inputs
            # check if corresponding .yaml file exists 
            assert os.path.exists(yaml_file.replace('.yaml','.nrrd')), f"Matching .nrrd file does not exist: {yaml_file}"
            
        
        # must provide input as list 
        files = [yaml_file]
    else: 
# os.path.isdir(args.input):
        assert os.path.exists(inputs), "Specified directory does not exist. Please check path"
        
        # find all .yaml files 
        inputs = inputs + "/" if not inputs.endswith("/") else inputs
        files = glob.glob(inputs + "*.yaml")
        
        assert files, f"No .yaml files found in the specified directory:\n{inputs}"
        
        # check if corresponding .nrrd exist 
        for file in files: 
            path, ext = os.path.splitext(file)
            assert os.path.exists(path+'.nrrd'), f"Matching .nrrd file to .yaml file does not exist {path+'.nrrd'}"
           
    
    newnames = []
   
    for file in files: 
            
        # check the yaml file if the order of echoes is correct
        proceed, original_echo_order,multi_echo_file = check_echo_ordering(file)

        # placeholder for a file path for a newly converted file
        newname=None
        
        # proceed with correcting the file 
        if proceed: 
            # get the index with which to re-order the echoes 
            correcting_index = get_correcting_index(original_echo_order)

            nrrd_file, im, header = load_nrrd(file,original_echo_order)

            imnew = sort_nrrd(im,correcting_index)


            # if the output directory is different, change the savedir  
            if outdir is not None:
                outdir = outdir + "/" if not outdir.endswith("/") else outdir
                assert os.path.exists(outdir), "Specified output directory does not exist. Please specify full correct path to output directory."
                assert os.path.dirname(outdir) != os.path.dirname(nrrd_file), f"if explicit output path is specified, it has to be different from the path to input files. Else, do not specify output path and file will be written to same directory"
                newname = outdir+os.path.basename(nrrd_file)
            else: 
                # rename the file (only executed if explicit outdir is specified)
                newname = nrrd_file.replace('.nrrd','_echoes_reordered.nrrd')
            
            # check if file already exists: 
            if os.path.exists(newname): 
                base = os.path.basename(newname)
                print(f"WARNING: A re-ordered .nrrd for this specific file already exists.\n{base}\nProceed? y/n")
                proceed_to_save = input()
                if proceed_to_save.lower() == 'n':
                    continue
                

            # save the new file
            nrrd.write(newname,imnew,header=header)
            print(f"Saved the new file to: {newname}\n")
        else: 
            # if args.outdir is explicitly specified and the file is a multi echo file - copy them to folder that folder
            # i.e. copy file with correct and incorrect echo ordering 
            if multi_echo_file and outdir is not None: 
                # copy nrrd file 
                shutil.copyfile(file, outdir+"/" + os.path.basename(file))
                # copy yaml file 
                file = file.replace(".yaml", ".nrrd")
                shutil.copyfile(file, outdir+"/" + os.path.basename(file))
        newnames.append(newname)

    return newnames

    
if __name__ == '__main__': 

    args = load_args()

    _ = process(args.input, args.outdir)