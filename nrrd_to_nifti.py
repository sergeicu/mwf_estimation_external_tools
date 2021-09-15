import os 
import argparse
import glob 
import sys 

import nrrd 
import nibabel as nb 
import numpy as np

def load_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', type=str, required=True,help='directory or file') 
    parser.add_argument('--cpmg', type=str, choices=['1','0'], default='1', help='process CPMG - will check for orientation of echoes')
    args = parser.parse_args()

    return args 

def convertCPMG(im):
    # reorder directions 
    dims = im.shape
    if dims[0] < 33: # assuming that this is MWF data, the 32 echoes will be put to the front. Any other kind of MR acquisition would have the first dimension much larger than 33 (since smallest matrix is typically 64x64)
        print(f"Changing order of axes to: ")
        if len(dims) == 4: 
            im = np.transpose(im,(1,2,3,0))    
        elif len(dims) == 3: 
            im = np.transpose(im,(1,2,0))    
        print(im.shape)
    
    return im

if __name__ == '__main__': 

    args = load_args()
    
    # check if file or folder 
    if os.path.isfile(args.input): 
        # check if extension is .nrrd 
        assert args.input.endswith('.nrrd'), "Please specify path to .nrrd file, or to a directory"    
        
        # save as list for processing
        files = [args.input]
    else: 
        assert os.path.exists(args.input), "Specified directory does not exist. Please check path"
        
        args.input = args.input + "/" if not args.input.endswith("/") else args.input
        files = glob.glob(args.input + "*.nrrd")
        
        assert files, f"No .nrrd files found in the specified directory:\n{args.input}"
        
    
    for file in files: 
        print(f"{file}")
        # load nrrd and save as .nii.gz
        im, _ = nrrd.read(file)
        
        if args.cpmg == '1':
            im = convertCPMG(im)
        
        img = nb.Nifti1Image(im, np.eye(4))
        nii_file = file.replace('.nrrd','.nii.gz')
        nb.save(img, nii_file)
        print(f"Saved {nii_file}")
    print('Done.')