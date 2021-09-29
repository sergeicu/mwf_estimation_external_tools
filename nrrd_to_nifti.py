"""Updated version. Allows preservation of the header information"""


import os 
import argparse
import glob 
import sys 
import copy 

import nrrd 
import nibabel as nb 
import numpy as np

import svtools as sv

def load_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', type=str, required=True,help='directory or file') 
    parser.add_argument('-r','--remove_header', action="store_true", help='does not preserve the header of the file. Akin to older version of the function') 
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

def convertCPMG2(f,save_intermediate=False):

    im, hdr = nrrd.read(f)
    hdr2=copy.copy(hdr)

    assert im.shape[0] == 32
    hdr2['sizes'] = hdr2['sizes'][1:]
    hdr2['kinds'] = hdr['kinds'][1:]
    hdr2['space directions'] = hdr['space directions'][1:]
    hdr2['dimension'] = 3

    basename = os.path.basename(f).replace(".nrrd", "")
    savedir = os.path.dirname(f) + "/" + basename + "/"
    os.makedirs(savedir, exist_ok=True)
    names = []
    
    print(f"Splitting echoes")
    for i in range(0,32):
        print(f"{i}/32")
        echo_i = im[i,:,:,:] if im.shape[0]==32 else im[:,:,:,i]
        adj = "0" if i<10 else ""
        savename = savedir + basename + adj + str(i) + ".nrrd"
        nrrd.write(savename,data=echo_i,header=hdr2)
        names.append(savename)

    # now let's convert all to .nii.gz
    print(f"Converting to .nii.gz")
    for i,n in enumerate(names): 
        print(f"{i}/32")
        sv.crl_convert_format(n, ".nii.gz",verbose=False)

    # now let's merge all .nii.gz 
    print(f"Merging .nii.gz echoes")
    final_output = savedir + basename + ".nii.gz"
    cmd = ["fslmerge", "-t", final_output]
    names_nii = [n.replace(".nrrd", ".nii.gz") for n in names]
    cmd.extend(names_nii)
    sv.execute(cmd)

    if not save_intermediate:
        for n in names:
            sv.execute(["rm", "-rf", n])
        for n in names_nii:
            sv.execute(["rm", "-rf", n])                    


    return final_output









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
        
        if not args.remove_header:
            
            # convert the .nrrd into .nii.gz (and flip dimensions from [echoes, x,y,z] to [x,y,z,echoes]) WHILE preserving the heade
            nii_file = convertCPMG2(file)
            
        else:
            # load nrrd and save as .nii.gz
            im, _ = nrrd.read(file)

            if args.cpmg == '1':
                im = convertCPMG(im)

            img = nb.Nifti1Image(im, np.eye(4))
            nii_file = file.replace('.nrrd','.nii.gz')
            nb.save(img, nii_file)
    print(f"Saved {nii_file}")
    print('Done.')