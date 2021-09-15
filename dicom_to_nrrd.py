""" Convert DICOM to NRRD for CPMG files. 

    Usage: 
            python dicom_2_nrrd.py --dicomdir <fullpath to dicom directory> --outdir <fullpath>
            
    The outputs are written to 'NRRD' folder for each dicomdir. 
            
            

"""

import os 
import shutil 
import glob 
import argparse 

import svtools as sv

def load_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--dicomdir', nargs='+', type=str, required=True,help='one or more DICOM directory paths') 
    args = parser.parse_args()

    return args 

def convert_DICOMS(paths):
        
    # run convert mechanism     
    for i, newpath in enumerate(paths): 
        
        print(f"\n\n\n{i}/{len(paths)}: {newpath}")

        outputdir = newpath + "/NRRD/"
        os.makedirs(outputdir,exist_ok=True)
        
        # convert 
        cmd = ["/opt/el7/pkgs/crlDcm/crlDcmConvert", "--volumeStack", "3D-MC", newpath, outputdir]
        sv.execute(cmd)
    
    print("Conversion complete.")

    
if __name__ == '__main__':
    
    # get args 
    args = load_args()
    
    # wrap into a list if single directory 
    if isinstance(args.dicomdir, str):
        args.dicomdir = [args.dicomdir]
    
    # process
    convert_DICOMS(args.dicomdir)
    
    
    
    