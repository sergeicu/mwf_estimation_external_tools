# Warning: all PHI data was substituted with <REMOVED> flag 

# Patient details
MRN:<REMOVED>
Date:<REMOVED>

# Get code 
git clone git@github.com:sergeicu/mwf_estimation_external_tools.git
cd mwf_estimation_external_tools/
codedir=$(pwd)

# Initialize correct conda environment 
# (prereqs: see email sent to Onur and Sila on 20210914 with title: "To use my conda (and newer python)")
# this loads the correct version of python and code dependencies
conda activate tch2_yml 

# Create new directory 
rootdir=/fileserver/fastscratch/serge/s20210815_mwf_example_processing_for_onur  # root folder in which all the processed data will be stored 
experiment=<REMOVED> # subfolder that identifies current experiment (e.g. today's date)
mkdir $rootdir/$experiment
cd $rootdir/$experiment
cp -r $codedir/* . # copy all our code here

# GET DICOMS 
MRN=<REMOVED>
DATE=<REMOVED> 
savedir=$rootdir/$experiment/
bash retrieve2.sh $MRN $DATE $savedir

# DICOM2NRRD
session=<REMOVED> # a child folder that retrieve.sh script creates when it fetches dicoms. This corresponds to MRI session id 
scan=7_myelin_1x1x4 # name of a particular scan that we are interested in
scandir=$rootdir/$experiment/$MRN/$session//$scan
python dicom_to_nrrd.py --dicomdir $scandir

# FIX YAML 
nrrd=$scandir/NRRD/
python fix_yaml.py --input $nrrd

# NRRD2NIFTI **
nrrd_fixed=$scandir/NRRD/
python nrrd_to_nifti.py --input $nrrd_fixed

# (if not installed) install julia language on your machine (via conda, or yum install), then: 
# add DECAES.js REPL to your your julia language packages: 1. type "julia" in terminal 2. type "]" inside julia to enter REPL 3. type "add DECAES" to download package

# run julia
julia_out=$scandir/NRRD/out_julia
mkdir $julia_out
cp -r $scandir/NRRD/*.nii.gz $julia_out/
python run_julia.py --input $julia_out/*.nii.gz

# (if not downloaded) download anima docker image (prereqs: installed docker software, and sudo access to the machine)
name=sergeicu/anima:latest
sudo docker pull $name

# run anima 
anima_out=$scandir/NRRD/out_anima
mkdir $anima_out
cp -r $scandir/NRRD/*.nii.gz $anima_out/
localfolder=$anima_out # will grab all .nii.gz images
sudo docker run -it --rm -v $localfolder:/data sergeicu/anima:latest /bin/bash
python run_anima.py /data/
exit

# ** = do not use direct `dicom to nifti` converter tool because the header is not preserved correctly. 
