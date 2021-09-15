"""

Helper tools for loading, converting and plotting MR images

Author: 
    serge.vasylechko@tch.harvard.edu 

Date: 
    2020-07-23

"""

import subprocess 
import os 
import json 
import sys
import glob 
import copy 
import shutil
import argparse
import re
import fnmatch
import math 
from inspect import currentframe, getframeinfo

import nrrd
import pickle
import git
import nibabel as nb
import numpy as np 
import matplotlib.pyplot as plt 
import yaml 
from yaml2object import YAMLObject



# -----------
# NEW - please add and comment
# -----------


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)



# -----------
# Execute in bash 
# -----------


def execute(cmd,sudo=False, remote=False):
    """Execute commands in bash and print output to stdout directly"""

    if sudo:
        cmd = ["sudo"]+cmd
    if remote:
        print(" ".join(cmd))
    else:
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='') # process line here

        if p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, p.args)

# -----------
# File conversion 
# -----------

def crl_convert_format(image, format_out,dirout=None, verbose = True,debug=False): 
    """
    Python wrapper for crlConvertBetweenFileFormats that converts between .nrrd, .nii.gz, .vtk 

    Args: 
        image (str):      path to image 
        format_out (str): enum to .nrrd, .vtk, .nii.gz, .nii
        [dirout] (str):   convert into a specific directory 
        [verbose] (str):  print the command being executed 
    Returns: 
        str: Path to converted file 

    """    
    """py wrapper for crlConvertBetweenFileFormats tool"""
    
    #print("WARNING: crlConvertBetweenFileFormats is NOT currently working when invoked from the terminal")
    #crlConvertFormat="/opt/el7/pkgs/crkit/nightly/20160503/bin/crlConvertBetweenFileFormats"
    crlConvertFormat="/opt/el7/pkgs/crkit/2021/crkit-master/bin/crlConvertBetweenFileFormats"
    #    crlConvertFormat="/opt/x86_64/pkgs/crkit/nightly/20170107/crkit/bin/crlConvertBetweenFileFormats"
    #    crlConvertFormat="/opt/x86_64/pkgs/crkit/march-2009/bin/crlConvertBetweenFileFormats" (if necessary)

    # explicitly search for the format_in (to avoid errors like before)
    if image.endswith('.nii.gz'):
        format_in = '.nii.gz'
    elif image.endswith('.nii'):
        format_in = '.nii'
    elif image.endswith('.nrrd'):
        format_in = '.nrrd'
    elif image.endswith('.vtk'):
        format_in = '.vtk'
    else: 
        print("FILE FORMAT IS WRONG")
    if dirout: 
        dirout = dirout+'/' if not dirout.endswith('/') else dirout
        d,f = os.path.split(image)
        imageout=dirout+f
    else:
        imageout=image
    #cmd = [crlConvertFormat, "-in", image, "-out", imageout.replace(format_in, format_out)]
    cmd = [crlConvertFormat, image, imageout.replace(format_in, format_out)]
    if debug: 
        print(" ".join(cmd))
    if verbose:
        print(f"converting: {image} from {format_in} to {format_out}")
    #subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    execute(cmd)
    return imageout.replace(format_in, format_out)

def vtk2nrrd(path): 
    """
    vtk2nrrd and nrrd2vtk functions convert between file formats.

    Functions can handle multiple input types. 
    If input is:
        - directory - convert all .vtk/.nrrd files in a directory 
        - list of file paths - convert all files in a list 
        - file path - convert a single file 

    Args: 
        path (str): directory, file path or list of file paths  

    """    
    if os.path.isdir(path): 
        vtks = glob.glob(path+".vtk")
        for vtk in vtks: 
            crl_convert_format(vtk, '.nrrd')
    elif os.path.isfile(path): 
        crl_convert_format(path, '.nrrd', verbose=True)
    elif path is list: 
        for vtk in path: 
            crl_convert_format(vtk, '.nrrd')
    else: 
        sys.exit("Incorrect input to vtk2nrrd() ")
    
def nrrd2vtk(path):
    """
    vtk2nrrd and nrrd2vtk functions convert between file formats.

    Functions can handle multiple input types. 
    If input is:
        - directory - convert all .vtk/.nrrd files in a directory 
        - list of file paths - convert all files in a list 
        - file path - convert a single file 

    Args: 
        path (str): directory, file path or list of file paths 
 

    """    
    if os.path.isdir(path): 
        vtks = glob.glob(path+".nrrd")
        for vtk in vtks: 
            crl_convert_format(vtk, '.vtk')
    elif os.path.isfile(path): 
        crl_convert_format(path, '.vtk', verbose=True)
    elif path is list: 
        for vtk in path: 
            crl_convert_format(vtk, '.vtk')
    else: 
        sys.exit("Incorrect input to nrrd2vtk() ")    

# -----------
# JSON / PICKLE
# -----------

def read_from_json(filename):
    """
    Read JSON file 

    Args: 
        filename (str): path to file 
    Returns: 
        data (dict): json file 

    """    
    with open(filename,'r') as f: 
        data = json.load(f)
    return data        

def write_to_json(dictionary,filename):
    """
    Write JSON file 

    Args: 
        filename (str): path to file 
    Returns: 
        data (dict): json file 

    """    
    with open(filename,'w') as f: 
        json.dump(dictionary,f)

def pickle_dump(filename, pickle_object):
    # wrapper around the pickle dump method (to avoid the clutter) 
    pickle.dump(pickle_object,open(filename,'wb'))
    
def pickle_load(filename):
    # wrapper around the pickle load method (to avoid the clutter) 
    return pickle.load(open(filename,'rb'))


def read_yaml_as_object(path):
    
    assert os.path.exists(path)
    assert path.endswith('.yaml')
    obj = YAMLObject('DefaultConfig', (object,), {'source': path, 'namespace': 'defaults'})

    return obj

def read_yaml_as_dict(path):

    assert os.path.exists(path)
    assert path.endswith('.yaml')
    #opt = YAMLObject('DefaultConfig', (object,), {'source': path, 'namespace': 'defaults'})
    with open(path, 'r') as stream:
        try:
            dictionary = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return dictionary

# -----------
# Strings 
# -----------


def get_nums_from_string(string):
    """
    Get all numbers from a string. 
    
    Args: 
        string (str): string to be formatted 
    Returns: 
        rr (list): list all numbers 
    
    Example: 
        # Format is [(<string>, <expected output>), ...]
        ss = [("apple-12.34 ba33na fanc-14.23e-2yapple+45e5+67.56E+3",
               ['-12.34', '33', '-14.23e-2', '+45e5', '+67.56E+3']),
              ('hello X42 I\'m a Y-32.35 string Z30',
               ['42', '-32.35', '30']),
              ('he33llo 42 I\'m a 32 string -30', 
               ['33', '42', '32', '-30']),
              ('h3110 23 cat 444.4 rabbit 11 2 dog', 
               ['3110', '23', '444.4', '11', '2']),
              ('hello 12 hi 89', 
               ['12', '89']),
              ('4', 
               ['4']),
              ('I like 74,600 commas not,500', 
               ['74,600', '500']),
              ('I like bad math 1+2=.001', 
               ['1', '+2', '.001'])]

        for s, r in ss:
            rr = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", s)
            if rr == r:
                print('GOOD')
            else:
                print('WRONG', rr, 'should be', r)    
    
    """
    return re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", string)

def strip_trailing_whitespace(filename):
    """\
    strip trailing whitespace from file
    usage: strip_trailing_whitespace.py <file>
    """

    content = ''
    outsize = 0
    inp = outp = filename
    with open(inp, 'rb') as infile:
        content = infile.read()
    with open(outp, 'wb') as output:
        for line in content.splitlines():
            newline = line.rstrip(" \t")
            outsize += len(newline) + 1
            output.write(newline + '\n')

    print("Done. Stripped %s bytes." % (len(content)-outsize))


# -----------
# Argparse 
# -----------

def init_argparse(): 
    from argparse import Namespace
    args = Namespace()
    return args 

def print_args(args):
    for arg in vars(args):
        print(arg, getattr(args, arg))
    

# -----------
# Search files
# -----------


def find_files(pattern, path):
    """
    Find all files that match a naming pattern. Search subdirectories.
    
    Args:
        pattern (str): pattern to match with a `*` 
        path (str):    directory to search in. Will search in subdirectories also. 
        
    Returns: 
        result (list): list of files that match the pattern
    
    Example: 
        find('*.txt', '/path/to/dir')
    """
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


# -----------
# Dictionaries
# -----------


def vars_to_dict(varnamelist):
    for i in varnamelist:
        d[i] = locals()[i]
    return d 

# -----------
# Plot images (python)
# -----------

def plot_slice(v,figtitle,slice_ = None, figsize=(5,5)):
    # plot a single image 
    
    # input `v` should be a numpy array with the order of dims of (x,y,z)
 
    if v.ndim != 3 and v.ndim != 2: 
        print(f" Shape of input is {v.shape}") 
        print(f" Number of dims is {v.ndim}")         
        sys.exit('Make sure that number of dims is 2 or 3')

    if v.ndim == 2: 
        # repmat the slice 3 times so it is plottable 
        v = np.repeat(v[:, :, np.newaxis], 3, axis=2)        
        
    if slice_ is None: 
        slice_ = v.shape[-1]//2 
        print(f"Slice {slice_}")

    ax = plt.figure(figsize=figsize)
    plt.imshow(v[:,:,slice_],cmap='gray')
    plt.suptitle(figtitle, fontsize=16,y=1.01)
    plt.axis('off')
    plt.tight_layout()
    
def plot_slices(v,figtitle,slices = None, figsize=(15,15),columns=3):
    # plot a single image 
    
    # input `v` should be a numpy array with the order of dims of (x,y,z)
    
    # slices can be either 
        # [bottom_slice, top_slice] 
        # list of slices 
        # none (default is to plot 10 slices around the middle)
 
    if v.ndim != 3: 
        print(f" Shape of input is {v.shape}") 
        print(f" Number of dims is {v.ndim}")         
        sys.exit('Make sure that number of image dims is 3')

    if slices is None: 
        slices = [v.shape[-1]//2-5, v.shape[-1]//2 +5]
        slices = list(range(slices[0],slices[1]))
    elif len(slices)==2: 
        if slices[1]-slices[0]>=2:
            slices = list(range(slices[0],slices[1]))
        elif slices[1]-slices[0]<1:
            sys.exit("Error. Please check slices you are plotting")        
    elif len(slices)<2: 
        sys.exit("Error. Please plot at least 6 slices")

    if columns <2: 
        sys.exit("n of columns must be greater than 1")   
    
    # check if specified slice exceed actual number of slices in the image 
    if any([i>v.shape[-1] for i in slices]): 
        sys.exit("Please specify correct slices. There are only %s slices in image" % (v.shape[-1]))


    L = len(slices)

    cols = columns
    #    rows = math.ceil(L/cols) if math.ceil(L/cols)>1 else 2 
    rows = L//cols
    fig, axs = plt.subplots(rows, cols,figsize=figsize)    
    fig.suptitle(figtitle, fontsize=16,y=1.01)
    k = 0 
    if rows != 1: 
        for i in range(rows): 
            for j in range(cols): 
                # if statement that will 
                if k<len(slices):
                    im = v[:,:,slices[k]]
                    slicenum = 'slice'+str(slices[k])
                else: 
                    # plot empty images
                    im = np.zeros_like(v[:,:,0])            
                    slicenum = ' '          
                axs[i,j].imshow(im,cmap='gray')
                axs[i, j].title.set_text(slicenum)
                axs[i, j].axis('off')
                k=k+1
    else: 
        for j in range(cols): 
            # if statement that will 
            if k<len(slices):
                im = v[:,:,slices[k]]
                slicenum = 'slice'+str(slices[k])
            else: 
                # plot empty images
                im = np.zeros_like(v[:,:,0])            
                slicenum = ' '          
            axs[j].imshow(im,cmap='gray')
            axs[j].title.set_text(slicenum)
            axs[j].axis('off')
            k=k+1            
    plt.tight_layout()  
 
def plot_ims(figtitle, root,compare_dirs,image_types,suffix='',slice_=None,figsize=(10,10),adjust_title=0,compare_dir_names=None,image_type_names=None,adjust_contrast=False,adjust_scale=2,mask=None):
    
    
    # compare images from different directories by supplying the directories and image types to plot 

    # choose which images to plot
    image_type_defaults = {'S0':'meanB0_1',
                           'D':'meanADC_1',
                           'P':'meanPER_1',
                           'Pf':'meanPER_FRAC_1',
                           'PF':'meanPER_FRAC_1',
                           'b0':'b0_averaged',
                           'b50':'b50_averaged',
                           'b100':'b100_averaged',
                           'b200':'b200_averaged',
                           'b400':'b400_averaged',
                           'b600':'b600_averaged',
                           'b800':'b800_averaged'}
    
    image_names = []
    for i in image_types: 
        if i in image_type_defaults.keys(): 
            image_names.append(image_type_defaults[i]+'.nrrd')
        else: 
            image_names.append(i+'.nrrd')
    

    # custom name for the compare dirs 
    if compare_dir_names is None:
        compare_dir_names = compare_dirs

    if image_type_names is None:
        image_type_names = image_types         
        
    ims_loaded = []
    names = []
    im_contrast = np.zeros((len(compare_dirs),len(image_names),2))  # for looking up contrast 
    # get all images 
    for i, compare_dir in enumerate(compare_dirs):
        print(compare_dir)
        for j,file_name in enumerate(image_names): 
            f = root + "/" + compare_dir + "/" + suffix + "/" + file_name
            im_loaded,_ = nrrd.read(f)
            if f.endswith('meanPER_FRAC_1.nrrd'):
                # remove zero values from roar 
                im_loaded[im_loaded<0] = 0 
            
            # adjust contrast 
            im_mean = np.mean(im_loaded[im_loaded!=0])
            im_std = np.std(im_loaded[im_loaded!=0])    
            if adjust_contrast:
                vmin = 0
                vmax = im_mean + adjust_scale*im_std
            else: 
                vmin = None
                vmax = None
            im_contrast[i,j] = [vmin,vmax]
            
            ims_loaded.append(im_loaded)
            names.append(compare_dir_names[i] + "\n" + image_type_names[j])            
            
    if slice_ is None:
        slice_ = ims_loaded[0].shape[-1]//2                    
        
    columns = len(image_types)    
    if columns <2: 
        sys.exit("n of columns must be greater than 1")          
        
    # plot all images 
    L = len(ims_loaded)
    cols = columns
    rows = L//cols
    
    
    fig, axs = plt.subplots(rows, cols,figsize=figsize)  
    fig.suptitle(figtitle, fontsize=16,y=adjust_title)    
    k = 0 
    if rows != 1: 
        for i in range(rows): 
            for j in range(cols): 
                vmax = im_contrast[i,j,1] if not np.isnan(im_contrast[i,j,1]) else None                   
                axs[i,j].imshow(ims_loaded[k][:,:,slice_],cmap='gray',vmin=0, vmax=vmax)
                axs[i, j].title.set_text(names[k])
                axs[i, j].axis('off')
                k=k+1
    else: 
        for j in range(cols): 
            vmax = im_contrast[i,j,1] if not np.isnan(im_contrast[i,j,1]) else None                        
            axs[j].imshow(ims_loaded[k][:,:,slice_],cmap='gray',vmin=0, vmax=vmax)
            axs[j].title.set_text(names[k])
            axs[j].axis('off')
            k=k+1
            
    plt.tight_layout()     

def plot_bvals(v,figtitle,slice_ = None, figsize=(10,10),columns=3,adjust_title=0):
    # plot first 6 bvalues from a list of b-value images 
    
    # input `v` should be: 
        # a list of images
        # or 
        # a numpy array with the order of dims of (bvals, x,y,z)

    # get `v` via svtools.get_dwi_images(path,'_averaged')
    
    if isinstance(v,list): 
        v = np.array(v)
    elif isinstance(v,np.ndarray): 
        if v.shape[0] != 7: 
            if v.shape[-1] == 7: 
                print('Warning: Transposing dims (x,y,z,bvals) -> (bvals,x,y,z)')
                v = np.transpose(v,(-1,0,1,2))
            else:
                print(f" Shape of input is {v.shape}") 
                sys.exit('Make sure that first dimension is the b-value order')

    if slice_ is None: 
        slice_ = v.shape[-1]//2 
        
    L = v.shape[0]
    cols = columns
    rows = L//cols 
    
    if columns <2: 
        sys.exit("n of columns must be greater than 1")       
    
    fig, axs = plt.subplots(rows, cols,figsize=figsize)    
    fig.suptitle(figtitle, fontsize=16,y=adjust_title)
    k = 0 
    bvals = [0,50,100,200,400,600,800]
    if rows != 1: 
        for i in range(rows): 
            for j in range(cols): 
                axs[i,j].imshow(v[k,:,:,slice_],cmap='gray')
                axs[i, j].title.set_text('B'+str(bvals[k]))
                axs[i, j].axis('off')
                k=k+1
    else: 
        for j in range(cols): 
            axs[j].imshow(v[k,:,:,slice_],cmap='gray')
            axs[j].title.set_text('B'+str(bvals[k]))
            axs[j].axis('off')
            k=k+1        
    plt.tight_layout()
   
def plot_params(v,figtitle,slice_ = None, figsize=(10,10)):
    # plot 4 IVIM parameters 
    
    # input `v` should be a dictionary that contains 'S0','D','P','P_f' entries

    # get `v` via svtools.get_ivim_images(path)
    
    if slice_ is None: 
        slice_ = v['S0'].shape[-1]//2 
        
    # remove all values below zero for P_f 
    v['P_f'][v['P_f']<0]=0
    
    L = 4
    cols = 2
    rows = L//cols 
    fig, axs = plt.subplots(rows, cols,figsize=figsize)    
    fig.suptitle(figtitle+' : slice'+str(slice_), fontsize=16,y=1.01)
    axs[0,0].imshow(v['S0'][:,:,slice_],cmap='gray')
    axs[0,1].imshow(v['D'][:,:,slice_],cmap='gray')
    axs[1,0].imshow(v['P'][:,:,slice_],cmap='gray')
    axs[1,1].imshow(v['P_f'][:,:,slice_],cmap='gray')
    axs[0, 0].title.set_text('S0')
    axs[0, 1].title.set_text('D')
    axs[1, 0].title.set_text('P')
    axs[1, 1].title.set_text('P_f')    
    axs[0, 0].axis('off')
    axs[0, 1].axis('off')    
    axs[1, 0].axis('off')
    axs[1, 1].axis('off')
    plt.tight_layout()

# -----------
# Saving input - source/args/git
# -----------

def save_source(thisfile_handle, savedir,justthisfile=False): 
    
    """
    Save copy of this script into a specified directory. Save args (if specified).

    Args: 
        thisfile_handle (sys.argv[0] object) - necessary to be invoked from the file that's running it 
        savedir (str): path to dir
        args: variable to save - [str, list, npy, dict, numpy.ndarray, arparse.Namespace]
        args_savetype: 
    """        
    savedir = savedir + "/" + 'source/' if not savedir.endswith('/') else savedir + 'source/'
    os.makedirs(savedir, exist_ok=True)

    if justthisfile:
        # save source 
        this_script = sys.argv[0]
        base,name = os.path.split(this_script)
        shutil.copy(this_script,savedir+name)
    else:
        thisdir = os.path.dirname(os.path.abspath(__file__))
        files = glob.glob(thisdir+'/*.py')
        _ = [shutil.copyfile(i,savedir+os.path.basename(i)) for i in files]

    print(f"Saved source to {savedir}")    

def save_args(savedir, args=None,args_savetype=None):
    """
    Save args to specified directory.

    Args: 
        savedir (str): path to dir
        args: variable to save - [dict, sys.argv, arparse.Namespace]
        args_savetype: type type to save - [npy,pkl,txt]. If not specified, isinstance(args,dict) or isinstance(args,argparse.Name) -> saved to json
    """     
    if not os.path.exists(savedir):
        os.makedirs(savedir)
        print('Created directory')
        print(savedir)

    savedir = savedir + "/" if not savedir.endswith('/') else savedir 
    print('Args will be saved to:')
    print(savedir)
    saved = True

    if isinstance(args,dict):
        # save to json 
        if args_savetype != 'npy' and args_savetype != 'pkl':
            # check every dict value for numpy and save to list instead (else there is error)
            for k in args: 
                if isinstance(args[k], np.ndarray):
                    args[k].tolist()
            write_to_json(args,savedir+"args.json") # save args into .json    
        # save to pkl
        elif args_savetype == 'pkl':
            pickle_dump(savedir+"args.pkl", args) # save args into .pkl
        # save to npy 
        elif args_savetype == 'npy': 
            for k in args: 
                np.save(savedir+k+'.npy', args[k]) 
    elif isinstance(args,argparse.Namespace):
        # save argparse to .json 
        write_to_json(vars(args),savedir+"args.json") # save args into .json        
    elif args=='argv':
        # save argv to .txt
        with open("args.txt",'w') as f: # save argv into .txt 
            f.write(' '.join(sys.argv))
    else:
        saved = False
        print('WARNING: args were NOT saved.')
        print('args types to save can only be a dict, argparse.Namespace object or `argv`')

    if saved:
        print(f"Saved args")    
  
def save_git_status(savedir):
    """
    Save sha hash of current git head for the repository and git status

    Args: 
        savedir (str): path to dir
    """     
    if not os.path.exists(savedir):
        os.makedirs(savedir)
        print('Created directory')
        print(savedir)

    savedir = savedir + "/" if not savedir.endswith('/') else savedir 
    print('Git hash and status will be saved to:')
    print(savedir)
    saved = True
    
    # prep git hash 
    dictionary = {'sha':'','status':''}
    repo = git.Repo(search_parent_directories=True)
    dictionary['sha'] = repo.head.object.hexsha
    # prep git status 
    dictionary['status'] = repo.git.status()
    
    # save 
    write_to_json(dictionary,savedir+"git_status.json") # save args into .json  
    
    if saved:
        print(f"Saved git head sha and status")    

def is_git_repo(path):
    try:
        _ = git.Repo(path, search_parent_directories=True).git_dir
        return True
    except git.exc.InvalidGitRepositoryError:
        return False
        
        
# -----------
# Image viewers (bash)
# -----------

def nrrd_temp(im,suffix="",savedir=None,header=None):
    """Save an array to .nrrd file and view output quickly"""
    
    # check correct slices 
    assert im.ndim == 3 or im.ndim == 2, "image must be 3D or 2D"
    
    im = np.nan_to_num(im)
    
    if savedir is not None: 
        assert os.path.exists(savedir), "directory does not exist"
        d = savedir + '/' if not savedir.endswith("/") else savedir
    else: 
        d = os.getcwd() + '/' 
    
    savename = "TEST_REMOVE_" + suffix + ".nrrd"
    if header is not None: 
        nrrd.write(d+savename,im,header=header)
    else: 
        nrrd.write(d+savename,im)
    
    print(f"Saved to: {d+savename}")
    print(f"cd {d}")
    print(f"itksnap -g {savename}")
    print(f"mv {d+savename} $trash")        

def temp_itksnap(array, ref_nii=None, compare_niis=None,custom_name=None):
    
    # fetch ref_nii if available 
    if ref_nii is not None:
        assert os.path.exists(ref_nii)
        ref_imo = nb.load(ref_nii)
    elif compare_niis is not None:
        # fetch ref from first nii 
        assert isinstance(compare_niis, list), "compare_nii must be a list of paths "
        assert [os.path.exists(i) for i in compare_niis]
        ref_imo = nb.load(compare_niis[0])
    else:
        ref_imo = None

    # save array to temp directory (which will then get removed)
    if ref_imo is not None:
        imnewo = nb.Nifti1Image(array,affine=ref_imo.affine, header=ref_imo.header)
    else:
        imnewo = nb.Nifti1Image(array,affine=np.eye(4))
    if custom_name is not None:
        savename = '/home/ch215616/trash/' + 'temp_itksnap'+'_'+custom_name+'.nii.gz'
    else:
        savename = '/home/ch215616/trash/' + 'temp_itksnap.nii.gz'
    nb.save(imnewo,savename)
    
    print(f"Saved temp image to:{savename}")


    files = [savename]
    if compare_niis is not None:
        files.extend(compare_niis)
    
    itksnap(files)


def itksnap(ims, seg=None, remote=False):
    """plot list of images in itksnap 

    Args: 
        ims (list): list of paths to file 
        seg (str): (optional) path to segmentation file
        remote (bool): (optional) if True, prints the path to executing itksnap on rayan, for use via Terminal in VNCserver 

    Returns: 
        opens itksnap window 

    """

    def _check_inputs(ims,seg=None):
        # verify that inputs are correct

        # check that ims are supplied as string and not numpy arrays 
        assert (isinstance(ims,list) or isinstance(ims, str)), f"'Ims' must be a list of strings or a string"

        # if string provided, turn into list
        if isinstance(ims,str):
            ims = [ims]

        # check that every path exists 
        for im_path in ims:
            assert os.path.exists(im_path), f"Path does not exist:\n{im_path}"
        # check segmentation 
        if seg is not None: 
            assert os.path.exists(seg), f"Segmentation does not exist:\n{im_path}"

        return ims # return modified ims)

    # check inputs 
    ims = _check_inputs(ims,seg=seg)

    # grab first image and add to command list
    ref_im = ims.pop(0)
    cmd = ['itksnap', '-g', ref_im]

    # add more images if supplied 
    if ims:
        cmd.append('-o')
        cmd.extend(ims) # extend with added values

    # add segmentation if supplied
    if seg is not None:
        cmd.append('-s')
        cmd.append(seg)

    # end with ampersand so that we don't suspend ipython terminal 
    cmd.append('&')

    if not remote:
        # execute command in bash
        execute(cmd)
    else: 
        # print path to terminal instead so it can be executed via Terminal VNCserver, rather than directly in browser window (when remote working)
        out = ' '.join(cmd)
        print(out)
        return out  
    
def test_itksnap():
	""" Tests sv.itksnap() function"""

	# Specify example files 
	rootdir = '/fileserver/abd/serge/IVIM_data/IVIM_data/all_cases2/Case100/'
	impath1 = rootdir + 'average6/b0_averaged.nrrd'
	impath2 = rootdir + 'average6/b800_averaged.nrrd'
	seg = rootdir + 'segmentation.nrrd'


	# Single image as list 
	itksnap(ims=[impath1])
	print('test1 passed')

	# Single image as string 
	itksnap(ims=impath1)
	print('test2 passed')

	# Two images 
	itksnap(ims=[impath1,impath2])
	print('test3a passed')

	# Single image as string with segmentation 
	itksnap(ims=impath1, seg=seg)
	print('test3b passed')

	# Single image as list with segmentation 
	itksnap(ims=[impath1], seg=seg)
	print('test4 passed')

	# Two images with segmentation 
	itksnap(ims=[impath1, impath2], seg=seg)
	print('test5 passed')

	# Incorrect path for single image
	itksnap(ims=[impath1+'adsfad'])
	print('test6 passed')

	# Incorrect path for second image
	itksnap(ims=[impath1, impath2+'adsfad'])
	print('test7 passed')

	# Incorrect path for segmentation
	itksnap(ims=[impath1],seg=seg+'adf')
	print('test8 passed')
	
	# Remote working printout 
	itksnap(ims=[impath1, impath2], seg=seg, remote=True)
	print('test9 passed')

# -----------
# VS code IDE 
# -----------

def launch_json(s):

    """Convert Terminal command into a launch .json string"""

    # check that there aren't any $ signs in the strings as we need the full path
    assert '$' not in s, f"Please eliminate all $ signs by passing the string through 'echo' command first in Terminal."

    print("Paste the following into launch.json with \"args\": [OUTPUT_LIST] without 'python' and 'COMMAND_NAME' ")
    print(' ')

    print(json.dumps(s.split(' ')))


# -----------
# IVIM project related 
# -----------

def get_ivim_images(rootdir):
    """
    Load 4 IVIM parameter images from a directory. 

    Most IVIM tools* output 4 parameter files in a specific naming convention. 
    This script assumes that naming convention and loads all 4 files into a dict. 
    * such as FBM, DIPY and ROAR methods, used Vasylechko et al. 2020 MICCAI submission. 

    Args: 
        rootdir (str): path to directory that contains 4 IVIM parameter files. 
                            meanB0_1.nrrd 
                            meanADC_1.nrrd 
                            meanPER_1.nrrd 
                            meanPER_FRAC_1.nrrd
    Returns: 
        params_p (dict): dictionary with 4 parameter images of the IVIM 
        header (:obj:`collections.OrderedDict`):   header of the .nrrd file that was used to load the files 

    """    
    params_p = {'S0':None, 'D':None, 'P':None, 'P_f':None}
    params_p['S0'],header = nrrd.read(rootdir+"meanB0_1.nrrd")
    params_p['D'],_ = nrrd.read(rootdir+"meanADC_1.nrrd")
    params_p['P'],_ = nrrd.read(rootdir+"meanPER_1.nrrd")
    params_p['P_f'],_ = nrrd.read(rootdir+"meanPER_FRAC_1.nrrd")
    return params_p, header

def get_dwi_images(rootdir,suffix,bvalues=None):
    """
    Loads DWI images (raw MRI signal) into a list as numpy arrays. 

    Args: 
        rootdir (str):    path to directory that contains the DWI images at multiple bvalues 
        [suffix] (str):   each file is assumed to be prefixed with some identifier. 
                             E.g. b100_averaged.nrrd -> '_averaged' (default value)
        [bvalues] (list): list of bvalues that will be used to load the images. If not specified, the default values will be used. 

    Returns: 
        images (list):    bvalue images in a list 
        header (:obj:`collections.OrderedDict`): header of the .nrrd file that was used to load the files 

    """
    if not bvalues: 
        bvalues = [0,50,100,200,400,600,800]
    if not suffix: 
        suffix = "_averaged"
    images = []
    for bval in bvalues: 
        im, header = nrrd.read(rootdir+"b"+str(bval)+suffix+".nrrd")
        images.append(im)
    return images, header

def remove_outliers(data, max_deviations=2):
    """
    Remove outliers in a numpy array. 
    
    Args: 
        data (:obj:`numpy.ndarray`): data in a n-dimensional numpy array 
        [max_deviations] (int):      deviations from the mean for which to remove outliers. Defaults to 2. 
    Returns: 
        no_outliers (:obj:`numpy.ndarray`): data with removed outliers

    """
    
    mean = np.mean(data)
    standard_deviation = np.std(data)
    distance_from_mean = abs(data - mean)
    not_outlier = distance_from_mean < max_deviations * standard_deviation
    no_outliers = data[not_outlier]
    return no_outliers

def get_label_masks2(segmented_image,labels):     
    """ 
    Create binary mask for each ROI in a segmentation image as a dictionary. 
    
    Args: 
        segmented_image (str): path to a segmentation 
        labels (dict):         dictionary denoting the name and value for labels. 
                                e.g. dict({'kidneys':1,'spleen':2, 'liver':3})
    Returns: 
        label_masks (list):    list of numpy.ndarray objects as binary masks for each ROI 
        header (:obj:`collections.OrderedDict`): header of the .nrrd file that was used to load the files 
    """
    def _create_label_mask(segmented_image,label):
        # get label and turn into mask  
        mask = copy.deepcopy(segmented_image)
        mask = np.squeeze(mask)      
        mask[mask!=label] = 0  #all values not matching the label set to zero 
        mask[mask==label] = 1  #all values matching the label set to one 
        return mask    

    # load nii or nrrd
    if segmented_image.endswith('.nii.gz') or segmented_image.endswith('.nii'): 
        im = nb.load(segmented_image).get_fdata()
        header = nb.load(segmented_image).header
    elif segmented_image.endswith('.nrrd'): 
        im, header = nrrd.read(segmented_image)
    label_masks = {label:_create_label_mask(im, label_id) for label,label_id in labels.items()}
    return label_masks, header

def get_label_masks(segmented_image,labels=None):     
    """ 
    Create binary mask for each ROI in a segmentation image 
    
    Args: 
        segmented_image (str): path to a segmentation 
        labels (dict):         dictionary denoting the name and value for labels. 
                                e.g. dict({'kidneys':1,'spleen':2, 'liver':3})
    Returns: 
        label_masks (list):    list of numpy.ndarray objects as binary masks for each ROI 
        header (:obj:`collections.OrderedDict`): header of the .nrrd file that was used to load the files 
    """
    def _create_label_mask(segmented_image,label):
        # get label and turn into mask  
        mask = copy.deepcopy(segmented_image)
        mask = np.squeeze(mask)      
        mask[mask!=label] = 0  #all values not matching the label set to zero 
        mask[mask==label] = 1  #all values matching the label set to one 
        return mask    
    
    if not labels:
        labels = dict({'kidneys':1,'spleen':2, 'liver':3})
    # load nii or nrrd
    if segmented_image.endswith('.nii.gz') or segmented_image.endswith('.nii'): 
        im = nb.load(segmented_image).get_fdata()
        header = nb.load(segmented_image).header
    elif segmented_image.endswith('.nrrd'): 
        im, header = nrrd.read(segmented_image)
    label_masks = [_create_label_mask(im, label_id) for label,label_id in labels.items()]
    return label_masks, header

def ivimFBMMRFEstimator(bvalsFiles_average,mask,out_dir,fit_model,iterations=0,debug = False): 
    """
    Wrapper for FBM DIPY algorithm in python for a single image. 
    
    C++ implementation of FBM algorithms can be found here:
        /fileserver/abd/bin/ivimFBMMRFEstimator
    Note that the algorithm generates .vtk files. 
    
    Args: 
        bvalsFiles_average (str): path to a .txt file that indicates the full path to all of the 7 b-value images
        mask (str): path to a binary segmentation image in .vtk format 
        out_dir (str): path to a directory to save the files 
        fit_model (str): DIPY (BOBYQA) or FBM (spatially regularized BOBYQA)
        
    
    """
    func = "/fileserver/abd/bin/ivimFBMMRFEstimator"
    log = ""
    if fit_model == 'FBM':
        if iterations == 0: 
            iterations = 500 # set iterations to 500 if iterations are not supplied with the model (aka if they regress to default of zero)    
    cmd = [func,"--optMode","FBM","-n",str(7),"-i",bvalsFiles_average,"-g",str(iterations),"-o",out_dir,"-m", mask,log]
    if debug: 
        print(' '.join(cmd))
    #subprocess.call(cmd)      
    execute(cmd)
    
def write_bvalsFileNames_average(signaldir,extension='.vtk'):
    """create bvalFilenames_average .txt files required for running DIPY/FBM model 
    
    Args: 
        signaldir (str): path to directory which contains the acquired b-value files (whether geometrically averaged or not) in the form 'b0_averaged.vtk', etc 
        extension (str): specify whether the filesnames are .nrrd or .vtk (default)
    Returns: 
        savedir (str): directory to which the bvalsFileNames.txt file was saved. 

    WARNING: .txt file will be saved to the same directory where the images are 
    NB Advanced user warning: this is different to bvalsFilename.txt which is required to run geometric averaging operation. Do not confused the two. 

    
    """
    savedir = signaldir + "/bvalsFileNames.txt"
    lines = []
    
    with open(savedir,'w') as f:
        for bval in [0,50,100,200,400,600,800]:
            fullpath=signaldir+"b"+str(bval)+"_averaged"+extension
            lines.append('\t'.join([str(bval),fullpath+"\n"]))
        f.writelines(lines) 
    return savedir


def ivim_eqtn_image(parameters): 
    # get signal estimate from IVIM parameter estimates
    S0 = parameters['S0'][...,np.newaxis] #expand with one new axis so that b* parameter can be broadcast
    D = parameters['D'][...,np.newaxis]
    P = parameters['P'][...,np.newaxis]
    P_f = parameters['P_f'][...,np.newaxis]
    x,y,z,t = S0.shape
    bvals = np.array([0,50,100,200,400,600,800])
    bvals_ = np.tile(bvals,(x,y,z,1))
    
    signal_image = S0*(P_f*np.exp(-1*bvals_*P)+(1-P_f)*np.exp(-1*bvals_*D))
    signal_image = np.moveaxis(signal_image,3,0)
    return signal_image


def ivim_eqtn_voxel(parameters): 
    # put parameter estimates of a voxel through the forward model to obtain the estimated signal 
    S0 = parameters['S0']
    D = parameters['D']
    P = parameters['P']
    P_f = parameters['P_f']
    bvals = np.array([0,50,100,200,400,600,800])
    signal = S0*(P_f*np.exp(-1*bvals*P)+(1-P_f)*np.exp(-1*bvals*D))
    return signal
    

# -----------
# Misc
# -----------


def msg(msg='', buffer=3,abspath=True):
    
    """ 
    This function is supposed to be used as a custom assert helper function. Allows for interactive python usage. 
    
    Example: 
        from IPython import embed; embed 
        
        assert 5==6, embed(header=sv.msg('Something went wrong'))
        
    This funciton prints line number and filename where a given piece of code was executed, together with a customer message.
    
    """
    
    def _get_line(line, buffer=3):
        """Get specific lines of code from a filename at a given line number. 
        
        Specify 'buffer' to get more lines on either side of the given line. 
        """
        assert os.path.exists(filename), f"{filename}"
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        
        lines_subset = lines[(line-1)-buffer:(line-1)+buffer]
        
        lines_subset = ' '.join(lines_subset)
        
        return lines_subset
    
    filename = getframeinfo(currentframe().f_back).filename # get filename 
    line = currentframe().f_back.f_lineno # gets line number 
    code = _get_line(line, buffer=buffer) # gets the code lines 
    
    if not abspath: 
        filename = os.path.basename(filename) # print only the filename, not the full path 
    
    output = f"Error on line: {line} of {filename}. \n\n{code}\n\n"
    
    if msg:
        # add extra message
        output = output + f"{msg}"

    return output



# #########################################################################################################################
# -------------------------------------------------------------------------------------------------------------------------
# #########################################################################################################################
#
# 
# 
# SECTION 2 - WIP 
#
#
#
# #########################################################################################################################
# -------------------------------------------------------------------------------------------------------------------------
# #########################################################################################################################



# -----------
# WIP 
# -----------

def sort_glob_slices(files,refpattern='0000'):
    
    #unfinished due to pressing deadline - return later 
    """Fetch files via `glob` and sort them in ascending order numerically by slice number
    
    Args:
        refpattern (str): reference pattern to match the slices to, for example '0000' would use this slice as reference and sort all files numerically in this substring
        
        files (str): list of files found by `glob` function 
    
    """      
    
    """ Debug (bash)
    
    rootdir='/home/ch215616/abd/mwf_data/mwf_maps_all_patients/volunteer/synthnet/simple_unet/'
    files = glob.glob(rootdir+'*.nii.gz')
    
    references = [f for f in files if refpattern in f]
    
    for r in references:
        slices = glob.glob(r.replace(refpattern,"*"))
        
        # now extract biggest and smallest number in the 'refpattern', then iterate over these numbers and fetch files
    
    
    
    """
    pass
            

def fix_random_seed():

	# Seed value
	# Apparently you may use different seed values at each stage
	seed_value= 0

	# 1. Set the `PYTHONHASHSEED` environment variable at a fixed value
	import os
	os.environ['PYTHONHASHSEED']=str(seed_value)

	# 2. Set the `python` built-in pseudo-random generator at a fixed value
	import random
	random.seed(seed_value)

	# 3. Set the `numpy` pseudo-random generator at a fixed value
	import numpy as np
	np.random.seed(seed_value)

	# 4. Set the `tensorflow` pseudo-random generator at a fixed value
	import tensorflow as tf
	tf.random.set_seed(seed_value)
	# for later versions: 
	# tf.compat.v1.set_random_seed(seed_value)

	# 5. Configure a new global `tensorflow` session
	from keras import backend as K
	session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
	sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
	K.set_session(sess)
	# for later versions:
	# session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
	# sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
	# tf.compat.v1.keras.backend.set_session(sess)

	# https://stackoverflow.com/questions/32419510/how-to-get-reproducible-results-in-keras
