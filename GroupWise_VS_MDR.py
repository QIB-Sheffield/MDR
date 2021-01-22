"""
Spyder Editor

This is a temporary script file.
"""
###############################################################################
### IMPORTS:

from collections import OrderedDict
import glob
from zipfile import ZipFile
from tqdm import tqdm
import os
import shutil
import skimage.io as io
import numpy as np
import SimpleITK as sitk
import subprocess
import sys  
import pydicom
from scipy.optimize import curve_fit
import time
import cv2
from numpy import trapz
from scipy.signal import argrelextrema
import copy
import pandas as pd
import multiprocessing as mp
from scipy.stats import iqr
np.set_printoptions(threshold=sys.maxsize)

###############################################################################
### GENERIC FUNCTIONS:

def unzip_patients(path):
    """
    This function will find all the .zip files in the given "path" folder
    and will extract each folder contained in the zip in "path". Once the .zip
    folders are extracted both them and the parent folders will be deleted
    hence only folders with prefix Leeds_Patients_ will be kept.
    """
    zip_folders = glob.glob(path+'/*.zip')
    print('\nExtracting patients images from .zip files: \n')
    if zip_folders:
        for file in tqdm(zip_folders):
            zip_file = ZipFile(file, 'r')
            zip_file.extractall(path)
            names_list = zip_file.namelist()
            temp_name = names_list[-1].split('/')[0:2]
            file_to_extract = '/'.join(temp_name)+'/'
            print(file_to_extract, path)
            zip_file.extract(names_list[-1], path)
            zip_file.close()
            os.remove(file)
            print('\n', file.split('\\')[-1], ' was extracted and deleted\n')
    else:
        print('No .zip folders were found in the given path')
    
    # Deleting the parent folders in order to keep only the Leeds_Patient_... 
    # folders
    os.chdir(path)
    parent_folders = os.listdir()
    for folder in parent_folders:
        os.chdir(folder)
        for inner_folder in os.listdir():
            shutil.move(path+'\\'+folder+'\\'+inner_folder, path+'\\'+inner_folder)
        os.chdir('..')
        shutil.rmtree(folder)
            
def create_slice_folders(sequence, slices):
    """
    This function will create a number of folders equal to "slices" 
    (in the current directory) with name in the form of: 
    "sequence"_slice_1, "sequence"_slice_2 etc 
    where "sequence" can be DTI, DCE etc
    """
    for i in range(slices):
        if not os.path.isdir(sequence+'_slice_'+str(i+1)):
            os.mkdir(sequence+'_slice_'+str(i+1))
        else:
            print(sequence+'_slice_'+str(i+1), ' already exists')
        
def move_images_to_slice_folders(sequence, dcm_files_found, slices, num_of_files, folder):
    """
    This function will take the images in the form of .dcm and will re-arrange
    them according to the slice they came from.
    It assumes that the slice_folders are already created and they are at 
    the same file with the "dcm_files_found".
    """
    repetitions = num_of_files/slices
    assert int(repetitions) == repetitions, "The number of .dcm images found in %s can not be divided to the number of slices expected" %os.getcwd()    
    images_path = os.getcwd()
    for image_name in dcm_files_found:
        print(image_name.lower())
        if image_name.lower().endswith('.dcm'):
            idx = int(image_name.split('-')[2])
        elif image_name.lower().endswith('.ima'):
            idx = int(image_name.lower().split('.')[3])            
        for check_slice in reversed(range(1, slices+1)):
            if (idx - check_slice) % slices == 0:
                slice_folder = sequence + '_slice_' + str(check_slice)
                shutil.move(images_path+'\\'+image_name, images_path+'\\'+slice_folder+'\\'+image_name)
    assert len(os.listdir())==slices, "Some images were not re-arranged to slice_folders check %s" %images_path 
   
def load_mhd3d(three_d_mhd_path):
    """
    This function loads a single 3d mhd file found in three_d_mhd_path
    """
    three_d_mhd = io.imread(three_d_mhd_path, plugin='simpleitk')
    return three_d_mhd
    
def find_and_load_as_mhd(folder_path, suffix, crop_flag = True):
    """
    This function will find all the paths of the files ending with 'suffix' 
    found in 'folder_path', then they get sorted and returned all in a single
    variable called three_d_mhd. If crop_flag is left to True automatic
    cropping will take place, if it is set to False no cropping takes place.
    """
    if suffix == '.dcm':
        img_paths = find_and_sort_paths(folder_path, suffix)
        if not img_paths:
            img_paths = find_and_sort_paths(folder_path, '.ima') 
    elif suffix == '.mhd':
        img_paths = glob.glob(folder_path + '\*' + suffix)
    three_d_mhd = []
    for idx, dcm_path in enumerate(img_paths):
        current_image = io.imread(dcm_path, plugin='simpleitk').squeeze()
        if crop_flag:
            cropped_image = cv2.resize(current_image, dsize=(192, 192))
        else:
            cropped_image = current_image
        if idx==0:
            three_d_mhd = np.zeros(np.shape(cropped_image) + (len(img_paths),))
        three_d_mhd[:, :, idx] = cropped_image
    return three_d_mhd
  
def elastix_registration(moving_image_path, fixed_image_path, output_dir, parameters_file):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    cmd = [ 'elastix', '-m', moving_image_path, '-f', fixed_image_path, '-out', output_dir, '-p', parameters_file]
    try:
        subprocess.check_call(cmd)
    except:
        print ('Image registration failed')
        print (sys.exc_info())

def rewrite_mhd(mhd_path, new_raw_name):
    """
    This function searches in the "mhd_path" to find the end result of 
    elastix registration process with name result.0.raw and transforms it
    into an .mhd with name: "new_raw_name"
    """
    with open(mhd_path, 'r') as file :
        filedata = file.read()
    filedata = filedata.replace('result.0.raw', new_raw_name)                        
    # Write the file out again
    with open(mhd_path, 'w') as file:
        file.write(filedata)
        
def clean_up_elastix_files(directory):
    """
    Delete files created by elastix (every file except .mhd and .raw) found 
    in 'directory'.
    """
    os.remove(directory + '\\elastix.log')
    os.remove(directory + '\\IterationInfo.0.R0.txt')
    os.remove(directory + '\\IterationInfo.0.R1.txt')
    try:
        os.remove(directory + '\\IterationInfo.0.R2.txt')
    except:
        pass
    try:
        os.remove(directory + '\\IterationInfo.0.R3.txt')
    except:
        pass

def clean_up_all_files_in_folder(directory):
    """
    Delete all files found in output_dir.
    """
    files = glob.glob(directory + '\*')
    for file in files:
        os.remove(file)
  
def find_and_sort_paths(slice_path, suffix):
    """
    It finds all paths in 'slice_path' ending with 'suffix' and gets them sorted
    Attention: in case .dcm is given but no folder ending in .dcm is found it searches for .ima
    """
    slice_dcms = glob.glob(slice_path + '\*' + suffix)
    if not slice_dcms and suffix=='.dcm':
        slice_dcms = glob.glob(slice_path + '\*' + '.ima')
    indexes = []
    for dcm in slice_dcms:
        if dcm.lower().endswith('.dcm'):
            idx = int(dcm.split('-')[2])
        elif dcm.lower().endswith('.ima'):
            idx = int(dcm.split('.')[3])
        indexes.append(idx)
    combined = list(zip(slice_dcms, indexes))
    combined.sort(key = lambda this: this[1])
    sorted_dcms, _ = zip(*combined)
    return sorted_dcms

def create_MoCoMo_folders(MoCoMo_slice_path):
    """
    Checks if the needed folders for MoCoMo already exist and if they don't
    exist it creates them.
    """
    if not os.path.isdir(MoCoMo_slice_path):
        os.mkdir(MoCoMo_slice_path)
    if not os.path.isdir(MoCoMo_slice_path+'\\Original'):
        os.mkdir(MoCoMo_slice_path+'\\Original')
    if not os.path.isdir(MoCoMo_slice_path+'\\Fitted'):
        os.mkdir(MoCoMo_slice_path+'\\Fitted')
    if not os.path.isdir(MoCoMo_slice_path+'\\BSplines_Registered_1'):
        os.mkdir(MoCoMo_slice_path+'\\BSplines_Registered_1')
    
def handle_original_mhd(slice_path, suffix = '.dcm', crop_flag = True):
    '''
    This function takes as input 'slice_path' where images ending with 'suffix'
    are expected to be found and loaded. It returns a matrix of size: (width*height, timepoints) 
    also returns the shape ("initial_shape") of the images found in slice_path.
    '''
    initial_images = find_and_load_as_mhd(slice_path, suffix, crop_flag)
    if len(initial_images)==0:
        initial_images = find_and_load_as_mhd(slice_path, '.ima', crop_flag)        
    initial_shape = np.shape(initial_images)    
    for i in range(1, initial_shape[-1]+1):
        if sequence=='DCE':
            sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(cv2.resize(initial_images[:,:,i-1], dsize=(192, 192))),sitk.sitkFloat64), MoCoMo_slice_path + '\\Original\\' + str(format(i, "03")) + '.mhd')
        else:
            sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(initial_images[:,:,i-1]),sitk.sitkFloat64), MoCoMo_slice_path + '\\Original\\' + str(format(i, "03")) + '.mhd')
            
    initial_images = initial_images.reshape((initial_shape[0]*initial_shape[1]), initial_shape[2])
    return initial_images, initial_shape

def transformix_deformation_field(output_dir, transform_parameters_file):
    cmd = [ 'transformix', '-def', 'all', '-out', output_dir, '-tp', transform_parameters_file]
    try:
        subprocess.check_call(cmd)
    except:
        print ('Transformix failed')
        print (sys.exc_info())

def MoCoMo_elastix_registration_wrapper(moving_images_paths, fixed_images_paths, output_dir, parameters_file):
    '''
    Wrapper function for elastix_registration it takes as inputs:
    'moving_images_paths': the full paths for the fixed images of the pairwise registrations 
    'fixed_images_paths': the full paths for the moving images of the pairwise registrations
    'output_dir': the full path for the folder where the output of the registration outcomes will be saved
    'parameters_file': the file of the parameters used for the specific elastix registration
    '''
    moving_images_paths = sorted(moving_images_paths)
    fixed_images_paths = sorted(fixed_images_paths)
    for i in range(len(moving_images_paths)):
        moving_image_path = moving_images_paths[i]
        fixed_image_path = fixed_images_paths[i]
    
        elastix_registration(moving_image_path, fixed_image_path, output_dir, parameters_file)
        
        # Renaming, change mhd content and delete
        os.rename(output_dir+'\\'+'result.0.raw', output_dir+'\\'+moving_image_path.split('\\')[-1].replace('.mhd', '.raw'))
        os.rename(output_dir+'\\'+'result.0.mhd', output_dir+'\\'+moving_image_path.split('\\')[-1])
        os.rename(output_dir+'\\'+'TransformParameters.0.txt', output_dir+'\\'+'TransformParameters.0.txt'.replace('0','{0:03d}'.format(i+1)))
        rewrite_mhd(output_dir+'\\'+moving_image_path.split('\\')[-1], moving_image_path.split('\\')[-1].replace('.mhd', '.raw'))
        
        clean_up_elastix_files(output_dir)
            
def cropping_dimensions(sequence):
    if sequence=='T1':
        rows_start = 52
        rows_end = 52+260
        cols_start = 62
        cols_end = 62+260
        
    elif sequence=='DTI':
        rows_start = 21
        rows_end = 21+130
        cols_start = 21
        cols_end = 21+130
        
    elif sequence=='DCE':
        rows_start = 21
        rows_end = 21+130
        cols_start = 31
        cols_end = 31+130
    
    return rows_start, rows_end, cols_start, cols_end

def draw_the_grid(img, lines_per_dimension=12):
    maximum_value = np.max(img)
    shape = np.shape(img)
    num = int(shape[0]/lines_per_dimension)
    for i in range(lines_per_dimension+1):
        cv2.line(img, (int(i*num), 0),(int(i*num), int(shape[0])), (int(maximum_value)), thickness=2)
        cv2.line(img, (0, int(i*num)), (int(shape[0]), int(i*num)), (int(maximum_value)), thickness=2)
    return img

def create_grid(images, full_output_path):
    images = copy.deepcopy(np.moveaxis(images, 0, -1))
    for i in range(np.shape(images)[-1]):
        images[:,:,i] = draw_the_grid(images[:,:,i])
    images = np.moveaxis(images, -1, 0)
    sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(images),sitk.sitkFloat64), full_output_path)       

def change_elastix_naming_result(full_path_mhd, str_to_be_replaced, new_name):
    counter_file = open(full_path_mhd, 'r+')
    content_lines = []
    for line in counter_file:
        if str_to_be_replaced in line:
            line = line.replace(str_to_be_replaced, new_name)
        content_lines.append(line)

    counter_file.seek(0)
    counter_file.truncate()
    counter_file.writelines(content_lines)
    counter_file.close()

##############################################################################
### T1 FUNCTIONS:  

def read_inversion_times(sorted_dcms_paths):
    """
    This function reads the inversion times for T1 sequences.
    It takes as argument a sorted list of dcm paths "sorted_dcms_paths"
    and returns the correspondent list of inversion times
    """
    inversion_times = []
    for dcm_path in sorted_dcms_paths:
        dataset = pydicom.dcmread(dcm_path)
        inversion_times.append(dataset.InversionTime)
    return inversion_times
    
def reverse_T1_signals(images):
    """
    It expects a [h*w, t] shaped matrix as input
    This function finds the local minimum (null point) for each T1 signal of 
    the image and reverses the values for all the points before the minimum 
    (to properly reverse them you must first change the signs of all the 
    values before the minimum and then add the value of the minimum itself). 
    The first point is not considered as possible candidate for the null point).
    """
    possible_null_point = images[:, 1:]
    value_to_add = []
    # Having left the first element out of the search for the minimum you need
    # to add 1 to the col index 
    col = np.add(np.argmin(possible_null_point, axis=1),1)
    new_images = np.zeros_like(images)
    for i in tqdm(range(np.shape(images)[0])):
        new_images[i, 0:col[i]] = -images[i, 0:col[i]]+ images[i, col[i]]
        value_to_add.append(images[i, col[i]])
        new_images[i, col[i]:] = images[i, col[i]:]
    return new_images, col, value_to_add

def exp_func(x, a, b, T1):
    """
    Exponential function used to perform T1-fitting
    """
    return a - b * np.exp(-x/T1)

def T1_fitting_for_parallel(signal, inversion_times, p0, lb, ub):
    try:
        popt, pcov = curve_fit(exp_func, xdata = inversion_times, ydata = signal, p0 = p0, bounds = (lb, ub), method = 'trf') 
        success = True
    except:
        # Here the ub gets overrided
        # For GroupWise:
        popt = [0, 0, 1] #when calculating T1_map it will assign 0 to this pixel
        success = False
    fitted = []
    for it in inversion_times:
        fitted.append(exp_func(it, *popt))

    # T1_estimated = T1 * ( b/a - 1):
    if success:
        T1_estimated = popt[2] * ((popt[1] / popt[0]) -1)
    else:
        T1_estimated = 0
    return fitted, T1_estimated, popt[2], popt[1], popt[0]
    
def parallel_T1_fitting(images, output_path, inversion_times, initial_shape, save_flag=True):
    '''
    images: is of shape [h*w, t], while initial_shape: shows the actual 3d shape of the image stack
    '''
    # Perform the reverse for T1:
    p0 = [np.max(images), np.max(images)-np.min(images), 50]
    images, inversion_points, value_to_add = reverse_T1_signals(images)
    
    # T1-Fitting:                  
    #p0 = [np.max(images), np.max(images)-np.min(images), 50]
    lb = [0, 0, 0]
    ub = [np.inf, np.inf, 2000]
      
    # Run in parallel:
    indices = range(initial_shape[0]*initial_shape[1])
    # construct pool
    pool = mp.Pool(mp.cpu_count())
    print('Parallel Processing for T1 fitting ...' )
    T1_estimated = np.zeros((initial_shape[0]*initial_shape[1],))
    T1_aparent = np.zeros((initial_shape[0]*initial_shape[1],))
    B = np.zeros((initial_shape[0]*initial_shape[1],))
    A = np.zeros((initial_shape[0]*initial_shape[1],))
    fitted = np.zeros((initial_shape[0]*initial_shape[1], initial_shape[2]))
    args = [(images[idx, :], inversion_times, p0, lb, ub) for idx in indices]
    results = pool.starmap(T1_fitting_for_parallel, args)        

    for index, result in zip([idx for idx in tqdm(range(len(indices)))], results):
        T1_estimated[index] = result[1]
        fitted[index,:] = result[0]
        T1_aparent[index] = result[2]
        B[index] = result[3]
        A[index] = result[4]
    fitted = np.array(fitted)
    for i in range(np.shape(fitted)[0]):
        fitted[i, 0:inversion_points[i]] = -fitted[i, 0:inversion_points[i]] + value_to_add[i] 
    fitted = fitted.reshape(initial_shape)
    
    T1_estimated = np.array(T1_estimated)
    pool.close()
    
    if save_flag and output_path!='':
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        for i in range(1, initial_shape[-1]+1):
            sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(fitted[:,:,i-1]),sitk.sitkFloat64), output_path + str(format(i, "03")) + '.mhd')
    return fitted.reshape(initial_shape), T1_estimated.reshape((initial_shape[0], initial_shape[1])), [T1_aparent.reshape((initial_shape[0], initial_shape[1])), B.reshape((initial_shape[0], initial_shape[1])), A.reshape((initial_shape[0], initial_shape[1]))], images.reshape(initial_shape)
    
###############################################################################
### DTI FUNCTIONS:  
    
def read_dicom_tags_DTI(sorted_dcms_paths):
    """
    This function reads the b-values,original b-vectors and the patient orientation for DTI sequences.
    It takes as argument a sorted list of dcm paths "sorted_dcms_paths"
    and returns: 
    the correspondent list of b-values "b_values"
    the original b-vectors "b_Vec_original"
    the patient's orientation "image_orientation_patient"
    """
    b_values = []
    b_Vec_original = []
    image_orientation_patient = []
    pixel_spacing = []
    slice_thickness = []
    for dcm_path in sorted_dcms_paths:
        dataset = pydicom.dcmread(dcm_path)
        b_values.append(dataset[0x19, 0x100c].value)
        b_Vec_original.append(dataset[0x19, 0x100e].value)
        image_orientation_patient.append(dataset.ImageOrientationPatient)
        pixel_spacing.append(dataset.PixelSpacing)
        slice_thickness.append(dataset.SliceThickness)
    return b_values, b_Vec_original, image_orientation_patient, pixel_spacing, slice_thickness

def mono_exp_model(x, b):
    return x[:, 6] * np.exp(np.matmul(-x[:,1:6], b))
    
def DTI_fitting(im, b, thresh_val, method='linear'):
    sz = np.shape(im)
    # handle a vector with the wrong orientation
    if len(sz)==2 and sz[1]==1:
        im = np.transpose(im)
        sz = np.shape(im)

    # reshape to a matrix
    im_matrix = np.reshape(im, (-1,sz[-1]))
    
    # get the mask
    mask = np.array(thresh_val, dtype=bool)

    # take only the voxels inside the mask
    mask = np.reshape(mask,(-1, sz[-1]))
    I = im_matrix[mask]
    
    if not np.all(np.isreal(I)):
        print('Some voxels are complex. Taking magnitude.')
        I = np.abs(I)
        
    # take the log of the image to linearise the equation
    abs_I = np.abs(I)
    imlog = np.ma.log(abs_I)
    imlog = np.reshape(imlog, (-1, sz[-1]))
    
    # Sort all b matrices in to a vector Bv=[Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz];
    Bv = np.vstack((b[0,0,:],
                    2*b[0,1,:],
                    2*b[0,2,:],
                    b[1,1,:],
                    2*b[1,2,:],
                    b[2,2,:]))
    
    Bv = np.transpose(Bv)
    
    # Add another column to Bv to handle the constant term:
    # Slog = Bv * M + log(S0)
    # becomes:
    # Slog = [Bv, -1] * [M; -log(S0)]
    minus_one_column = -np.ones((np.shape(Bv)[0]))
    Bv_new = np.c_[Bv, minus_one_column]
    assert method=='linear', "Wrong method as argument!!!"
    M = np.linalg.lstsq(Bv_new, -imlog.T, rcond=None)
    M = M[0].T
    
    M[np.isnan(M)]=0
    M[np.isinf(M)]=0
    ### Initialize Variables
    FA_mask = np.empty(sz[0]*sz[1])
    ADC_mask = np.empty(sz[0]*sz[1])

    start = time.time()
    for i in range(np.shape(M)[0]):
        
        # The DiffusionTensor (Remember it is a symetric matrix,
        # thus for instance Dxy == Dyx)
        # DiffusionTensor=[Mi[0] Mi[1] Mi[2]; Mi[1] Mi[3] Mi[4]; Mi[2] Mi[4] Mi[5]]
        DiffusionTensor = np.zeros((3,3))
        DiffusionTensor[0][0] = M[i, 0]
        DiffusionTensor[0][1] = M[i, 1]
        DiffusionTensor[0][2] = M[i, 2]
        DiffusionTensor[1][0] = M[i, 1]
        DiffusionTensor[1][1] = M[i, 3]
        DiffusionTensor[1][2] = M[i, 4]
        DiffusionTensor[2][0] = M[i, 2]
        DiffusionTensor[2][1] = M[i, 4]
        DiffusionTensor[2][2] = M[i, 5]

        # Calculate the eigenvalues and vectors, and sort the 
        # eigenvalues from small to large
        [EigenValues, EigenVectors]=np.linalg.eig(DiffusionTensor);
        if np.sum(EigenValues)!=0:
            EigenValues, EigenVectors = zip(*sorted(zip(EigenValues, EigenVectors)))
        
        # Regulating of the eigen values (negative eigenvalues are
        # due to noise and other non-idealities of MRI)
        EigenValues=np.abs(EigenValues)

        # Apparent Diffuse Coefficient
        ADCv = np.mean(EigenValues)

        # FA definition:
        denominator = np.sqrt(EigenValues[0]**2+EigenValues[1]**2+EigenValues[2]**2)
        if denominator == 0:
            FA_mask[i] = np.nan
        else:    
            FA_mask[i]=np.sqrt(1.5)*(np.sqrt((EigenValues[0]-ADCv)**2+(EigenValues[1]-ADCv)**2+(EigenValues[2]-ADCv)**2)/denominator)
        ADC_mask[i]=ADCv

    stop = time.time()
    print('Took: ', stop-start)
    
    Bv_new_times_M_new = np.moveaxis(np.matmul(Bv_new, M.T),0,-1).reshape(sz)
    Fitted = np.exp(-Bv_new_times_M_new)
    
    return M, Bv_new, FA_mask, ADC_mask, Fitted

def DTI_fitting_wrapper(corrected_images, output_path, b_values, bVec_original, image_orientation_patient, pixel_spacing, slice_thickness, save_flag=True):
    sz = np.shape(corrected_images)
    
    for i in range(len(image_orientation_patient)-1):
        assert image_orientation_patient[i] == image_orientation_patient[i+1], "Error in image_orientation_patient for DTI"
    R1 = image_orientation_patient[0][3:6]
    R1 = [-float(x) for x in R1]
    R2 = image_orientation_patient[0][0:3]
    R2 = [-float(x) for x in R2]
    R3 = np.cross(R1, R2)
    R = np.vstack((R1, R2, R3))
    bVec = np.dot(R, np.array(bVec_original).T).T
    for i in range(len(pixel_spacing)-1):
        assert pixel_spacing[i] == pixel_spacing[i+1], "Error in pixel_spacing in DTI"
    x_res, y_res = pixel_spacing[0]
    x_res = float(x_res)
    y_res = float(y_res)
    for i in range(len(slice_thickness)-1):
        assert slice_thickness[i] == slice_thickness[i+1], "Error in slice_thickness in DTI"
    
    ### Mask
    sz_mask = np.shape(corrected_images)
    mask = np.ones(sz_mask)
    ### Fitting
    B = np.zeros((3, 3, len(b_values)))
    for idx_b in range(len(b_values)):
        B[:, :, idx_b] = np.outer(np.outer(b_values[idx_b], bVec[idx_b,:].T), bVec[idx_b,:])

    [M_new, Bv_new, fa, adc, Fitted] = DTI_fitting(corrected_images, B, mask, 'linear')
    
    if save_flag and output_path!='':
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        for i in range(1, sz[-1]+1):
            sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(Fitted[:,:,i-1]),sitk.sitkFloat64), output_path + str(format(i, "03")) + '.mhd')
         
    return Fitted, M_new.T[:, :], Bv_new[:, :], fa, adc

##############################################################################    
### DCE FUNCTIONS:    

def load_txt(full_path_txt):
    counter_file = open(full_path_txt, 'r+')
    content_lines = []
    for cnt, line in enumerate(counter_file):
        content_lines.append(line)
    x_values_index = content_lines.index('X-values\n')
    assert (content_lines[x_values_index+1]=='\n')
    y_values_index = content_lines.index('Y-values\n')
    assert (content_lines[y_values_index+1]=='\n')
    time = list(map(lambda x: float(x), content_lines[x_values_index+2 : y_values_index-1]))
    aif = list(map(lambda x: float(x), content_lines[y_values_index+2 :]))
    return aif, time

def Integral_Trapezoidal_Rule_initial(x,time):
    first_pass=[]
    for tm in range(len(time)):
        first_pass.append(trapz(x[0:tm+1],time[0:tm+1]))
    return first_pass

def Integral_Trapezoidal_Rule_second(first_pass,time):
    second_pass=[]
    for tm in range(len(time)):
        second_pass.append(trapz(first_pass[0:tm+1],time[0:tm+1]))
    return second_pass

def Linear_Least_Squares_2CFM(ct,time,aif, timepoint = 39,Hct = 0.45):
    ct0 = np.mean(ct[0:timepoint])
    aif0 =  np.mean(aif[0:timepoint])
    ct_new = ct-ct0
    aif_new = (aif-aif0)/(1-Hct)

    #initialization of matrix A
    A=np.array([0,0,0,0]) 

    first_pass_ct_new=Integral_Trapezoidal_Rule_initial(ct_new,time)
    first_pass_aif_new=Integral_Trapezoidal_Rule_initial(aif_new,time)
    second_pass_ct_new=Integral_Trapezoidal_Rule_second(first_pass_ct_new,time)
    second_pass_aif_new=Integral_Trapezoidal_Rule_second(first_pass_aif_new,time)
    for t in range(0,len(time)):
        A1_1=second_pass_ct_new[t]
        A1_2=first_pass_ct_new[t]
        A1_3=second_pass_aif_new[t]
        A1_4=first_pass_aif_new[t]
        A_next=np.array([-A1_1,-A1_2,A1_3,A1_4])
        A=np.vstack((A,A_next))    

    # Drop first row of matrix [A] which is full of zeros and was used for initialization purposes
    A=np.delete(A,(0),axis=0)
    
    X = np.linalg.lstsq(A,ct_new,rcond=None)[0]

    # Extract physical parameters
    alpha = X[0]
    beta = X[1]
    gamma = X[2]
    Fp = X[3]
    
    if alpha == 0 or Fp == 0 :
        Tp = 0
        PS = 0
        Te = 0
        fit = 0
    else:    
        T = gamma/(alpha*Fp)  
        det = np.square(beta)-4*alpha
        if det < 0 :
            Tp = beta/(2*alpha)
            Te = beta/(2*alpha)
        else:
            Tp = (beta - np.sqrt(np.square(beta)-4*alpha))/(2*alpha)
            Te = (beta + np.sqrt(np.square(beta)-4*alpha))/(2*alpha)
        if Te == 0:
            Fp = 0
            Tp = 0
            PS = 0
            fit = 0
        else:    
            PS = Fp*(T-Tp)/Te                       
            #params = [Fp, Tp, PS, Te]
            fit =X[0]*A[:,0] + X[1]*A[:,1] + X[2]*A[:,2] + X[3]*A[:,3]
            fit=ct0+fit
    
    BF = Fp/(1 -Hct)	#Blood Flow
    
    if isinstance(fit, int):
        fit = np.zeros_like(BF)

    return BF, fit, Fp, Tp, PS, Te 

def DCE_find_timepoint(aif):
    # find timepoint
    idx = np.argmax(aif)
    minima = argrelextrema(np.asarray(aif), np.less)[0]
    maxima = argrelextrema(np.asarray(aif), np.greater)[0]
    assert idx in maxima, "The global maximum was not found in the maxima list"
    wanted_minimum = np.nan
    for ix, val in enumerate(minima):
        if ix==len(minima)-1:
            continue
        elif val<idx and minima[ix+1]>idx:
            wanted_minimum = val
    assert wanted_minimum + 1 < idx, "Minimum + 1 position equals to Global Maximum"
    return wanted_minimum + 1  

#################################################################################    
### Evaluation metrics:

def coefficinet_of_variation(image_stack):
    """
    Takes as input an 'image_stack' (3d stack of 2d images)
    and calculates the coefficient of variation = 1/tSNR = std/mean
    Higher tSNR is better
    Lower coefficient of variation is better
    """
    means = np.nanmean(image_stack, axis=2)
    stds = np.nanstd(image_stack, axis=2)
    non_zero_means_mask = np.zeros(np.shape(stds), dtype=bool)
    non_zero_means = np.where(means!=0)
    non_zero_means_mask[non_zero_means[0], non_zero_means[1]]=True
    c_o_v = np.divide(stds, means, where=non_zero_means_mask)
    c_o_v[~non_zero_means_mask] = np.nan           
    return np.round(np.nanmean(c_o_v),3), np.round(np.nanstd(c_o_v),3), np.round(np.nanmedian(c_o_v),3), np.round(iqr(c_o_v, nan_policy='omit'),3), np.round(c_o_v,3)
 
    
def iqr_calculation(image_stack):
    iqrs = iqr(image_stack, axis=2, nan_policy='omit')
    return np.round(np.nanmean(iqrs),3), np.round(np.nanstd(iqrs),3), np.round(np.nanmedian(iqrs),3), np.round(iqr(iqrs, nan_policy='omit'),3)

        
def mean_of_std_2nd_derivative(image_stack):
    """
    Takes as input an 'image_stack' (3d stack of 2d images)
    and calculates the mean of the std of the 2nd derivative of the pixels
    timeseries.
    Lower mean_of_std_2nd_derivative is better
    """
    shape = np.shape(image_stack)
    image_stack = image_stack.reshape((shape[0]*shape[1], shape[2]))
    gradient_1 = np.zeros((shape[0]*shape[1], shape[2]-1))
    gradient_2 = np.zeros((shape[0]*shape[1], shape[2]-2))
    for i in range(shape[0]*shape[1]):
        gradient_1[i,:] = np.diff(image_stack[i,:])
        gradient_2[i,:] = np.diff(gradient_1[i,:])
    return np.round(np.nanmean(np.nanstd(gradient_2, axis=1)),3), np.round(np.nanstd(np.nanstd(gradient_2, axis=1)),3), np.round(np.nanmedian(np.nanstd(gradient_2, axis=1)),3), np.round(iqr(np.nanstd(gradient_2, axis=1), nan_policy='omit'),3), np.nanstd(gradient_2, axis=1).reshape((shape[0], shape[1]))
        
   
def chi_squared_normalised(observed, expected):
    shape = np.shape(observed)
    my_chi_squared = np.zeros((shape[0]*shape[1], ))
    observed = observed.reshape((-1, shape[-1]))
    expected = expected.reshape((-1, shape[-1]))    
    for i in range(shape[0]*shape[1]):
        my_chi_squared[i] = np.round(((observed[i,:]-expected[i,:])**2).sum()/(expected[i,:]**2).sum(),3)
    my_chi_squared = my_chi_squared.reshape((shape[0], shape[1]))
    return my_chi_squared


def calculate_all_metrics(img, sharpness):
    img = np.moveaxis(img, 0, -1)
    dictionary = OrderedDict()    
    dictionary['mean_coefficient_of_variation'], dictionary['std_coefficient_of_variation'],dictionary['median_coefficient_of_variation'], dictionary['iqr_coefficient_of_variation'], dictionary['coefficient_of_variation'] = coefficinet_of_variation(img)
    dictionary['mean_std_sec_der'], dictionary['std_std_sec_der'], dictionary['median_std_sec_der'], dictionary['iqr_std_sec_der'], dictionary['std_sec_der'] = mean_of_std_2nd_derivative(img)
    dictionary['mean_iqr'], dictionary['std_iqr'], dictionary['median_iqr'], dictionary['iqr_iqr'] = iqr_calculation(img)
    return dictionary

        
def sharpness(img):
    # not necessary function sharpness was not used in any evaluation metric
    gy, gx = np.gradient(img)
    gnorm = np.sqrt(gx**2 + gy**2)
    return np.mean(gnorm)

#################################################################################    
### GeneralClass (this will be inhereted to classes that will perform the registrations)

class GeneralClass:
    def __init__(self, sequence, slice_path, AIFs_PATH):
        self.sequence = sequence
        self.slice_path = slice_path
        self.sorted_dcms_paths = find_and_sort_paths(slice_path, '.dcm')
        # Attention the following must be set accordingly to the path of the user:
        self.patient_folder = self.slice_path.split('\\')[5]
        self.AIFs_PATH = AIFs_PATH

        if sequence == 'T1':
            self.resize_flag = False
            self.inversion_times = read_inversion_times(self.sorted_dcms_paths)
            self.outcome_str = 'T1_estimated_map'
        elif sequence == 'DTI':
            self.resize_flag = False
            self.b_values, self.bVec_original, self.image_orientation_patient, self.pixel_spacing, self.slice_thickness = read_dicom_tags_DTI(self.sorted_dcms_paths)
            self.outcome_str = 'FA_map'
        elif sequence == 'DCE':
            self.resize_flag = True
            self.aif, self.times = load_txt(self.AIFs_PATH + '\\patient_' + str(self.patient_folder[-3:]) + '\\' + 'AIF__2C Filtration__Curve.txt')
            # Hack because the aif and time txts given were always an element shorter, hence the last element is duplicated:
            self.aif.append(self.aif[-1])
            self.times.append(self.times[-1])
            self.timepoint = DCE_find_timepoint(self.aif)
            self.outcome_str = 'BF_map'
        else:
            raise Exception("The sequence given is not T1, DTI, or DCE")
    																		   
    def create_3d_mhd(self, folder_path, suffix='.dcm'):
        """
        This function takes as input a "folder_path" inside which files with 
        "suffix" are expected to be found. All the "suffix" files found inside 
        will be combined in a single 3D-mhd file that will be used for Group-Wise 
        registration.
        The 3D-mhd file will have the same name as the last folder of the 
        "folder_path" and will be saved in the "folder_path".
        """
        dcm_paths = find_and_sort_paths(folder_path, suffix)
        if not dcm_paths:
            dcm_paths = find_and_sort_paths(folder_path, '.ima')
        three_d_mhd = []
        for idx, dcm_path in enumerate(dcm_paths):
            current_image = io.imread(dcm_path, plugin='simpleitk').squeeze()
            if idx==0:
                original_shape = (len(dcm_paths),) + np.shape(current_image)
            if self.resize_flag:
                # self.resize_flag is True only for DCE which we want to resize to (192, 192)
                cropped_image = cv2.resize(current_image, dsize=(192, 192))
                if idx==0:
                    three_d_mhd = np.zeros((len(dcm_paths), 192, 192))
            else:
                cropped_image = current_image
                if idx==0:
                    three_d_mhd = np.zeros(original_shape)  
            three_d_mhd[idx, :, :] = cropped_image
        sitk.WriteImage(sitk.GetImageFromArray(three_d_mhd), folder_path + '\\' + folder_path.split('\\')[-1] + '.mhd') 
        return original_shape
    
    def parallel_DCE_Fitting(self, shape, images, output_path, save_flag):
        indices = range(shape[0]*shape[1])
        images = images.reshape((shape[0]*shape[1], shape[2]))
        # construct pool
        pool = mp.Pool(mp.cpu_count())
        print('Parallel Processing for DCE fitting ...' )
        blood_flow_map = np.zeros((shape[0]*shape[1],))
        Fp = np.zeros((shape[0]*shape[1],))
        Tp = np.zeros((shape[0]*shape[1],)) 
        PS = np.zeros((shape[0]*shape[1],))
        Te = np.zeros((shape[0]*shape[1],))
        fitted = np.zeros((shape[0]*shape[1], shape[2]))
        args = [(images[idx, :], self.times, self.aif, self.timepoint) for idx in indices]
        results = pool.starmap(Linear_Least_Squares_2CFM, args)        

        for index, result in zip([idx for idx in tqdm(range(len(indices)))], results):
            blood_flow_map[index] = result[0]
            if np.shape(result[1])==(shape[2],):
                fitted[index,:] = result[1]
            # For the cases where Linear_Least_Squares_2CFM return fit=0, fit gets converted to a vector full of 0
            elif np.shape(result[1])==():
                fitted[index,:] = np.zeros((shape[2]))
            Fp[index] = result[2]
            Tp[index] = result[3]
            PS[index] = result[4]
            Te[index] = result[5]
        
        blood_flow_map = np.array(blood_flow_map)
        blood_flow_map = blood_flow_map.reshape((shape[0],shape[1]))
        
        fitted = np.array(fitted)
        fitted = fitted.reshape(shape)
        
        pool.close()
        
        if save_flag and output_path!='':
            if not os.path.isdir(output_path):
                os.mkdir(output_path)
            for i in range(1, shape[-1]+1):
                sitk.WriteImage(sitk.Cast(sitk.GetImageFromArray(fitted[:,:,i-1]),sitk.sitkFloat64), output_path + str(format(i, "03")) + '.mhd')
        
        return blood_flow_map, fitted, [Fp.reshape((shape[0],shape[1])), Tp.reshape((shape[0],shape[1])), PS.reshape((shape[0],shape[1])), Te.reshape((shape[0],shape[1]))]
    
    def fitting(self, images_to_be_fitted, fitting_output_path, initial_shape, save_flag):
        """
        images_to_be_fitted: is expected to have shape [h, w, t]
        initial_shape: is expected to be in format [h, w, t]
        """
        start = time.time()
        if self.sequence=='T1':
            # Here the outcome is the T1_estimated:
            shape = np.shape(images_to_be_fitted)
            # For the T1 parallel fitting process the 3d-images will be reshaped properly 
            # since it expects [h*w, t]
            images_to_be_fitted = np.reshape(images_to_be_fitted,(shape[0]*shape[1], shape[2]))
            start = time.time()
            fitted, outcome, params, originals = parallel_T1_fitting(images_to_be_fitted, fitting_output_path, inversion_times = self.inversion_times, initial_shape = initial_shape, save_flag = save_flag)
                
        elif self.sequence=='DTI':
            # For the DTI fitting process the images will remain [h,w,t] as ensured by the following command:
            images_to_be_fitted = np.reshape(images_to_be_fitted, initial_shape)
            # Here the outcome is the FA_map:
            fitted, params_0, params_1, outcome, _ = DTI_fitting_wrapper(images_to_be_fitted, fitting_output_path, b_values = self.b_values,bVec_original = self.bVec_original, image_orientation_patient = self.image_orientation_patient, pixel_spacing = self.pixel_spacing, slice_thickness = self.slice_thickness, save_flag=save_flag)
            params = [params_0, params_1]
            outcome = outcome.reshape((initial_shape[0], initial_shape[1]))
            
        elif self.sequence=='DCE':
            # For the DCE parallel fitting process the images remain as is hence in [h,w,t] shape:
            dce_shape = np.shape(images_to_be_fitted) # effectively: dce_shape is the same as initial_shape and should be substituted
            outcome, fitted, params = self.parallel_DCE_Fitting(dce_shape, images_to_be_fitted, fitting_output_path, save_flag)
        else:
            raise Exception("The sequence given is not T1, DTI, or DCE")
        stop = time.time()
        time_taken = stop - start
        return fitted, outcome, np.round(time_taken/60, 3), params
    
    def registration_metrics_wrapper(self):
        np.save(self.output_dir + '\\Original_Fitted_' + self.technique_str, self.fitted_original.data)
        np.save(self.output_dir + '\\Final_Fitted_' + self.technique_str, self.fitted_final.data)
        print(self.sequence)
                              
        rows_start, rows_end, cols_start, cols_end = cropping_dimensions(self.sequence)
        original_map = copy.deepcopy(self.map_original)[rows_start:rows_end, cols_start:cols_end]
        corrected_map = copy.deepcopy(self.map_corrected)[rows_start:rows_end, cols_start:cols_end]
        observed_original = copy.deepcopy(self.original_mhd)[:, rows_start:rows_end, cols_start:cols_end]
        expected_original = copy.deepcopy(self.fitted_original)[rows_start:rows_end, cols_start:cols_end, :]
        
        if self.technique_str=='MoCoMo':
            observed_method = copy.deepcopy(self.result_mhd)[rows_start:rows_end, cols_start:cols_end, :]
            observed_method = np.moveaxis(observed_method, -1, 0)
            
        elif self.technique_str=='GroupWise_Huizinga':    
            observed_method = copy.deepcopy(self.result_mhd)[:, rows_start:rows_end, cols_start:cols_end]
        
        else:
            raise Exception("Unknown technique_str")
        expected_method = copy.deepcopy(self.fitted_final)[rows_start:rows_end, cols_start:cols_end, :]  
                                    
        
        # Original Metrics
        self.original_metrics_dict = calculate_all_metrics(observed_original, sharpness(original_map))
        original_dict = copy.deepcopy(self.original_metrics_dict)
        del original_dict['coefficient_of_variation']
        del original_dict['std_sec_der']
        df_original = pd.DataFrame.from_dict([original_dict])  
        
        # Chi-Original
        observed_original = np.moveaxis(observed_original, 0 ,-1)
        self.original_metrics_dict['chisquared_original'] =  chi_squared_normalised(observed_original, expected_original)
        temp_chi_squared = copy.deepcopy(self.original_metrics_dict['chisquared_original'])
        df_original['mean_chi_squared'] = np.round(np.nanmean(temp_chi_squared), 3)
        df_original['std_chi_squared'] = np.round(np.nanstd(temp_chi_squared), 3)
        df_original['median_chi_squared'] = np.round(np.nanmedian(temp_chi_squared), 3)
        df_original['iqr_chi_squared'] = np.round(iqr(temp_chi_squared), 3)
    
    
        # Method Metrics
        self.method_metrics_dict = calculate_all_metrics(observed_method, sharpness(corrected_map))
        method_dict = copy.deepcopy(self.method_metrics_dict)
        del method_dict['coefficient_of_variation']
        del method_dict['std_sec_der']
        df_method = pd.DataFrame.from_dict([method_dict])
        
        # Chi-Method
        observed_method = np.moveaxis(observed_method, 0 ,-1)                    
        self.method_metrics_dict['chisquared_' + self.technique_str] = chi_squared_normalised(observed_method, expected_method)
        temp_chi_squared_method = copy.deepcopy(self.method_metrics_dict['chisquared_' + self.technique_str])
        df_method['mean_chi_squared'] = np.round(np.nanmean(temp_chi_squared_method), 3)
        df_method['std_chi_squared'] = np.round(np.nanstd(temp_chi_squared_method), 3)
        df_method['median_chi_squared'] = np.round(np.nanmedian(temp_chi_squared_method), 3)
        df_method['iqr_chi_squared'] = np.round(iqr(temp_chi_squared_method), 3)
        
        # Save all metrics Original and Method in csv      
        df_total = pd.concat([df_original, df_method], keys=[self.patient_folder + '_' + 'Original', self.patient_folder + '_' + self.technique_str])

        csv_path = '\\'.join(self.slice_path.split('\\')[:-4])+ '\\' + self.current_slice + '_Metrics_' + self.technique_str + '_cropped.csv'
        
        if not os.path.isfile(csv_path):
            df_total.to_csv(csv_path)
        else: # else it exists so append without writing the header
            df_total.to_csv(csv_path, mode='a', header=False)
        
        return df_original['median_chi_squared'][0], df_method['median_chi_squared'][0]
    

class Huizinga(GeneralClass):
    
    def __init__(self, sequence, slice_path, AIFs_PATH, data_path, sequence_images_path, parameters_file):
        GeneralClass.__init__(self, sequence, slice_path, AIFs_PATH)
        self.current_slice = self.slice_path.split('\\')[-1]
        self.data_path = data_path
        self.sequence_images_path = sequence_images_path
        self.parameters_file = parameters_file + 'GroupWise_Huizinga_' + sequence + '.txt'
        self.img_path = self.data_path + '\\' + self.sequence_images_path + '\\' + self.current_slice + '\\' + self.current_slice + '.mhd'
        self.output_dir = self.data_path + '\\' + self.sequence_images_path + '\\' + self.current_slice + '_GroupWise_Huizinga'
        self.technique_str = 'GroupWise_Huizinga'
        
        #Create the 3d mhd file of the original images (it returns the shape:[t, h, w] which is the shape of the 3d_mhd file it stored):
        self.initial_shape = self.create_3d_mhd(self.slice_path, suffix='.dcm')

        # Load the 3d mhd file of the original images (it loads the 3d_mhd file with shape [t, h, w]:
        self.original_mhd = load_mhd3d(self.slice_path + '\\' + self.current_slice + '.mhd')
    
        
    def Huizinga_wrapper_function(self):
        # Perform GroupWise Registration:
        start = time.time()

        # http://elastix.bigr.nl/wiki/index.php/Par0039
        # For groupwise registration one should use the following command line call:
        # elastix -f <qMRI image> -m <qMRI image> -p <par filename> -out <output dir>
        # where <qMRI image> is the entire group of images that are acquired in a single quantitative MRI acquisition. 
        # Note that the fixed and moving image should be the same.
        elastix_registration(self.img_path, self.img_path, self.output_dir, self.parameters_file)
        stop = time.time()
        self.time_huizinga_registration = np.round((stop-start)/60, 3)
        
        # Load Huizinga corrected images (with shape [t, h, w]):
        self.result_mhd = load_mhd3d(self.output_dir + '\\' + 'result.0.mhd')
        
        # Create Unregistered map:
        self.fitted_original, self.map_original, self.time_original_fitting, self.params_original = self.fitting(np.rollaxis(self.original_mhd, 0, 3), '', (self.initial_shape[1], self.initial_shape[2], self.initial_shape[0]), False)
        np.save(self.slice_path + '\\' + self.sequence+ '_Original_' + self.outcome_str + '.npy', self.map_original)
        np.save(self.slice_path + '\\' + 'Fitting_Parameters_Original.npy', self.params_original)
        
        # Create GroupWise Registered map:
        self.fitted_final, self.map_corrected, self.time_huizinga_fitting, self.params_final = self.fitting(np.rollaxis(self.result_mhd, 0, 3), '', (self.initial_shape[1], self.initial_shape[2], self.initial_shape[0]), False)
        np.save(self.output_dir + '\\' + self.sequence+ '_GroupWise_Huizinga_' + self.outcome_str+'.npy', self.map_corrected)
        np.save(self.output_dir + '\\' + 'Fitting_Parameters_Huizinga.npy', self.params_final)
        
    
class MoCoMo(GeneralClass):
    
    def __init__(self, sequence, slice_path, sequence_images_path, data_path, MoCoMo_slice_path, crop_flag=False, emergency_stop=35, tolerance=1):
        self.technique_str = 'MoCoMo'
        GeneralClass.__init__(self, sequence, slice_path, AIFs_PATH)
        self.sequence_images_path = sequence_images_path
        self.current_slice = self.slice_path.split('\\')[-1]
        self.data_path = data_path
        self.sorted_dcms_paths = find_and_sort_paths(slice_path, '.dcm')
        self.MoCoMo_slice_path = MoCoMo_slice_path
        self.crop_flag = crop_flag
        self.emergency_stop = emergency_stop
        
        if self.sequence == 'T1':
            self.tolerance = tolerance
        elif self.sequence == 'DCE':
            # for DCE we take half because we have downsampled the original image
            self.tolerance = tolerance/ 2
        elif self.sequence == 'DTI':
            self.tolerance = tolerance
        
        handle_original_mhd(self.slice_path, suffix='.dcm', crop_flag=self.crop_flag)

        # Load the original images in a 3d stack (this will be needed for the fitting process)
        # Create the 3d mhd file of the original images: (it returns the shape:[t, h, w] which is the shape of the 3d_mhd file it stored)
        self.initial_shape = self.create_3d_mhd(self.slice_path, suffix='.dcm')
        # Load the 3d mhd file of the original images: (it loads the 3d_mhd file with shape [t, h, w])
        self.original_mhd = load_mhd3d(self.slice_path + '\\' + self.current_slice + '.mhd')

    
    def bsplines_registration(self, moving_images_paths, fixed_images_paths, number_of_registration):    
        moving_images_paths = glob.glob(moving_images_paths)
        fixed_images_paths = glob.glob(fixed_images_paths)
        output_dir = self.MoCoMo_slice_path+'\\BSplines_Registered_'+ str(number_of_registration)
        if self.sequence == 'DTI':
            parameters_file = self.data_path + '\\' + 'Elastix_Parameters_Files\\' + 'BSplines_DTI.txt'
        elif self.sequence == 'T1':
            parameters_file = self.data_path + '\\' + 'Elastix_Parameters_Files\\' + 'BSplines_T1.txt'
        elif self.sequence == 'DCE':
            parameters_file = self.data_path + '\\' + 'Elastix_Parameters_Files\\' + 'BSplines_DCE.txt'
        else:
            raise Exception("Unknown sequence")
        MoCoMo_elastix_registration_wrapper(moving_images_paths, fixed_images_paths, output_dir, parameters_file)
    
    
    @staticmethod
    def MoCoMo_load_deformation_field(path):
        # The deformation field for each 2d images has shape (h, w, 2) the last is equal to 2 one for each dimension x and y
        mhds = glob.glob(path+'\*.mhd')
        assert sorted(mhds)==mhds
        df = []
        for idx, mhd in enumerate(mhds):
            if idx==0:
                current_image = io.imread(mhd, plugin='simpleitk').squeeze()
                original_shape = np.shape(current_image) + (len(mhds),)
                df = np.zeros(original_shape)
                df[:, :, :, idx] = current_image
            else:    
                df[:, :, :, idx] = io.imread(mhd, plugin='simpleitk').squeeze()
        return df
    
    def MoCoMo_wrapper_function(self):
        df_eucl_distances = []
        self.MoCoMo_time_total = 0
        
        # First Fitting:
        # self.fitting expects: [h, w, t]
        self.fitted_original, self.map_original, self.time_original_fitting, self.params_original = self.fitting(np.rollaxis(self.original_mhd, 0, 3), 
                                                                                                                 self.MoCoMo_slice_path + '\\Fitted\\', 
                                                                                                                 (self.initial_shape[1], self.initial_shape[2], self.initial_shape[0]),
                                                                                                                 True)
        np.save(self.MoCoMo_slice_path + '\\' + 'Original_map.npy' , self.map_original)
        self.fitted_original.dump(self.MoCoMo_slice_path + '\\' + '0_Fitted')
        np.save(self.MoCoMo_slice_path + '\\' + '0_Fitted_Parameters.npy', self.params_original)
    
        # self.MoCoMo_time_total will be measured in seconds, while self.time_original_fitting is in minutes:
        self.MoCoMo_time_total = self.MoCoMo_time_total + self.time_original_fitting*60
    
        np.save(self.MoCoMo_slice_path + '\\Original\\' + 'Fitting_Parameters_Original.npy', self.params_original)
                    
            
        for i in range(self.emergency_stop):
            
            iterations = i + 1
            
            if not os.path.isdir(self.MoCoMo_slice_path+'\\BSplines_Registered_' + str(iterations)):
                os.mkdir(self.MoCoMo_slice_path+'\\BSplines_Registered_' + str(iterations))
            
            start = time.time()
            

            #BSplines Registration with Elastix:
            self.bsplines_registration(self.MoCoMo_slice_path + '\\Original\\'+'\*mhd',
                                       self.MoCoMo_slice_path + '\\Fitted\\'+'\*mhd',
                                       iterations)
            
            stop = time.time()
            # self.MoCoMo_time_total will be measured in seconds
            self.MoCoMo_time_total = self.MoCoMo_time_total + (stop - start)
            
            self.output_dir = self.MoCoMo_slice_path+'\\BSplines_Registered_' + str(iterations)
            previous_output_dir = self.MoCoMo_slice_path+'\\BSplines_Registered_' + str(iterations-1)
            
            self.result_mhd = find_and_load_as_mhd(self.output_dir, suffix='.mhd', crop_flag = self.crop_flag)
            self.fitted_final, self.map_corrected, self.time_fitting, self.params_final = self.fitting(self.result_mhd, 
                                                                                                       self.MoCoMo_slice_path + '\\Fitted\\', 
                                                                                                       (self.initial_shape[1], self.initial_shape[2], self.initial_shape[0]), 
                                                                                                       True)
            
            np.save(self.output_dir + 'Fitted_MoCoMo.npy', self.fitted_final.data)
            # self.MoCoMo_time_total will be measured in seconds, while self.time_original_fitting is in minutes:
            self.MoCoMo_time_total = self.MoCoMo_time_total + self.time_fitting*60
            
            np.save(self.output_dir + 'Fitting_Parameters_MoCoMo.npy', self.params_final)
            
            # Create deformation field for the BSplines registered images:
            self.MoCoMo_create_deformation_field()
                        
            # Decide if the next iteration is going to be calculated:
            if i==0:
                diction = OrderedDict()
                diction['largest_deformation'] = []
                median_improvement_list=[]
                deformation_csv_path = self.output_dir + 'largest_deformations.csv' 
                pass
            else:
                df_new = self.MoCoMo_load_deformation_field(self.output_dir + '\\deformation_field\\')
                df_previous = self.MoCoMo_load_deformation_field(previous_output_dir + '\\deformation_field\\')
                df_difference = df_new - df_previous
                df_difference_x_squared = np.square(df_difference[:,:,0,:].squeeze())
                df_difference_y_squared = np.square(df_difference[:,:,1,:].squeeze())
                diff_norm = np.sqrt(np.add(df_difference_x_squared, df_difference_y_squared))
                df_eucl_distances.append(diff_norm)
                
                maximum_deformations_per_pixel  = np.nanmax(diff_norm, axis=2)
                median_imporvement = np.median(maximum_deformations_per_pixel)
                median_improvement_list.append(median_imporvement)
                        
                print('MEDIAN OF LARGEST CHANGES: ', median_imporvement)
                print('\n\n')
                
                if median_imporvement <= self.tolerance:
                    np.save(self.MoCoMo_slice_path+'\\BSplines_Registered_' + str(iterations) + '\\'+'MoCoMo_'+self.outcome_str+'.npy', self.map_corrected)
                    break
                    
        diction['largest_deformation'] = median_improvement_list
        df_deformations = pd.DataFrame.from_dict([diction])
        df_deformations.to_csv(deformation_csv_path, index=False)
        
        return df_eucl_distances


    def MoCoMo_create_deformation_field(self):
        start = time.time()
        path = self.output_dir
        if not os.path.isdir(path+'\\deformation_field'):
            os.mkdir(path+'\\deformation_field')
        txts = glob.glob(path+'\*.txt')
        transform_parameters_files = []
        for txt in txts:
            if 'TransformParameters' in txt:
                transform_parameters_files.append(txt)
        assert sorted(transform_parameters_files)==transform_parameters_files
        for transform_parameter_file in transform_parameters_files:
            transformix_deformation_field(path+'\\deformation_field\\', transform_parameter_file)
            change_elastix_naming_result(path+'\\deformation_field\\deformationField.mhd', 'deformationField.raw', 'deformationField_' + transform_parameter_file.split('.txt')[0][-3:] + '.raw')
            os.rename(path+'\\deformation_field\\deformationField.mhd', path+'\\deformation_field\\deformationField_' + transform_parameter_file.split('.txt')[0][-3:] + '.mhd')
            os.rename(path+'\\deformation_field\\deformationField.raw', path+'\\deformation_field\\deformationField_' + transform_parameter_file.split('.txt')[0][-3:] + '.raw')
        stop = time.time()
        return start - stop
        
##############################################################################    
### INPUTS:
# The user will provide the full path to a folder containing the .zip files
# for all th patients that will be processed:
    
#TODO: Uncomment the following
#DATA_PATH = input("Please provide the full path to folder containg the .zip "
#                  "files for the patients that will be motion corrected "
#                  " Note that if DCE is to be processed as well it is exepected"
#                  " to find a folder with format AIF -> patient ### -> AIF__2C Filtration__Curve.txt"
#                  " that contains the AIFs for each patient, the parent folder AIF is "
#                  " expected to be found inside DATA_PATH")
    
#assert os.path.exists(DATA_PATH), "Please provide a valid path for the .zip folders"

# The user will choose to perform either Group-wise (Huizinga et al)
# or MoCoMo Registration     
#TODO: Uncomment the following
#TECHNIQUE = input("Please insert the number 1 for Group-wise Registration "
#                  "or number 2 for MoCoMo Registration ")
#assert TECHNIQUE in ['1', '2'], "Please insert a proper value for the motion correction technique that will be performed."

if __name__ == '__main__':
    # The sequences to be motion corrected
    POOL_OF_SEQUENCES = ['T1', 'DTI', 'DCE']
    SEQUENCES = ['T1', 'DTI','DCE']#  
    #assert set(SEQUENCES).issubset(set(POOL_OF_SEQUENCES)), " Please provide a proper list of sequences choosing from: %s" % POOL_OF_SEQUENCES
    
    # CORRESPONDANCE is a dictionary which contains as keys the names
    # of the sequences that will processed and as values a list:
    # the first element of the list is the numbre of the folder that contains the 
    # respective sequences
    # the second element of the list is the numbr of files expected to be found in
    # the folder
    # the third element is the number of slices for each sequence
    CORRESPONDANCE = OrderedDict()
    CORRESPONDANCE['T1'] = [19, 140, 5] 
    CORRESPONDANCE['DTI'] = [31, 4380, 30]    
    CORRESPONDANCE['DCE'] = [39, 2385, 9]
    
    ### MAIN:
    # Unzip to extract Leeds_Patients_ ... folders
    #DATA_PATH = r'C:\Users\md1jgra\Desktop\FotisMotionCorrectionPaper\ToRun'
    DATA_PATH = os.getcwd() 
    AIFs_PATH = DATA_PATH + r'\AIFs'
    TECHNIQUE = 1
    
    
    if TECHNIQUE == 1:
        technique_str = 'GroupWise_Huizinga'
    elif TECHNIQUE == 2:
        technique_str = 'MoCoMo'
    #unzip_patients(DATA_PATH)
    
    # Reorganizing files per sequence:
    for sequence in SEQUENCES:
        folder = CORRESPONDANCE[sequence][0]
        num_of_files = CORRESPONDANCE[sequence][1]
        slices = CORRESPONDANCE[sequence][2]
        
        os.chdir(DATA_PATH)
        print(DATA_PATH)
        patients_folders = os.listdir()
        for patient_folder in patients_folders:
            if patient_folder not in ['Leeds_Patient_4128001', 'Leeds_Patient_4128002','Leeds_Patient_4128003', 'Leeds_Patient_4128004', 'Leeds_Patient_4128005', 'Leeds_Patient_4128006', 'Leeds_Patient_4128007','Leeds_Patient_4128008','Leeds_Patient_4128009', 'Leeds_Patient_4128010']:
                continue
            #TODO: Scratch the following line completely
            #patient_folder = 'Leeds_Patient_4128010'
            # Change directory to find the sequence images:
            sequence_images_path = patient_folder + '\\' + str(folder) + '\\DICOM'
            print(sequence_images_path)
            os.chdir(DATA_PATH + '\\' + sequence_images_path)
            
            # Make sure that no acquisition is missing:
            dcm_files_found = glob.glob("*.dcm")
            if not dcm_files_found:
                dcm_files_found = glob.glob("*.ima")
            #TODO: Uncomment this
            #assert len(dcm_files_found)==num_of_files#, "Did not find the expect number of images in %s " %os.getcwd()
            
            # Create slice_folders for the images to be arranged: 
            create_slice_folders(sequence, slices)
            
            #Re-arrange images in their respective slice_folders:
            if dcm_files_found:
                move_images_to_slice_folders(sequence, dcm_files_found, slices, num_of_files, folder)
            
            for slice in range(1, slices+1):
                current_slice = sequence + '_slice_' + str(slice)
                #TODO: Delete this to work for all slices and sequences
                if sequence == 'T1' and current_slice not in [sequence + '_slice_3']:
                    continue
                elif sequence == 'DCE' and current_slice not in [sequence + '_slice_5']:
                    continue
                elif sequence == 'DTI' and current_slice not in [sequence + '_slice_15']:
                    continue
                
                slice_path = DATA_PATH + '\\' + sequence_images_path + '\\' + current_slice
                
                if TECHNIQUE == 1:
                    
                    HuiReg = Huizinga(sequence, 
                                      slice_path, 
                                      AIFs_PATH,
                                      DATA_PATH,
                                      sequence_images_path,
                                      parameters_file = str(os.path.join(DATA_PATH, 'Elastix_Parameters_Files') + '\\'))
                    HuiReg.Huizinga_wrapper_function()
                    HuiReg.registration_metrics_wrapper()
                                       
                   
                elif TECHNIQUE ==2:
                     # Create needed folders:
                    MoCoMo_slice_path = DATA_PATH + '\\' + sequence_images_path + '\\' + current_slice + '_' + technique_str
                    create_MoCoMo_folders(MoCoMo_slice_path)
                    
                    # Perform MoCoMo:
                    # tolerance is set to 1 pixel
                    myReg = MoCoMo(sequence, slice_path, sequence_images_path, DATA_PATH, MoCoMo_slice_path, crop_flag = False, emergency_stop = 35, tolerance = 1)
                    df_eucl_distances = myReg.MoCoMo_wrapper_function()
                    myReg.registration_metrics_wrapper()
                                   
            #break            