# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:02:00 2019

@author: wsy88
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 12:17:33 2018

@author: vaskalab
This code is for analyzing reconstructed PET images of different PSF kernels, different iteration #, different noise realizations for 
1) spatial variability (background noise)
which is average of std of background of multiple realizations 
2) emsemble noise 
which is the std of the mean of background ROI of different noise relizations
3) contrast reconvery 
which is the (average of selected lesion intensity / average of the bkg ROI - 1) / (true contrast which is 9 in our case)


"""


import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from collections import defaultdict

string_a1 = 'E:\\newnew_article_results\\newresults_inproposal\\8spheres_truePSF_1mmvoxel_new_'
string_a2 = '\\truePSF_1mmvoxel_'
#string_b1 = 'E:\\newnew_article_results\\newresults_inproposal\\8spheres_kernel1_1mmvoxel_new_'
#string_b2 = '\\kernel1_1mmvoxel_'
#string_c1 = 'E:\\newnew_article_results\\newresults_inproposal\\8spheres_kernel2_1mmvoxel_new_'
#string_c2 = '\\kernel2_1mmvoxel_'
#string_d1 = 'E:\\newnew_article_results\\newresults_inproposal\\8spheres_kernel3_1mmvoxel_new_'
#string_d2 = '\\kernel3_1mmvoxel_'
#string_e1 = 'E:\\newnew_article_results\\newresults_inproposal\\8spheres_kernel4_1mmvoxel_new_'
#string_e2 = '\\kernel4_1mmvoxel_'
string_f1 = 'E:\\newnew_article_results\\newresults_inproposal\\8spheres_mlem_1mmvoxel_new_'
string_f2 = '\\nopsf_1mmvoxel_'


#string_a1 = 'E:\\newnew_article_results\\PSF-MLEM_newphantom_lowercontrasttobkgratio\\trueLSOdensity\\12spheres_'
#string_a2 = '\\12spheres_'
#string_b1 = 'E:\\newnew_article_results\\MLEM_newphantom_lowercontrasttobkgratio\\12spheres_'
#string_b2 = '\\12spheres_'
#string_c1 = 'E:\\newnew_article_results\\PSF-MLEM_newphantom_lowercontrasttobkgratio\\5timesLSOdensity\\12spheres_'
#string_c2 = '\\12spheres_'

string_list = []
std_mean_iteration_20mm_2mm = []
mean_bias_iteration_20mm_2mm = []
std_mean_ROI_20mm_2mm = []
image_roughness_iteration_20mm_2mm = []

std_mean_iteration_20mm_3mm = []
mean_bias_iteration_20mm_3mm = []
std_mean_ROI_20mm_3mm = []
image_roughness_iteration_20mm_3mm = []

std_mean_iteration_20mm_5mm = []
mean_bias_iteration_20mm_5mm = []
std_mean_ROI_20mm_5mm = []
image_roughness_iteration_20mm_5mm = []

std_mean_iteration_20mm_10mm = []
mean_bias_iteration_20mm_10mm = []
std_mean_ROI_20mm_10mm = []
image_roughness_iteration_20mm_10mm = []

std_mean_iteration_40mm_2mm = []
mean_bias_iteration_40mm_2mm = []
std_mean_ROI_40mm_2mm = []
image_roughness_iteration_40mm_2mm = []

std_mean_iteration_40mm_3mm = []
mean_bias_iteration_40mm_3mm = []
std_mean_ROI_40mm_3mm = []
image_roughness_iteration_40mm_3mm = []

std_mean_iteration_40mm_5mm = []
mean_bias_iteration_40mm_5mm = []
std_mean_ROI_40mm_5mm = []
image_roughness_iteration_40mm_5mm = []

std_mean_iteration_40mm_10mm = []
mean_bias_iteration_40mm_10mm = []
std_mean_ROI_40mm_10mm = []
image_roughness_iteration_40mm_10mm = []


#avg_std_ROI_iteration = []
contrast_recovery_iteration = []
#average_image_iteration = []

string_list.append((string_a1, string_a2))
#string_list.append((string_b1, string_b2))
#string_list.append((string_c1, string_c2))
#string_list.append((string_d1, string_d2))
#string_list.append((string_e1, string_e2))
string_list.append((string_f1, string_f2))
#methods = ['truePSF', 'no-PSF'] #
methods = ['truePSF','kernel1','kernel2', 'kernel3', 'kernel4','no-PSF']

def sphericalROI(img, center, slice_range, radius):
    return [img[z, x, y]  for x in range(0,153) for y in range(0,153) for z in slice_range if math.sqrt(((x+0.5)*1.0-center[0])**2+((y+0.5)*1.0-center[1])**2 + ((z+0.5)*1.0-center[2])**2)  <= radius]    

#def sphereROI(img, center, slice_range, radius):
 #   return [img[z, x, y]  for x in range(0,153) for y in range(0,153) for z in slice_range if math.sqrt(((x+0.5)*1.0-center[0])**2+((y+0.5)*1.0-center[1])**2 +((z+0.5)*1.0-center[2])**2) <= radius]    

def circularROI(img, center, slice, radius):
    return [img[slice, x, y] for x in range(0,153) for y in range(0,153) if math.sqrt(((x+0.5)*1.0-center[0])**2+((y+0.5)*1.0-center[1])**2) <= radius]    

def circularMask(img, center,  radius):
    for x in range(0,153):
        for y in range(0,153):
            if math.sqrt(((x+0.5)*1.0-center[0])**2+((y+0.5)*1.0-center[1])**2) <= radius:
                img[x,y] = 1               
    return img

for iteration in [10,50,90,130,170,210,250]: #range(10,210,10):  [10,30,50,70,90,110,130,150,170,190,210,230,250]
    filename_array = []
    #filename_array_b = []
    img_array = []
    #mean_ROI_array = []
    std_ROI_array_20mm_2mm = []
    mean_std_array_20mm_2mm = []
    mean_ROI_array_20mm_2mm = []
    bias_ROI_array_20mm_2mm = []
    std_mean_mean_ROI_20mm_2mm = []
    mean_mean_ROI_20mm_2mm = []
    std_mean_ROI_20mm_2mm = []
    
    std_ROI_array_20mm_3mm = []
    mean_std_array_20mm_3mm = []
    mean_ROI_array_20mm_3mm = []
    bias_ROI_array_20mm_3mm = []
    std_mean_mean_ROI_20mm_3mm = []
    mean_mean_ROI_20mm_3mm = []
    std_mean_ROI_20mm_3mm = []
    
    std_ROI_array_20mm_5mm = []
    mean_std_array_20mm_5mm = []
    mean_ROI_array_20mm_5mm = []
    bias_ROI_array_20mm_5mm = []
    std_mean_mean_ROI_20mm_5mm = []
    mean_mean_ROI_20mm_5mm = []
    std_mean_ROI_20mm_5mm = []
 
    std_ROI_array_20mm_10mm = []
    mean_std_array_20mm_10mm = []
    mean_ROI_array_20mm_10mm = []
    bias_ROI_array_20mm_10mm = []
    std_mean_mean_ROI_20mm_10mm = []
    mean_mean_ROI_20mm_10mm = []
    std_mean_ROI_20mm_10mm = []
      
    std_ROI_array_40mm_2mm = []
    mean_std_array_40mm_2mm = []
    mean_ROI_array_40mm_2mm = []
    bias_ROI_array_40mm_2mm = []
    std_mean_mean_ROI_40mm_2mm = []
    mean_mean_ROI_40mm_2mm = []
    std_mean_ROI_40mm_2mm = []
    
    std_ROI_array_40mm_3mm = []
    mean_std_array_40mm_3mm = []
    mean_ROI_array_40mm_3mm = []
    bias_ROI_array_40mm_3mm = []
    std_mean_mean_ROI_40mm_3mm = []
    mean_mean_ROI_40mm_3mm = []
    std_mean_ROI_40mm_3mm = []
    
    std_ROI_array_40mm_5mm = []
    mean_std_array_40mm_5mm = []
    mean_ROI_array_40mm_5mm = []
    bias_ROI_array_40mm_5mm = []
    std_mean_mean_ROI_40mm_5mm = []
    mean_mean_ROI_40mm_5mm = []
    std_mean_ROI_40mm_5mm = []
 
    std_ROI_array_40mm_10mm = []
    mean_std_array_40mm_10mm = []
    mean_ROI_array_40mm_10mm = []
    bias_ROI_array_40mm_10mm = []
    std_mean_mean_ROI_40mm_10mm = []
    mean_mean_ROI_40mm_10mm = []
    std_mean_ROI_40mm_10mm = []
    
       
    contrast_recovery_methods = []
    
    #fig = plt.figure(figsize=(15,4))
    #ax0 = fig.add_subplot(121)
    #ax1 = fig.add_subplot(122)
    for string_1, string_2 in string_list:
        filename_subarray = []
        for realization in range(1,11):
            filename_subarray.append(string_1+str(realization)+string_2+str(iteration)+'.v')
        filename_array.append(filename_subarray)
            #filename_array_b.append(string_b+str(realization)+'\\newcontrast12spheres_0.7mmDOI_2timesLSOdensity_'+str(realization)+'_'+str(iteration)+'.v')
    for indexing, filename_subarraynew in enumerate(filename_array):
        img_subarray = []
        mean_bkgROI = []
        mean_ROI_20mm_2mm = []
        bias_ROI_20mm_2mm = []       
        std_ROI_20mm_2mm = []
        std_mean_ROI_20mm_2mm = []
        
        mean_ROI_20mm_3mm = []
        bias_ROI_20mm_3mm = []       
        std_ROI_20mm_3mm = []
        std_mean_ROI_20mm_3mm = []

        mean_ROI_20mm_5mm = []
        bias_ROI_20mm_5mm = []       
        std_ROI_20mm_5mm = []
        std_mean_ROI_20mm_5mm = []

        
        mean_ROI_20mm_10mm = []
        bias_ROI_20mm_10mm = []       
        std_ROI_20mm_10mm = []
        std_mean_ROI_20mm_10mm = []
        
        mean_ROI_40mm_2mm = []
        bias_ROI_40mm_2mm = []       
        std_ROI_40mm_2mm = []
        std_mean_ROI_40mm_2mm = []
        
        
        mean_ROI_40mm_3mm = []
        bias_ROI_40mm_3mm = []       
        std_ROI_40mm_3mm = []
        std_mean_ROI_40mm_3mm = []
        
        mean_ROI_40mm_5mm = []
        bias_ROI_40mm_5mm = []       
        std_ROI_40mm_5mm = []
        std_mean_ROI_40mm_5mm = []
        
        mean_ROI_40mm_10mm = []
        bias_ROI_40mm_10mm = []       
        std_ROI_40mm_10mm = []
        std_mean_ROI_40mm_10mm = []
        bkgROI_realization = []
        
        
        
        contrast_recovery_realization = []
        for index, filename in enumerate(filename_subarraynew):
            print(filename)
           # if filename == 'E:\\newnew_article_results\\PSF-MLEM-newphantom\\newcontrast_phantom_mlem\\newcontrast_phantom_mlem_10\\newcontrast12spheres_mlem_10_250.v':
            #    stop = 1
                
            f = open(filename,'rb')
            data = np.fromfile(f,dtype=np.float32)
            img = 26.57*np.reshape(data,[47,153,153])   #used to be 19.58
            bkgROI = []
            #ROI = img[5:15,55:75,55:75]  #this is just simple rectangular ROI
            #what if I want to have circular ROI?
            #ROI = circularROI(img,[51.2,66.2], range(19,20), 10)  
            for slice in [10,15,30,35]:   #48 background ROIs
                bkgROI.append(circularROI(img,[61.5,76.5],slice,5))
                bkgROI.append(circularROI(img,[46.5,76.5],slice,5))
                bkgROI.append(circularROI(img,[31.5,76.5],slice,5))
                bkgROI.append(circularROI(img,[76.5,91.5],slice,5))
                bkgROI.append(circularROI(img,[76.5,106.5],slice,5))
                bkgROI.append(circularROI(img,[76.5,121.5],slice,5))
                bkgROI.append(circularROI(img,[91.5,76.5],slice,5))
                bkgROI.append(circularROI(img,[106.5,76.5],slice,5))
                bkgROI.append(circularROI(img,[121.5,76.5],slice,5))
                bkgROI.append(circularROI(img,[76.5,61.5],slice,5))
                bkgROI.append(circularROI(img,[76.5,46.5],slice,5))
                bkgROI.append(circularROI(img,[76.5,31.5],slice,5))
           
            bkgROI_realization.append(bkgROI)
            bkgROI_sum = sum([sum(bkgROI_elem) for bkgROI_elem in bkgROI])
            bkgROI_len = sum([len(bkgROI_elem) for bkgROI_elem in bkgROI])
            mean_bkgROI_value = bkgROI_sum/bkgROI_len
            
            
           # bkgROI = circularROI(img,[76.5,76.5,12.5],range(5,15), 30)
            mean_bkgROI.append(mean_bkgROI_value)
            # ROI = circularROI(img,[61.2,76.2],range(19,20),1.5)
           
            img_subarray.append(img[23,:,:])
            #Basically I am drawing the ROI based on the real size of the lesions
            """
            ROI_15mm_3mm = circularROI(img,[61.2,76.2,17.4],range(12,18),1.5)
            ROI_15mm_4mm = circularROI(img,[76.2,91.2,17.4],range(12,18),2)
            ROI_15mm_5mm = circularROI(img,[91.2,76.2,17.4],range(12,18),2.5)
            ROI_15mm_6mm = circularROI(img,[76.2,61.2,17.4],range(12,18),3)
            ROI_30mm_3mm = circularROI(img,[46.2,76.2,17.4],range(12,18),1.5)
            ROI_30mm_4mm = circularROI(img,[76.2,106.2,17.4],range(12,18),2)
            ROI_30mm_5mm = circularROI(img,[106.2,76.2,17.4],range(12,18),2.5)
            ROI_30mm_6mm = circularROI(img,[76.2,46.2,17.4],range(12,18),3)
            ROI_45mm_3mm = circularROI(img,[31.2,76.2,17.4],range(12,18),1.5)
            ROI_45mm_4mm = circularROI(img,[76.2,121.2,17.4],range(12,18),2)
            ROI_45mm_5mm = circularROI(img,[121.2,76.2,17.4],range(12,18),2.5)
            ROI_45mm_6mm = circularROI(img,[76.2,31.2,17.4],range(12,18),3)
            
            
            """
           #mask = circularMask(mask,[96.5,76.5],2.5) #2mm
           #mask = circularMask(mask,[56.5,76.5],2.5) #3mm
           #mask = circularMask(mask,[76.5,56.5],2.5) #10mm
           #mask = circularMask(mask,[76.5,96.5],2.5) #5mm
            
            
            
            ROI_20mm_2mm = sphericalROI(img,[96.5,76.5,23.5],range(18,30),1)   #used to be 23,24
            ROI_20mm_3mm = sphericalROI(img,[56.5,76.5,23.5],range(18,30),1.5)
            ROI_20mm_5mm = sphericalROI(img,[76.5,96.5,23.5],range(18,30),2.5)
            ROI_20mm_10mm = sphericalROI(img,[76.5,56.5,23.5],range(18,30),5)
            ROI_40mm_2mm = sphericalROI(img,[116.5,76.5,23.5],range(18,30),1)
            ROI_40mm_3mm = sphericalROI(img,[36.5,76.5,23.5],range(18,30),1.5)
            ROI_40mm_5mm = sphericalROI(img,[76.5,116.5,23.5],range(18,30),2.5)
            ROI_40mm_10mm = sphericalROI(img,[76.5,36.5,23.5],range(18,30),5)

          
            
            std_ROI_20mm_2mm.append(np.std(ROI_20mm_2mm))
            mean_ROI_20mm_2mm.append(np.mean(ROI_20mm_2mm))
            std_mean_ROI_20mm_2mm.append(np.std(ROI_20mm_2mm)/np.mean(ROI_20mm_2mm))
            bias_ROI_20mm_2mm.append((-8*3.7+np.mean(ROI_20mm_2mm))/29.6*100)
            
            std_ROI_20mm_3mm.append(np.std(ROI_20mm_3mm))
            mean_ROI_20mm_3mm.append(np.mean(ROI_20mm_3mm))
            std_mean_ROI_20mm_3mm.append(np.std(ROI_20mm_3mm)/np.mean(ROI_20mm_3mm))
            bias_ROI_20mm_3mm.append((-8*3.7+np.mean(ROI_20mm_3mm))/29.6*100)
            
            std_ROI_20mm_5mm.append(np.std(ROI_20mm_5mm))
            mean_ROI_20mm_5mm.append(np.mean(ROI_20mm_5mm))
            std_mean_ROI_20mm_5mm.append(np.std(ROI_20mm_5mm)/np.mean(ROI_20mm_5mm))
            bias_ROI_20mm_5mm.append((-8*3.7+np.mean(ROI_20mm_5mm))/29.6*100)
            
            std_ROI_20mm_10mm.append(np.std(ROI_20mm_10mm))
            mean_ROI_20mm_10mm.append(np.mean(ROI_20mm_10mm))
            std_mean_ROI_20mm_10mm.append(np.std(ROI_20mm_10mm)/np.mean(ROI_20mm_10mm))
            bias_ROI_20mm_10mm.append((-8*3.7+np.mean(ROI_20mm_10mm))/29.6*100)
            
            std_ROI_40mm_2mm.append(np.std(ROI_40mm_2mm))
            mean_ROI_40mm_2mm.append(np.mean(ROI_40mm_2mm))
            std_mean_ROI_40mm_2mm.append(np.std(ROI_40mm_2mm)/np.mean(ROI_40mm_2mm))
            bias_ROI_40mm_2mm.append((-8*3.7+np.mean(ROI_40mm_2mm))/29.6*100)
            
            std_ROI_40mm_3mm.append(np.std(ROI_40mm_3mm))
            mean_ROI_40mm_3mm.append(np.mean(ROI_40mm_3mm))
            std_mean_ROI_40mm_3mm.append(np.std(ROI_40mm_3mm)/np.mean(ROI_40mm_3mm))
            bias_ROI_40mm_3mm.append((-8*3.7+np.mean(ROI_40mm_3mm))/29.6*100)
            
            std_ROI_40mm_5mm.append(np.std(ROI_40mm_5mm))
            mean_ROI_40mm_5mm.append(np.mean(ROI_40mm_5mm))
            std_mean_ROI_40mm_5mm.append(np.std(ROI_40mm_5mm)/np.mean(ROI_40mm_5mm))
            bias_ROI_40mm_5mm.append((-8*3.7+np.mean(ROI_40mm_5mm))/29.6*100)
            
            std_ROI_40mm_10mm.append(np.std(ROI_40mm_10mm))
            mean_ROI_40mm_10mm.append(np.mean(ROI_40mm_10mm))
            std_mean_ROI_40mm_10mm.append(np.std(ROI_40mm_10mm)/np.mean(ROI_40mm_10mm))
            bias_ROI_40mm_10mm.append((-8*3.7+np.mean(ROI_40mm_10mm))/29.6*100)
            
            
            
            contrast_recovery_20mm_2mm = (np.mean(ROI_20mm_2mm)/np.mean(bkgROI) - 1)/7
            contrast_recovery_20mm_3mm = (np.mean(ROI_20mm_3mm)/np.mean(bkgROI) - 1)/7
            contrast_recovery_20mm_5mm = (np.mean(ROI_20mm_5mm)/np.mean(bkgROI) - 1)/7
            contrast_recovery_20mm_10mm = (np.mean(ROI_20mm_10mm)/np.mean(bkgROI) - 1)/7            
            contrast_recovery_40mm_2mm = (np.mean(ROI_40mm_2mm)/np.mean(bkgROI) - 1)/7
            contrast_recovery_40mm_3mm = (np.mean(ROI_40mm_3mm)/np.mean(bkgROI) - 1)/7
            contrast_recovery_40mm_5mm = (np.mean(ROI_40mm_5mm)/np.mean(bkgROI) - 1)/7
            contrast_recovery_40mm_10mm = (np.mean(ROI_40mm_10mm)/np.mean(bkgROI) - 1)/7 
            contrast_recovery_realization.append([contrast_recovery_20mm_2mm, contrast_recovery_20mm_3mm, contrast_recovery_20mm_5mm, contrast_recovery_20mm_10mm,
                                                  contrast_recovery_40mm_2mm, contrast_recovery_40mm_3mm, contrast_recovery_40mm_5mm, contrast_recovery_40mm_10mm])
    
        contrast_recovery_methods.append(contrast_recovery_realization)  #for different realizations
        img_array.append(img_subarray)
        mean_ROI_array_20mm_2mm.append(mean_ROI_20mm_2mm)
        bias_ROI_array_20mm_2mm.append(bias_ROI_20mm_2mm)
        std_ROI_array_20mm_2mm.append(std_ROI_20mm_2mm)
        mean_std_array_20mm_2mm.append(std_mean_ROI_20mm_2mm)
        
        
        mean_ROI_array_20mm_3mm.append(mean_ROI_20mm_3mm)
        bias_ROI_array_20mm_3mm.append(bias_ROI_20mm_3mm)
        std_ROI_array_20mm_3mm.append(std_ROI_20mm_3mm)
        mean_std_array_20mm_3mm.append(std_mean_ROI_20mm_3mm)

        
        mean_ROI_array_20mm_5mm.append(mean_ROI_20mm_5mm)
        bias_ROI_array_20mm_5mm.append(bias_ROI_20mm_5mm)
        std_ROI_array_20mm_5mm.append(std_ROI_20mm_5mm)
        mean_std_array_20mm_5mm.append(std_mean_ROI_20mm_5mm)

        
        mean_ROI_array_20mm_10mm.append(mean_ROI_20mm_10mm)
        bias_ROI_array_20mm_10mm.append(bias_ROI_20mm_10mm)
        std_ROI_array_20mm_10mm.append(std_ROI_20mm_10mm)
        mean_std_array_20mm_10mm.append(std_mean_ROI_20mm_10mm)
        
        mean_ROI_array_40mm_2mm.append(mean_ROI_40mm_2mm)
        bias_ROI_array_40mm_2mm.append(bias_ROI_40mm_2mm)
        std_ROI_array_40mm_2mm.append(std_ROI_40mm_2mm)
        mean_std_array_40mm_2mm.append(std_mean_ROI_40mm_2mm)
        
        mean_ROI_array_40mm_3mm.append(mean_ROI_40mm_3mm)
        bias_ROI_array_40mm_3mm.append(bias_ROI_40mm_3mm)
        std_ROI_array_40mm_3mm.append(std_ROI_40mm_3mm)
        mean_std_array_40mm_3mm.append(std_mean_ROI_40mm_3mm)
        
        mean_ROI_array_40mm_5mm.append(mean_ROI_40mm_5mm)
        bias_ROI_array_40mm_5mm.append(bias_ROI_40mm_5mm)
        std_ROI_array_40mm_5mm.append(std_ROI_40mm_5mm)
        mean_std_array_40mm_5mm.append(std_mean_ROI_40mm_5mm)
        
        mean_ROI_array_40mm_10mm.append(mean_ROI_40mm_10mm)
        bias_ROI_array_40mm_10mm.append(bias_ROI_40mm_10mm)
        std_ROI_array_40mm_10mm.append(std_ROI_40mm_10mm)
        mean_std_array_40mm_10mm.append(std_mean_ROI_40mm_10mm)
        
        #mean_std_array.append([mean_ROI_array, std_ROI_array]) 
        
    contrast_recovery_iteration.append(contrast_recovery_methods)
    if iteration == 250:
       img_avg = [np.mean(image_subarray,axis=0) for image_subarray in img_array] #average image for each method for a fixed iteration
       
    std_mean_ROI_20mm_2mm = [np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_2mm]
    mean_mean_ROI_20mm_2mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_2mm] #different methods mean fixed iteration
    mean_bias_ROI_20mm_2mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_20mm_2mm]  
    std_mean_mean_ROI_20mm_2mm = [std_mean_ROI_20mm_2mm, mean_mean_ROI_20mm_2mm]
    std_mean_iteration_20mm_2mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_20mm_2mm)])  
    mean_bias_iteration_20mm_2mm.append(mean_bias_ROI_20mm_2mm)
    image_roughness_iteration_20mm_2mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_20mm_2mm])
    
    
    std_mean_ROI_20mm_3mm=[np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_3mm]
    mean_mean_ROI_20mm_3mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_3mm] #different methods mean fixed iteration
    mean_bias_ROI_20mm_3mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_20mm_3mm]
    std_mean_mean_ROI_20mm_3mm = [std_mean_ROI_20mm_3mm, mean_mean_ROI_20mm_3mm]
    std_mean_iteration_20mm_3mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_20mm_3mm)])  
    mean_bias_iteration_20mm_3mm.append(mean_bias_ROI_20mm_3mm)
    image_roughness_iteration_20mm_3mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_20mm_3mm])

    
    std_mean_ROI_20mm_5mm = [np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_5mm]
    mean_mean_ROI_20mm_5mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_5mm] #different methods mean fixed iteration
    mean_bias_ROI_20mm_5mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_20mm_5mm]
    std_mean_mean_ROI_20mm_5mm = [std_mean_ROI_20mm_5mm, mean_mean_ROI_20mm_5mm]
    std_mean_iteration_20mm_5mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_20mm_5mm)])  
    mean_bias_iteration_20mm_5mm.append(mean_bias_ROI_20mm_5mm)
    image_roughness_iteration_20mm_5mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_20mm_5mm])
   
    std_mean_ROI_20mm_10mm=[np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_10mm]
    mean_mean_ROI_20mm_10mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_20mm_10mm] #different methods mean fixed iteration
    mean_bias_ROI_20mm_10mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_20mm_10mm]
    std_mean_mean_ROI_20mm_10mm = [std_mean_ROI_20mm_10mm, mean_mean_ROI_20mm_10mm]
    std_mean_iteration_20mm_10mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_20mm_10mm)])  
    mean_bias_iteration_20mm_10mm.append(mean_bias_ROI_20mm_10mm)
    image_roughness_iteration_20mm_10mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_20mm_10mm])
   
    std_mean_ROI_40mm_2mm = [np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_2mm]
    mean_mean_ROI_40mm_2mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_2mm] #different methods mean fixed iteration
    mean_bias_ROI_40mm_2mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_40mm_2mm]
    std_mean_mean_ROI_40mm_2mm = [std_mean_ROI_40mm_2mm, mean_mean_ROI_40mm_2mm]
    std_mean_iteration_40mm_2mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_40mm_2mm)])  
    mean_bias_iteration_40mm_2mm.append(mean_bias_ROI_40mm_2mm)
    image_roughness_iteration_40mm_2mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_40mm_2mm])
    
    std_mean_ROI_40mm_3mm=[np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_3mm]
    mean_mean_ROI_40mm_3mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_3mm] #different methods mean fixed iteration
    mean_bias_ROI_40mm_3mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_40mm_3mm]
    std_mean_mean_ROI_40mm_3mm = [std_mean_ROI_40mm_3mm, mean_mean_ROI_40mm_3mm]
    std_mean_iteration_40mm_3mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_40mm_3mm)])  
    mean_bias_iteration_40mm_3mm.append(mean_bias_ROI_40mm_3mm)
    image_roughness_iteration_40mm_3mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_40mm_3mm])
   
    std_mean_ROI_40mm_5mm=[np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_5mm]
    mean_mean_ROI_40mm_5mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_5mm] #different methods mean fixed iteration
    mean_bias_ROI_40mm_5mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_40mm_5mm]
    std_mean_mean_ROI_40mm_5mm = [std_mean_ROI_40mm_5mm, mean_mean_ROI_40mm_5mm]
    std_mean_iteration_40mm_5mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_40mm_5mm)])  
    mean_bias_iteration_40mm_5mm.append(mean_bias_ROI_40mm_5mm)
    image_roughness_iteration_40mm_5mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_40mm_5mm])
    
    std_mean_ROI_40mm_10mm = [np.std(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_10mm]
    mean_mean_ROI_40mm_10mm = [np.mean(mean_ROI_elem) for mean_ROI_elem in mean_ROI_array_40mm_10mm] #different methods mean fixed iteration
    mean_bias_ROI_40mm_10mm = [np.mean(bias_ROI_elem) for bias_ROI_elem in bias_ROI_array_40mm_10mm]
    std_mean_mean_ROI_40mm_10mm = [std_mean_ROI_40mm_10mm, mean_mean_ROI_40mm_10mm]
    std_mean_iteration_40mm_10mm.append([std_mean_ROI_elem/mean_mean_ROI_elem * 100 for std_mean_ROI_elem, mean_mean_ROI_elem in zip(*std_mean_mean_ROI_40mm_10mm)])  
    mean_bias_iteration_40mm_10mm.append(mean_bias_ROI_40mm_10mm)
    image_roughness_iteration_40mm_10mm.append([np.mean(mean_std_elem) for mean_std_elem in mean_std_array_40mm_10mm])
    
    
   
    ####
    #std_mean_ROI_array = [std_ROI_array, mean_ROI_array]
    #avg_std_ROI = [[std_ROI_elem/mean_ROI_elem for std_ROI_elem, mean_ROI_elem in zip(*[std_ROI_array_elem, mean_ROI_array_elem])] for std_ROI_array_elem, mean_ROI_array_elem in zip(*std_mean_ROI_array)]
    #mean_avg_std_ROI = [np.mean(avg_std_ROI_elem) for avg_std_ROI_elem in avg_std_ROI]
    #avg_std_ROI_iteration.append(mean_avg_std_ROI)
    
    
    #im0 = ax0.imshow(img_avg[20,:,:], cmap='jet')
    #ax0.set_xlabel('x')
    #ax0.set_ylabel('y')
    #ax0.set_title('Avg trueLSOdensity')
    #fig.colorbar(im0,ax = ax0)
    #im0.set_clim(0,2)
    
        
    
    #plt.figure(figsize=(10,3.5))
    #im1 = ax1.imshow(img_avg[20,:,:], cmap='jet')
    #ax1.set_xlabel('x')
    #ax1.set_ylabel('y')
    #ax1.set_title('Avg 2timesLSOdensity')
    #fig.colorbar(im1,ax = ax1)
    #im1.set_clim(0,2)



contrast_recovery_iteration_zipzip = [[zip(*contrast_recovery_iteration_zip_elem_elem) for contrast_recovery_iteration_zip_elem_elem in contrast_recovery_iteration_zip_elem] for contrast_recovery_iteration_zip_elem in zip(*contrast_recovery_iteration)]
allreconmethods = []
allreconmethods_std = []
for contrast_zipzip_elem in contrast_recovery_iteration_zipzip:
   allesionsinalliterations = []
   allesionsinalliterations_std = []
   for contrast_zipzip_elemelem in contrast_zipzip_elem:
       allesionsinoneiteration = []
       allesionsinoneitreation_std = []
       for eachlesionwithallrealizations in contrast_zipzip_elemelem:
          allesionsinoneiteration.append(np.mean(eachlesionwithallrealizations)) #this saves mean contrast recovery from multiple realizations of 12 lesions only in our iteration
          allesionsinoneitreation_std.append(np.std(eachlesionwithallrealizations))
       allesionsinalliterations.append(allesionsinoneiteration) #this saves 50 iterations of each iteration
       allesionsinalliterations_std.append(allesionsinoneitreation_std)
   allreconmethods.append(allesionsinalliterations) #this saves different recon methods
   allreconmethods_std.append(allesionsinalliterations_std)
   
contrast_recovery_ordered = list([zip(*allesionsinoneiteration_cp) for allesionsinoneiteration_cp in allreconmethods])
contrast_recovery_ordered_std = list([zip(*allesionsinoneiteration_std_cp) for allesionsinoneiteration_std_cp in allreconmethods_std])
contrast_recovery_ordered_std_zip = list(zip(*contrast_recovery_ordered_std))
iteration = [10,50,90,130,170,210,250]

#iteration = range(10,210,10)
#average contrast recovery of multiple noise realizations

index = 0;
title = ['contrast recovery of 2mm at 20mm', 'contrast recovery of 3mm at 20mm', 'contrast recovery of 5mm at 20mm', 'contrast recovery of 10mm at 20mm', 
         'contrast recovery of 2mm at 40mm', 'contrast recovery of 3mm at 40mm', 'contrast recovery of 5mm at 40mm', 'contrast recovery of 10mm at 40mm'] 
legend_properties = {'weight':'normal'}
for contrast_recovery_ordered_elem in zip(*contrast_recovery_ordered):
    contrast_recovery_ordered_std_elem = contrast_recovery_ordered_std_zip[index]
   # plt.subplot(4,3,index+1)
    plt.figure(figsize=(5,3))
    plt.rc('font', weight='normal')
   # plt.plot(iteration,contrast_recovery_ordered_elem[0],'r.-', iteration,contrast_recovery_ordered_elem[1],'g.-',iteration,contrast_recovery_ordered_elem[2],'b.-', iteration, contrast_recovery_ordered_elem[3],'m.-', iteration, contrast_recovery_ordered_elem[4],'y.-')
    plt.errorbar(iteration,contrast_recovery_ordered_elem[0], yerr=contrast_recovery_ordered_std_elem[0], fmt='r.-', ecolor='r', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[1], yerr=contrast_recovery_ordered_std_elem[1], fmt='g.-', ecolor='g', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[2], yerr=contrast_recovery_ordered_std_elem[2], fmt='b.-', ecolor='b', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[3], yerr=contrast_recovery_ordered_std_elem[3], fmt='m.-', ecolor='m', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[4], yerr=contrast_recovery_ordered_std_elem[4], fmt='y.-', ecolor='y', capsize=1, elinewidth=0.5, markeredgewidth=1)
#    plt.errorbar(iteration,contrast_recovery_ordered_elem[5], yerr=contrast_recovery_ordered_std_elem[5], fmt='k.-', ecolor='k', capsize=1, elinewidth=0.5, markeredgewidth=1)    
    plt.title(title[index],fontweight = 'normal')
    plt.xlabel('iterations',fontweight = 'normal')
    plt.ylabel('contrast recovery',fontweight = 'normal')
    plt.axis([0, 250, 0, 1.5],fontweight = 'normal')
    plt.legend(('truePSF', 'noPSF'),loc = 'best', fontsize = 'xx-small',prop = legend_properties, ncol = 3)
    index = index + 1
  


"""
index = 0;
title = ['contrast recovery of 3mm at 15mm', 'contrast recovery of 4mm at 15mm', 'contrast recovery of 5mm at 15mm', 'contrast recovery of 6mm at 15mm', 
         'contrast recovery of 3mm at 30mm', 'contrast recovery of 4mm at 30mm', 'contrast recovery of 5mm at 30mm', 'contrast recovery of 6mm at 30mm', 
         'contrast recovery of 3mm at 45mm', 'contrast recovery of 4mm at 45mm', 'contrast recovery of 5mm at 45mm', 'contrast recovery of 6mm at 45mm']
legend_properties = {'weight':'normal'}
for contrast_recovery_ordered_elem in zip(*contrast_recovery_ordered):
    contrast_recovery_ordered_std_elem = contrast_recovery_ordered_std_zip[index]
   # plt.subplot(4,3,index+1)
    plt.figure(figsize=(5,3))
    plt.rc('font', weight='normal')
   # plt.plot(iteration,contrast_recovery_ordered_elem[0],'r.-', iteration,contrast_recovery_ordered_elem[1],'g.-',iteration,contrast_recovery_ordered_elem[2],'b.-', iteration, contrast_recovery_ordered_elem[3],'m.-', iteration, contrast_recovery_ordered_elem[4],'y.-')
    plt.errorbar(iteration,contrast_recovery_ordered_elem[0], yerr=contrast_recovery_ordered_std_elem[0], fmt='r.-', ecolor='r', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[1], yerr=contrast_recovery_ordered_std_elem[1], fmt='g.-', ecolor='g', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[2], yerr=contrast_recovery_ordered_std_elem[2], fmt='b.-', ecolor='b', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[3], yerr=contrast_recovery_ordered_std_elem[3], fmt='m.-', ecolor='m', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[4], yerr=contrast_recovery_ordered_std_elem[4], fmt='y.-', ecolor='y', capsize=1, elinewidth=0.5, markeredgewidth=1)
    plt.errorbar(iteration,contrast_recovery_ordered_elem[5], yerr=contrast_recovery_ordered_std_elem[5], fmt='k.-', ecolor='k', capsize=1, elinewidth=0.5, markeredgewidth=1)    
    plt.title(title[index],fontweight = 'normal')
    plt.xlabel('iterations',fontweight = 'normal')
    plt.ylabel('contrast recovery',fontweight = 'normal')
    plt.axis([0, 250, 0, 1.5],fontweight = 'normal')
    plt.legend(('truePSF', 'kernel1','kernel2', 'kernel3', 'kernel4', 'noPSF'),loc = 'best', fontsize = 'xx-small',prop = legend_properties, ncol = 3)
    index = index + 1
  
"""    

#noise

#plt.figure(figsize=(5,3.5))
std_mean_iteration_unpack_20mm_2mm = list(zip(*std_mean_iteration_20mm_2mm))
mean_bias_iteration_unpack_20mm_2mm = list(zip(*mean_bias_iteration_20mm_2mm))
image_roughness_iteration_unpack_20mm_2mm = list(zip(*image_roughness_iteration_20mm_2mm))
image_roughness_iteration_unpack_20mm_2mm  = np.multiply(image_roughness_iteration_unpack_20mm_2mm, 100)
std_mean_iteration_unpack_20mm_3mm = list(zip(*std_mean_iteration_20mm_3mm))
mean_bias_iteration_unpack_20mm_3mm = list(zip(*mean_bias_iteration_20mm_3mm))
image_roughness_iteration_unpack_20mm_3mm = list(zip(*image_roughness_iteration_20mm_3mm))
image_roughness_iteration_unpack_20mm_3mm = np.multiply(image_roughness_iteration_unpack_20mm_3mm, 100)
std_mean_iteration_unpack_20mm_5mm = list(zip(*std_mean_iteration_20mm_5mm))
mean_bias_iteration_unpack_20mm_5mm = list(zip(*mean_bias_iteration_20mm_5mm))
image_roughness_iteration_unpack_20mm_5mm = list(zip(*image_roughness_iteration_20mm_5mm))
image_roughness_iteration_unpack_20mm_5mm = np.multiply(image_roughness_iteration_unpack_20mm_5mm, 100)
std_mean_iteration_unpack_20mm_10mm = list(zip(*std_mean_iteration_20mm_10mm))
mean_bias_iteration_unpack_20mm_10mm = list(zip(*mean_bias_iteration_20mm_10mm))
image_roughness_iteration_unpack_20mm_10mm = list(zip(*image_roughness_iteration_20mm_10mm))
image_roughness_iteration_unpack_20mm_10mm = np.multiply(image_roughness_iteration_unpack_20mm_10mm, 100)

std_mean_iteration_unpack_40mm_2mm = list(zip(*std_mean_iteration_40mm_2mm))
mean_bias_iteration_unpack_40mm_2mm = list(zip(*mean_bias_iteration_40mm_2mm))
image_roughness_iteration_unpack_40mm_2mm = list(zip(*image_roughness_iteration_40mm_2mm))
image_roughness_iteration_unpack_40mm_2mm  = np.multiply(image_roughness_iteration_unpack_40mm_2mm, 100)
std_mean_iteration_unpack_40mm_3mm = list(zip(*std_mean_iteration_40mm_3mm))
mean_bias_iteration_unpack_40mm_3mm = list(zip(*mean_bias_iteration_40mm_3mm))
image_roughness_iteration_unpack_40mm_3mm = list(zip(*image_roughness_iteration_40mm_3mm))
image_roughness_iteration_unpack_40mm_3mm = np.multiply(image_roughness_iteration_unpack_40mm_3mm, 100)
std_mean_iteration_unpack_40mm_5mm = list(zip(*std_mean_iteration_40mm_5mm))
mean_bias_iteration_unpack_40mm_5mm = list(zip(*mean_bias_iteration_40mm_5mm))
image_roughness_iteration_unpack_40mm_5mm = list(zip(*image_roughness_iteration_40mm_5mm))
image_roughness_iteration_unpack_40mm_5mm = np.multiply(image_roughness_iteration_unpack_40mm_5mm, 100)
std_mean_iteration_unpack_40mm_10mm = list(zip(*std_mean_iteration_40mm_10mm))
mean_bias_iteration_unpack_40mm_10mm = list(zip(*mean_bias_iteration_40mm_10mm))
image_roughness_iteration_unpack_40mm_10mm = list(zip(*image_roughness_iteration_40mm_10mm))
image_roughness_iteration_unpack_40mm_10mm = np.multiply(image_roughness_iteration_unpack_40mm_10mm, 100)












"""
avg_std_ROI_iteration_unpack = zip(*avg_std_ROI_iteration)
plt.figure(figsize=(5,3))
plt.plot(iteration, std_mean_iteration_unpack[0],'r.-', iteration, std_mean_iteration_unpack[1],'g.-', iteration, std_mean_iteration_unpack[2], 'b.-', iteration, std_mean_iteration_unpack[3], 'm.-', iteration, std_mean_iteration_unpack[4], 'y.-', iteration, std_mean_iteration_unpack[5], 'k.-')
plt.title('of 20th image slice')
plt.xlabel('iteration')
plt.ylabel('emsemble noise')
"""

#I will comment out the following

index = 0;
title = ['2mm lesions', '3mm lesions', '5mm lesions', '10mm lesions'] 
#plt.axis([0, 1.2, 2, 10],fontweight = 'normal')
        # 'ensemble noise vs CRC of 3mm at 30mm', 'ensemble noise vs CRC of 4mm at 30mm', 'ensemble noise vs CRC of 5mm at 30mm', 'ensemble noise vs CRC of 6mm at 30mm', 
        # 'ensemble noise vs CRC of 3mm at 45mm', 'ensemble noise vs CRC of 4mm at 45mm', 'ensemble noise vs CRC of 5mm at 45mm', 'ensemble noise vs CRC of 6mm at 45mm']
"""
fig, axs = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True, figsize=(10,6))
#plt.xlabel('ensemble noise/%',fontweight = 'normal')
#plt.ylabel('CRC', fontweight = 'normal')
legend_properties = {'weight':'normal'}
for contrast_recovery_ordered_elem in zip(*contrast_recovery_ordered):
    contrast_recovery_ordered_std_elem = contrast_recovery_ordered_std_zip[index]
   # plt.subplot(4,3,index+1)  
    ax = axs[index/4, index%4]
    plt.rc('font', weight='normal')
   # plt.plot(iteration,contrast_recovery_ordered_elem[0],'r.-', iteration,contrast_recovery_ordered_elem[1],'g.-',iteration,contrast_recovery_ordered_elem[2],'b.-', iteration, contrast_recovery_ordered_elem[3],'m.-', iteration, contrast_recovery_ordered_elem[4],'y.-')
    ax.errorbar(list(ensemble_noise_matrix[:,0]),contrast_recovery_ordered_elem[0], yerr=contrast_recovery_ordered_std_elem[0], fmt='r.-', ecolor='r', capsize=1, elinewidth=0.5, markeredgewidth=1)
    ax.errorbar(list(ensemble_noise_matrix[:,1]),contrast_recovery_ordered_elem[1], yerr=contrast_recovery_ordered_std_elem[1], fmt='g.-', ecolor='g', capsize=1, elinewidth=0.5, markeredgewidth=1)
    ax.errorbar(list(ensemble_noise_matrix[:,2]),contrast_recovery_ordered_elem[2], yerr=contrast_recovery_ordered_std_elem[2], fmt='b.-', ecolor='b', capsize=1, elinewidth=0.5, markeredgewidth=1)
    ax.errorbar(list(ensemble_noise_matrix[:,3]),contrast_recovery_ordered_elem[3], yerr=contrast_recovery_ordered_std_elem[3], fmt='m.-', ecolor='m', capsize=1, elinewidth=0.5, markeredgewidth=1)
    ax.errorbar(list(ensemble_noise_matrix[:,4]),contrast_recovery_ordered_elem[4], yerr=contrast_recovery_ordered_std_elem[4], fmt='y.-', ecolor='y', capsize=1, elinewidth=0.5, markeredgewidth=1)
    ax.errorbar(list(ensemble_noise_matrix[:,5]),contrast_recovery_ordered_elem[5], yerr=contrast_recovery_ordered_std_elem[5], fmt='k.-', ecolor='k', capsize=1, elinewidth=0.5, markeredgewidth=1)    
    #ax.set_title(title[index],fontweight = 'normal')
    if index/4 == 0 and index%4 == 0:
       ax.set_title(title[0],fontweight = 'normal')
    if index/4 == 0 and index%4 == 1:
       ax.set_title(title[1],fontweight = 'normal')
    if index/4 == 0 and index%4 == 2:
       ax.set_title(title[2],fontweight = 'normal')
    if index/4 == 0 and index%4 == 3:
       ax.set_title(title[3],fontweight = 'normal')
    if index/4 == 2:
       ax.set_xlabel('$EN_{hs}$/%',fontweight = 'normal')
    if index%4 == 0 and index/4 == 0:
       ax.set_ylabel('CRC at 15mm offset',fontweight = 'normal')
    if index%4 == 0 and index/4 == 1:
       ax.set_ylabel('CRC at 30mm offset',fontweight = 'normal')
    if index%4 == 0 and index/4 == 2:
       ax.set_ylabel('CRC at 45mm offset',fontweight = 'normal')        
    #ax.axis([1, 10, 0, 1.2],fontweight = 'normal')
    #ax.legend(('truePSF', 'kernel1', 'kernel2', 'kernel3', 'kernel4','noPSF'),loc = 'best', fontsize = 'xx-small',prop = legend_properties, ncol = 2)
    index = index + 1
"""
fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, figsize=(10,5))
plt.rc('font', weight='bold')


axs[0,0].plot(std_mean_iteration_unpack_20mm_2mm[0], mean_bias_iteration_unpack_20mm_2mm[0], 'r.-',  std_mean_iteration_unpack_20mm_2mm[1], mean_bias_iteration_unpack_20mm_2mm[1],'g.-', std_mean_iteration_unpack_20mm_2mm[2], mean_bias_iteration_unpack_20mm_2mm[2], 'b.-', std_mean_iteration_unpack_20mm_2mm[3], mean_bias_iteration_unpack_20mm_2mm[3], 'm.-', std_mean_iteration_unpack_20mm_2mm[4],  mean_bias_iteration_unpack_20mm_2mm[4], 'y.-')
axs[0,0].set_title('2mm lesions', fontweight = 'bold')
#axs[0,0].set_xlabel('ensemble noise in hot ROIs/%')
axs[0,0].set_ylabel('bias at 20mm offset/%',fontsize = 'small', fontweight = 'bold')
legend_properties = {'weight':'bold'}
axs[0,0].legend(('truePSF', 'kernel1','kernel2', 'kernel3', 'noPSF'), bbox_to_anchor=(4.8, 0.3), loc=2, fontsize = 'x-small',prop = legend_properties, ncol = 1)
axs[0,0].grid()
#axs[0,0].axis([1, 8, -100, 25],fontweight = 'normal')
#axs[0,0].legend(loc=2,  bbox_to_anchor=(1,0.815),shadow=True, ncol=2)

#bbox_to_anchor=(0., 1.2, 1., .25) loc = 3

axs[0,1].plot(std_mean_iteration_unpack_20mm_3mm[0], mean_bias_iteration_unpack_20mm_3mm[0], 'r.-', std_mean_iteration_unpack_20mm_3mm[1], mean_bias_iteration_unpack_20mm_3mm[1], 'g.-', std_mean_iteration_unpack_20mm_3mm[2], mean_bias_iteration_unpack_20mm_3mm[2], 'b.-', std_mean_iteration_unpack_20mm_3mm[3], mean_bias_iteration_unpack_20mm_3mm[3], 'm.-', std_mean_iteration_unpack_20mm_3mm[4], mean_bias_iteration_unpack_20mm_3mm[4], 'y.-')
axs[0,1].set_title('3mm lesions', fontweight = 'bold')
axs[0,1].grid()

#axs[0,0]..xlabel('ensemble noise in hot ROIs/%')
#axs[0,1].set_ylabel('bias/%')
#axs[0,1]..axis([1, 8, -100, 25],fontweight = 'normal')
#axs[0,1].legend(('truePSF', 'kernel1','kernel2', 'kernel3', 'kernel4', 'noPSF'), loc = 'best', fontsize = 'xxxx-small',prop = legend_properties, ncol = 2)

axs[0,2].plot(std_mean_iteration_unpack_20mm_5mm[0], mean_bias_iteration_unpack_20mm_5mm[0], 'r.-', std_mean_iteration_unpack_20mm_5mm[1], mean_bias_iteration_unpack_20mm_5mm[1], 'g.-',  std_mean_iteration_unpack_20mm_5mm[2], mean_bias_iteration_unpack_20mm_5mm[2], 'b.-',  std_mean_iteration_unpack_20mm_5mm[3], mean_bias_iteration_unpack_20mm_5mm[3],'m.-', std_mean_iteration_unpack_20mm_5mm[4], mean_bias_iteration_unpack_20mm_5mm[4], 'y.-',)
axs[0,2].set_title('5mm lesions', fontweight = 'bold')
axs[0,2].grid()

#axs[0,2].set_ylabel('bias/%')
#plt.axis([1, 8, -100, 25],fontweight = 'normal')
#axs[0,2].legend(('truePSF', 'kernel1','kernel2', 'kernel3', 'kernel4', 'noPSF'), loc = 'best', fontsize = 'xxxx-small',prop = legend_properties, ncol = 2)

axs[0,3].plot(std_mean_iteration_unpack_20mm_10mm[0], mean_bias_iteration_unpack_20mm_10mm[0], 'r.-', std_mean_iteration_unpack_20mm_10mm[1], mean_bias_iteration_unpack_20mm_10mm[1], 'g.-',  std_mean_iteration_unpack_20mm_10mm[2], mean_bias_iteration_unpack_20mm_10mm[2], 'b.-',  std_mean_iteration_unpack_20mm_10mm[3], mean_bias_iteration_unpack_20mm_10mm[3],'m.-', std_mean_iteration_unpack_20mm_10mm[4], mean_bias_iteration_unpack_20mm_10mm[4], 'y.-')
axs[0,3].set_title('10mm lesions', fontweight = 'bold')
axs[0,3].grid()

#axs[0,2].set_ylabel('bias/%')
#plt.axis([1, 8, -100, 25],fontweight = 'normal')
#axs[0,2].legend(('truePSF', 'kernel1','kernel2', 'kernel3', 'kernel4', 'noPSF'), loc = 'best', fontsize = 'xxxx-small',prop = legend_properties, ncol = 2)

axs[1,0].plot(std_mean_iteration_unpack_40mm_2mm[0], mean_bias_iteration_unpack_40mm_2mm[0], 'r.-', std_mean_iteration_unpack_40mm_2mm[1], mean_bias_iteration_unpack_40mm_2mm[1], 'g.-',  std_mean_iteration_unpack_40mm_2mm[2], mean_bias_iteration_unpack_40mm_2mm[2], 'b.-',  std_mean_iteration_unpack_40mm_2mm[3], mean_bias_iteration_unpack_40mm_2mm[3],'m.-', std_mean_iteration_unpack_40mm_2mm[4], mean_bias_iteration_unpack_40mm_2mm[4], 'y.-')
#axs[0,2].set_title('5mm lesions')
axs[1,0].set_ylabel('bia at 40mm offset/%', fontsize = 'small', fontweight = 'bold')
axs[1,0].grid()

#plt.axis([1, 8, -100, 25],fontweight = 'normal')
#axs[0,2].legend(('truePSF', 'kernel1','kernel2', 'kernel3', 'kernel4', 'noPSF'), loc = 'best', fontsize = 'xxxx-small',prop = legend_properties, ncol = 2)
axs[1,0].set_xlabel('$EN_{hs}$/%',fontweight = 'bold')
axs[1,1].plot(std_mean_iteration_unpack_40mm_3mm[0], mean_bias_iteration_unpack_40mm_3mm[0], 'r.-', std_mean_iteration_unpack_40mm_3mm[1], mean_bias_iteration_unpack_40mm_3mm[1], 'g.-',  std_mean_iteration_unpack_40mm_3mm[2], mean_bias_iteration_unpack_40mm_3mm[2], 'b.-',  std_mean_iteration_unpack_40mm_3mm[3], mean_bias_iteration_unpack_40mm_3mm[3],'m.-', std_mean_iteration_unpack_40mm_3mm[4], mean_bias_iteration_unpack_40mm_3mm[4], 'y.-')
#axs[0,2].set_title('5mm lesions')
#axs[1,1].set_ylabel('bia at 30mm offset/%')
axs[1,1].set_xlabel('$EN_{hs}$/%',fontweight = 'bold')
axs[1,1].grid()


axs[1,2].plot(std_mean_iteration_unpack_40mm_5mm[0], mean_bias_iteration_unpack_40mm_5mm[0], 'r.-', std_mean_iteration_unpack_40mm_5mm[1], mean_bias_iteration_unpack_40mm_5mm[1], 'g.-',  std_mean_iteration_unpack_40mm_5mm[2], mean_bias_iteration_unpack_40mm_5mm[2], 'b.-',  std_mean_iteration_unpack_40mm_5mm[3], mean_bias_iteration_unpack_40mm_5mm[3],'m.-', std_mean_iteration_unpack_40mm_5mm[4], mean_bias_iteration_unpack_40mm_5mm[4], 'y.-')

axs[1,2].set_xlabel('$EN_{hs}$/%',fontweight = 'bold')
axs[1,2].grid()


axs[1,3].plot(std_mean_iteration_unpack_40mm_10mm[0], mean_bias_iteration_unpack_40mm_10mm[0], 'r.-', std_mean_iteration_unpack_40mm_10mm[1], mean_bias_iteration_unpack_40mm_10mm[1], 'g.-',  std_mean_iteration_unpack_40mm_10mm[2], mean_bias_iteration_unpack_40mm_10mm[2], 'b.-',  std_mean_iteration_unpack_40mm_10mm[3], mean_bias_iteration_unpack_40mm_10mm[3],'m.-', std_mean_iteration_unpack_40mm_10mm[4], mean_bias_iteration_unpack_40mm_10mm[4], 'y.-')

axs[1,3].set_xlabel('$EN_{hs}$/%',fontweight = 'bold')
axs[1,3].grid()









fig, ax1 = plt.subplots(ncols=5, sharex='col', sharey='row', figsize=(20,12))
fig.subplots_adjust(hspace=-0.20, wspace=0.0)
index = 0
for img_avg_elem in img_avg:
    if index==0:
        #ax1[int(index/3),int(index%3)].set_title('truePSF')
        ax1[int(index%5)].text(50,10,'truePSF',color='white', fontsize=15)
        ax1[int(index%5)].axis('off')
        f = open('H:\\alex_backup\\BACKUP_important\\VersaPET_article\\truePSF_250th_new2.v','wb')
        f.write(img_avg_elem)
        f.close()
    elif index == 1:
        ax1[int(index%5)].text(50,10,'kernel1',color='white', fontsize=15)
        ax1[int(index%5)].axis('off')
        f = open('H:\\alex_backup\\BACKUP_important\\VersaPET_article\\kernel1_250th_new2.v','wb')
        f.write(img_avg_elem)
        f.close()
    elif index == 2:
        ax1[int(index%5)].text(50,10,'kernel2',color='white', fontsize=15)
        ax1[int(index%5)].axis('off')
        f = open('H:\\alex_backup\\BACKUP_important\\VersaPET_article\\kernel3_250th_new2.v','wb')
        f.write(img_avg_elem)
        f.close()
    elif index == 3:
        ax1[int(index%5)].text(50,10,'kernel3',color='white', fontsize=15)
        ax1[int(index%5)].axis('off')        
        f = open('H:\\alex_backup\\BACKUP_important\\VersaPET_article\\kernel4_250th_new2.v','wb')
        f.write(img_avg_elem)
        f.close()
    else:
        ax1[int(index%5)].text(50,10,'noPSF',color='white', fontsize=15)
        ax1[int(index%5)].axis('off')

       # mask = np.zeros((127,127),dtype=float)
        #mask = circularMask(mask,[61.2,76.2],1.5)
        #mask = circularMask(mask,[76.2,91.2],2)
        #mask = circularMask(mask,[91.2,76.2],2.5)
        #mask = circularMask(mask,[76.2,61.2],3)
        #mask = circularMask(mask,[46.2,76.2],1.5)
        #mask = circularMask(mask,[76.2,106.2],2)
        #mask = circularMask(mask,[106.2,76.2],2.5)
        #mask = circularMask(mask,[76.2,46.2],3)
        #mask = circularMask(mask,[31.2,76.2],1.5)
        #mask = circularMask(mask,[76.2,121.2],2)
        #mask = circularMask(mask,[121.2,76.2],2.5)
       # print "separate"
       # mask = circularMask(mask,[76.2,31.2],3)
        
        
           # ROI_15mm_4mm = circularROI(img,[76.2,91.2],range(19,20),2)
           # ROI_15mm_5mm = circularROI(img,[91.2,76.2],range(19,20),2.5)
           # ROI_15mm_6mm = circularROI(img,[76.2,61.2],range(19,20),3)
           # ROI_30mm_3mm = circularROI(img,[46.2,76.2],range(19,20),1.5)
           # ROI_30mm_4mm = circularROI(img,[76.2,106.2],range(19,20),2)
           # ROI_30mm_5mm = circularROI(img,[106.2,76.2],range(19,20),2.5)
           # ROI_30mm_6mm = circularROI(img,[76.2,46.2],range(19,20),3)
           # ROI_45mm_3mm = circularROI(img,[31.2,76.2],range(19,20),1.5)
           # ROI_45mm_4mm = circularROI(img,[76.2,121.2],range(19,20),2)
           # ROI_45mm_5mm = circularROI(img,[121.2,76.2],range(19,20),2.5)
           # ROI_45mm_6mm = circularROI(img,[76.2,31.2],range(19,20),3)
            
        mask = np.zeros((153,153),dtype=float)
        #mask = circularMask(mask,[96.5,76.5],2.5) #2mm
        #mask = circularMask(mask,[56.5,76.5],2.5) #3mm
        #mask = circularMask(mask,[76.5,36.5],5) #10mm
        #mask = circularMask(mask,[76.5,96.5],2.5) #5mm
        mask = circularMask(mask,[76.5,76.5],60)
        img_avg_elem = np.multiply(img_avg_elem,mask)
        f = open('H:\\alex_backup\\BACKUP_important\\VersaPET_article\\noPSF_250th_new2.v','wb')
        f.write(img_avg_elem)
        f.close()
        
    im1 = ax1[int(index%5)].imshow(img_avg_elem, cmap='jet')
    index = index+1
    im1.set_clim(0,27)

    #ax1[index/3,index%3].set_xlabel('x')
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.4)
cax = plt.axes([0.81, 0.135, 0.02, 0.22]) 
fig.colorbar(im1,ax = ax1, cax=cax)
    
    
  
