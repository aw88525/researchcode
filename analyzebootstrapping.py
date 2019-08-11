# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 11:49:27 2018

@author: Alex Wei @vaskalab
this code is used to analyze the background noise in different bootstrapped images and original images
and report the average background uptake value, noise and lesion uptake value, and then perform statistical
tests to showcase if there is any significant differences between the bootstrapped results and initial results

"""



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from collections import defaultdict
import pandas as pd

# first I need to import the data
# original

string_original = []
string_original.append('F:\\bootstrappingresults\\original\\VersaPET_8spherephantom_1\\8spheres_50.v')
string_original.append('F:\\bootstrappingresults\\original\\VersaPET_8spherephantom_2\\8spheres_50.v')
string_original.append('F:\\bootstrappingresults\\original\\VersaPET_8spherephantom_3\\8spheres_50.v')
string_original.append('F:\\bootstrappingresults\\original\\VersaPET_8spherephantom_4\\8spheres_50.v')
string_original.append('F:\\bootstrappingresults\\original\\VersaPET_8spherephantom_5\\8spheres_50.v')
string_original.append('F:\\bootstrappingresults\\original\\VersaPET_8spherephantom_6\\8spheres_50.v')

# bootstrapped from the original length

string_avg1lenboot = []
string_avg1lenboot.append('F:\\bootstrappingresults\\1timeslengthbootstrapping\\VersaPET_8spherephantom_avg1len_boot1\\8spheres_avg1len_50.v')
string_avg1lenboot.append('F:\\bootstrappingresults\\1timeslengthbootstrapping\\VersaPET_8spherephantom_avg1len_boot2\\8spheres_avg1len_50.v')
string_avg1lenboot.append('F:\\bootstrappingresults\\1timeslengthbootstrapping\\VersaPET_8spherephantom_avg1len_boot3\\8spheres_avg1len_50.v')
string_avg1lenboot.append('F:\\bootstrappingresults\\1timeslengthbootstrapping\\VersaPET_8spherephantom_avg1len_boot4\\8spheres_avg1len_50.v')
string_avg1lenboot.append('F:\\bootstrappingresults\\1timeslengthbootstrapping\\VersaPET_8spherephantom_avg1len_boot5\\8spheres_avg1len_50.v')
string_avg1lenboot.append('F:\\bootstrappingresults\\1timeslengthbootstrapping\\VersaPET_8spherephantom_avg1len_boot6\\8spheres_avg1len_50.v')

# bootstrapped from 2 times length

string_avg2lenboot = []
string_avg2lenboot.append('F:\\bootstrappingresults\\2timeslengthbootstrapping\\VersaPET_8spherephantom_avg2len_boot1\\8spheres_avg2len_50.v')
string_avg2lenboot.append('F:\\bootstrappingresults\\2timeslengthbootstrapping\\VersaPET_8spherephantom_avg2len_boot2\\8spheres_avg2len_50.v')
string_avg2lenboot.append('F:\\bootstrappingresults\\2timeslengthbootstrapping\\VersaPET_8spherephantom_avg2len_boot3\\8spheres_avg2len_50.v')
string_avg2lenboot.append('F:\\bootstrappingresults\\2timeslengthbootstrapping\\VersaPET_8spherephantom_avg2len_boot4\\8spheres_avg2len_50.v')
string_avg2lenboot.append('F:\\bootstrappingresults\\2timeslengthbootstrapping\\VersaPET_8spherephantom_avg2len_boot5\\8spheres_avg2len_50.v')
string_avg2lenboot.append('F:\\bootstrappingresults\\2timeslengthbootstrapping\\VersaPET_8spherephantom_avg2len_boot6\\8spheres_avg2len_50.v')

# bootstrapped from 5 times length

string_avg5lenboot = []
string_avg5lenboot.append('F:\\bootstrappingresults\\5timeslengthbootstrapping\\VersaPET_8spherephantom_avg5len_boot1\\8spheres_avg5len_50.v')
string_avg5lenboot.append('F:\\bootstrappingresults\\5timeslengthbootstrapping\\VersaPET_8spherephantom_avg5len_boot2\\8spheres_avg5len_50.v')
string_avg5lenboot.append('F:\\bootstrappingresults\\5timeslengthbootstrapping\\VersaPET_8spherephantom_avg5len_boot3\\8spheres_avg5len_50.v')
string_avg5lenboot.append('F:\\bootstrappingresults\\5timeslengthbootstrapping\\VersaPET_8spherephantom_avg5len_boot4\\8spheres_avg5len_50.v')
string_avg5lenboot.append('F:\\bootstrappingresults\\5timeslengthbootstrapping\\VersaPET_8spherephantom_avg5len_boot5\\8spheres_avg5len_50.v')
string_avg5lenboot.append('F:\\bootstrappingresults\\5timeslengthbootstrapping\\VersaPET_8spherephantom_avg5len_boot6\\8spheres_avg5len_50.v')

# bootstrapped from 10 times length

string_avg10lenboot = []
string_avg10lenboot.append('F:\\bootstrappingresults\\10timeslengthbootstrapping\\VersaPET_8spherephantom_avg10len_boot1\\8spheres_avg10len_50.v')
string_avg10lenboot.append('F:\\bootstrappingresults\\10timeslengthbootstrapping\\VersaPET_8spherephantom_avg10len_boot2\\8spheres_avg10len_50.v')
string_avg10lenboot.append('F:\\bootstrappingresults\\10timeslengthbootstrapping\\VersaPET_8spherephantom_avg10len_boot3\\8spheres_avg10len_50.v')
string_avg10lenboot.append('F:\\bootstrappingresults\\10timeslengthbootstrapping\\VersaPET_8spherephantom_avg10len_boot4\\8spheres_avg10len_50.v')
string_avg10lenboot.append('F:\\bootstrappingresults\\10timeslengthbootstrapping\\VersaPET_8spherephantom_avg10len_boot5\\8spheres_avg10len_50.v')
string_avg10lenboot.append('F:\\bootstrappingresults\\10timeslengthbootstrapping\\VersaPET_8spherephantom_avg10len_boot6\\8spheres_avg10len_50.v')


# defines a function that draws 2D ROI from image
def circularROI(img, center, slice, radius):
    return [img[slice, x, y] for x in range(0,127) for y in range(0,127) if math.sqrt(((x+0.5)*1.2-center[0])**2+((y+0.5)*1.2-center[1])**2) <= radius]    

def sphericalROI(img, center, slice_range, radius):
    return [img[z, x, y]  for x in range(0,127) for y in range(0,127) for z in slice_range if math.sqrt(((x+0.5)*1.2-center[0])**2+((y+0.5)*1.2-center[1])**2 + ((z+0.5)*1.2-center[2])**2)  <= radius] 

def imganalysis(file_holder):
    bkg_vector = []
    noise_vector = []
    _3mm_lesion_vector = []
    _5mm_lesion_vector = []
    _10mm_lesion_vector = []
    for filename in file_holder:
        f = open(filename,'rb')
        data = np.fromfile(f,dtype=np.float32)
        f.close()
        img = np.reshape(data,[39,127,127])   #used to be 19.58
        bkgROI = []
        for slice in [10,12,14,26,28,30]:   #48 background ROIs
            bkgROI.append(circularROI(img,[76.2,76.2],slice,20))
        
        ROI_20mm_3mm = sphericalROI(img,[56.2,76.2,23.4],range(18,30),1.5)
        ROI_20mm_5mm = sphericalROI(img,[76.2,96.2,23.4],range(18,30),2.5)
        ROI_20mm_10mm = sphericalROI(img,[76.2,56.2,23.4],range(18,30),5)
         
        _3mm_lesion_vector.append(np.mean(ROI_20mm_3mm))
        _5mm_lesion_vector.append(np.mean(ROI_20mm_5mm))     
        _10mm_lesion_vector.append(np.mean(ROI_20mm_10mm))  
        bkg_vector.append(np.mean(bkgROI))
        noise_vector.append(np.std(bkgROI)/np.mean(bkgROI))
        
        
        
    return bkg_vector, noise_vector, _3mm_lesion_vector, _5mm_lesion_vector, _10mm_lesion_vector
        
[bkg_original, noise_original, _3mm_lesion_original, _5mm_lesion_original, _10mm_lesion_original] = imganalysis(string_original)
[bkg_av1lenboot, noise_av1lenboot, _3mm_lesion_av1lenboot, _5mm_lesion_av1lenboot, _10mm_lesion_av1lenboot] = imganalysis(string_avg1lenboot)
[bkg_av2lenboot, noise_av2lenboot, _3mm_lesion_av2lenboot, _5mm_lesion_av2lenboot, _10mm_lesion_av2lenboot] = imganalysis(string_avg2lenboot)    
[bkg_av5lenboot, noise_av5lenboot, _3mm_lesion_av5lenboot, _5mm_lesion_av5lenboot, _10mm_lesion_av5lenboot] = imganalysis(string_avg5lenboot) 
[bkg_av10lenboot, noise_av10lenboot, _3mm_lesion_av10lenboot, _5mm_lesion_av10lenboot, _10mm_lesion_av10lenboot] = imganalysis(string_avg10lenboot)


bkg_data = pd.DataFrame({'original':bkg_original, 'av1lenboot':bkg_av1lenboot, 'av2lenboot':bkg_av2lenboot, 'av5lenboot':bkg_av5lenboot, 'av10lenboot':bkg_av10lenboot})

plt.subplots(5,1,figsize=(12,12))
plt.subplot(511)
pd.DataFrame(bkg_data).boxplot()
plt.title('bkg average value')
noise_data = pd.DataFrame({'original':noise_original, 'av1lenboot':noise_av1lenboot, 'av2lenboot':noise_av2lenboot, 'av5lenboot':noise_av5lenboot, 'av10lenboot':noise_av10lenboot})
plt.subplot(512)
pd.DataFrame(noise_data).boxplot()
plt.title('bkg noise value')

_3mm_lesion_data = pd.DataFrame({'original': _3mm_lesion_original, 'av1lenboot':_3mm_lesion_av1lenboot, 'av2lenboot':_3mm_lesion_av2lenboot, 'av5lenboot':_3mm_lesion_av5lenboot,  \
                            'av10lenboot':_3mm_lesion_av10lenboot})
plt.subplot(513)
pd.DataFrame(_3mm_lesion_data).boxplot()
plt.title('3mm lesion average value')

_5mm_lesion_data = pd.DataFrame({'original': _5mm_lesion_original, 'av1lenboot':_5mm_lesion_av1lenboot, 'av2lenboot':_5mm_lesion_av2lenboot, 'av5lenboot':_5mm_lesion_av5lenboot,  \
                            'av10lenboot':_5mm_lesion_av10lenboot})

plt.subplot(514)
pd.DataFrame(_5mm_lesion_data).boxplot()
plt.title('5mm lesion average value')

_10mm_lesion_data = pd.DataFrame({'original': _10mm_lesion_original, 'av1lenboot':_10mm_lesion_av1lenboot, 'av2lenboot':_10mm_lesion_av2lenboot, 'av5lenboot':_10mm_lesion_av5lenboot,  \
                            'av10lenboot':_10mm_lesion_av10lenboot})

plt.subplot(515)
pd.DataFrame(_10mm_lesion_data).boxplot()
plt.title('10mm lesion average value')


plt.subplots_adjust(left  = 0.125,  # the left side of the subplots of the figure
right = 0.9,    # the right side of the subplots of the figure
bottom = 0.1 ,  # the bottom of the subplots of the figure
top = 0.9 ,     # the top of the subplots of the figure
wspace = 0.2 ,  # the amount of width reserved for blank space between subplots,
               # expressed as a fraction of the average axis width
hspace = 0.5   # the amount of height reserved for white space between subplots,
)              # expressed as a fraction of the average axis height)