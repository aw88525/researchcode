# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 13:11:11 2019

@author: wsy88
"""
from sklearn.utils import resample
import numpy as np

#data = range(1,3072*3071/2)
"""
files = []
for i in range(10):
    f = open('C:\\Users\\wsy88\\Documents\\MATLAB\\Hist_8spheres_newmerged_'+str(i+1)+'_mlem.s','rb')
    files.append(f)
    f.close()
"""
listdata = []
#length = []
size = int(3072*3071/2)
number = 1
for k in range(number):    
    f = open('C:\\Users\\wsy88\\Documents\\MATLAB\\Hist_8spheres_newmerged_'+str(k+1)+'_mlem.s','rb')
    print("this is " + str(k+1) + " files!")
    histogram = np.fromfile(f, dtype = np.float32)
    f.close()
    for i in range(0, size):
        try:
            if(int(histogram[i]) != 0):
                for j in range(0, int(histogram[i])):
                    listdata.append(i)
            else:
                continue
        except:
            print("something is wrong!")
    #length.append(len(listdata))  
 

avg_list_length = int(len(listdata) / number)
for ind in range(3,11):
    boot = resample(listdata, replace = True, n_samples=avg_list_length, random_state = ind)
    newhistogram = np.zeros(size)
    for i in range(0, avg_list_length):
        try:
            newhistogram[boot[i]] = newhistogram[boot[i]] + 1
        except:
            print("something is wrong again!")
            
    f = open('C:\\Users\\wsy88\\Documents\\Visual Studio 2015\\Projects\\MLEM_parallel\\Hist_8spheres_newmerged_mlem_boot'+str(ind)+'.s','wb')
    np.array(newhistogram, dtype = np.float32).tofile(f)
    f.close()


