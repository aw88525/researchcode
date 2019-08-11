# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 12:22:40 2019

@author: wsy88
"""
"""
# The goal of this code is to implement the DOI blurring especially along the
# crystal centerline direction with 2 mm filter size note that for the points 
# close to the crystal end, the blurring may end up outside of the crystal volume 
# and then in such case we will place the points to the surfaces of the crystal
"""


import numpy as np
import pandas as pd



a = [line.rstrip('\n').split('\t') for line in open('C:\\Users\\wsy88\\Documents\\MATLAB\\crystal_position_halfring_surface_20cm_flatpanel.txt').readlines()]
b = [line.rstrip('\n').split('\t') for line in open('C:\\Users\\wsy88\\Documents\\MATLAB\\crystal_position_halfring_center_20cm_flatpanel.txt').readlines()]

surface = np.array(a, dtype = 'float')
center = np.array(b, dtype = 'float')

n_crystal = 8320 # this is the number of crystals for flat panel /dshape
id_array_1 = []
id_array_2 = []

for i in range(n_crystal):
    for j in range(i+1, n_crystal):
        id = 0
        id = i * (n_crystal - 1) - (i - 1) * i / 2
        id += j - i - 1
        id_array_1.append(i)
        id_array_2.append(j)

assert len(id_array_1) == n_crystal * (n_crystal - 1) / 2


id_file = 'C:\\Users\\wsy88\\Documents\\Visual Studio 2015\\Projects\\MLEM_Raytracing_wDOI\\MLEM_Raytracing_wDOI\\list_id_partial_paul_only3mm_2400s_wdoi.s'
xyz1_file = 'C:\\Users\\wsy88\\Documents\\Visual Studio 2015\\Projects\\MLEM_Raytracing_wDOI\\MLEM_Raytracing_wDOI\\list_xyz1_partial_paul_3mmonly_2400s.s'
xyz2_file = 'C:\\Users\\wsy88\\Documents\\Visual Studio 2015\\Projects\\MLEM_Raytracing_wDOI\\MLEM_Raytracing_wDOI\\list_xyz2_partial_paul_3mmonly_2400s.s'

f = open(id_file,'rb')
pairid = np.fromfile(f,dtype = np.int32)
f.close()
f = open(xyz1_file, 'rb')
xyz1 = np.fromfile(f, dtype = np.float32)
f.close()
f = open(xyz2_file,'rb')
xyz2 = np.fromfile(f, dtype = np.float32)
f.close()

xyz1_reshape = np.reshape(xyz1, [int(len(xyz1)/3), 3])
xyz2_reshape = np.reshape(xyz2, [int(len(xyz2)/3), 3])

xyz1_fixedpoint_surface = [surface[id_array_1[pairid_elem], :] for pairid_elem in pairid]
xyz1_fixedpoint_center = [center[id_array_1[pairid_elem], :] for pairid_elem in pairid]

xyz2_fixedpoint_surface = [surface[id_array_2[pairid_elem], :] for pairid_elem in pairid]
xyz2_fixedpoint_center = [center[id_array_2[pairid_elem], :] for pairid_elem in pairid]

xyz1_surface = np.array(xyz1_fixedpoint_surface)
xyz1_center = np.array(xyz1_fixedpoint_center)	
xyz2_surface = np.array(xyz2_fixedpoint_surface)
xyz2_center = np.array(xyz2_fixedpoint_center)

xyz1_diff = xyz1_center - xyz1_surface
xyz2_diff = xyz2_center - xyz2_surface

xyz1_diff2 = xyz1_reshape - xyz1_surface
xyz2_diff2 = xyz2_reshape - xyz2_surface

def compute_projection(a1, a2):
    return np.array(np.sum((a1*a2), axis = 1, keepdims = True) / (np.linalg.norm(a2, axis = 1, keepdims = True))) 


scalar_proj1 = compute_projection(xyz1_diff, xyz1_diff2)

#np.array(np.sum((xyz1_diff * xyz1_diff2), axis = 1, keepdims = True) / (np.linalg.norm(xyz1_diff2, axis = 1, keepdims = True)))

# now add blurring
scalar_proj1_blur = np.random.normal(scalar_proj1, 0.84)
scalar_proj1_blur[scalar_proj1_blur < 0] = 0
scalar_proj1_blur[scalar_proj1_blur > 14] = 14



hitpoints1 = scalar_proj1_blur * (xyz1_diff / np.linalg.norm(xyz1_diff, axis = 1, keepdims = True)) + xyz1_surface


scalar_proj2 = compute_projection(xyz2_diff, xyz2_diff2)
#np.array(np.sum((xyz2_diff * xyz2_diff2), axis = 1, keepdims = True) / (np.linalg.norm(xyz2_diff2, axis = 1, keepdims = True)))

#### Now add blurring
scalar_proj2_blur = np.random.normal(scalar_proj2, 0.84)
scalar_proj2_blur[scalar_proj2_blur < 0] = 0
scalar_proj2_blur[scalar_proj2_blur > 14] = 14


hitpoints2 = scalar_proj2_blur * (xyz2_diff / np.linalg.norm(xyz2_diff, axis = 1, keepdims = True)) + xyz2_surface


f1 = open('C:\\Users\\wsy88\\Documents\\Visual Studio 2015\\Projects\\MLEM_Raytracing_wDOI\\list_xyz1_doiblurring_dshape_3mmonly_2400s.s','wb')
f2 = open('C:\\Users\\wsy88\\Documents\\Visual Studio 2015\\Projects\\MLEM_Raytracing_wDOI\\list_xyz2_doiblurring_dshape_3mmonly_2400s.s', 'wb')
hitpoints1.astype(np.float32).tofile(f1)
hitpoints2.astype(np.float32).tofile(f2)
f1.close()
f2.close()

'''

for (int i = 0; i < n_crystal; i++) {
	for (int j = i + 1; j < n_crystal; j++) {
		int id = 0;
		id = i*(n_crystal - 1) - (i - 1)*i / 2;
		id += j - i - 1;
		index_array_1[id] = i;
		index_array_2[id] = j;
	}
}
'''