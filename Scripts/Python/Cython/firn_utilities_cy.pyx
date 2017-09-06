# Description: Cython implementation of firn_utilities
# Optional compiler directives (must be at top): #cython: boundscheck=False, wraparound=False, cdivision=True
#
# Author: Christian Steger, April 2016

# Load modules
from libc.math cimport isnan, NAN
import numpy as np

######################################################################################################################## 

cpdef firnproresamp(float[:,::1] data_in, float[:,::1] bound_in, float[::1] bound_out, bytes cons_quant):
    """Description: Conservative resampling of firn profiles to common grid
       Input:
       data_in:    numerical array with data (layer, time) (may contain sequential nan-elements at end)
       bound_in:   numerical array with input boundaries (layer + 1 (monotonously increasing), time)
       bound_out:  numerical array with output boundaries (monotonously increasing)
       cons_quant: character string with quantitiy that should be conserved. Either sum or weighted_mean
       Output:
       data_out:   numerical array with output data (layer, time) (may contain sequential nan-elements)
       Author: Christian Steger, March 2016"""
       
    cdef int N = data_in.shape[0]
    cdef int O = data_in.shape[1]
    cdef int P = bound_out.shape[0]
    cdef int k, m
    cdef int ind_max, ind_out, ind_in
    cdef float dist, dist_in, dist_out
    cdef bint flag_ini
    cdef float[:, :] data_out = np.empty((P - 1, O), dtype = np.float32)

    # Check input arguments
    if (cons_quant != "sum" and cons_quant != "weighted_mean"):
        print "Error: cons_quant is not a character string with either sum or weighted_mean"
        return
        
    data_out[:,:] = NAN
    
    # Resample data
    if (cons_quant == "weighted_mean"): # weighted mean is conserved
         for k in range(O): # loop over time steps
            ind_max = 0
            for m in range(N):
                if (not isnan(data_in[m,k])):
                    ind_max += 1
                else:
                    break
            if ((ind_max > 0) and (bound_out[P - 2] >= bound_in[0,k])): # short-circuit evaluation
                for m in range(P):
                    if (bound_out[m] >= bound_in[0,k]):
                        ind_out = m
                        break
                for m in range(ind_max + 1):
                    if (bound_in[m,k] <= bound_out[ind_out]):
                        ind_in = m
                    else:
                        break
                dist = bound_out[ind_out]
                flag_ini = True
                while ((dist < bound_out[P - 1]) and (bound_in[ind_max,k] >= bound_out[ind_out + 1])): # short-circuit evaluation
                    if (flag_ini):
                        data_out[ind_out,k] = 0.
                    dist_in = bound_in[(ind_in + 1),k] - dist
                    dist_out = bound_out[ind_out + 1] - dist
                    if (dist_in < dist_out):
                        data_out[ind_out,k] += data_in[ind_in,k] * dist_in / \
                                               (bound_out[ind_out + 1] - bound_out[ind_out])
                        ind_in += 1
                        dist = bound_in[ind_in,k]
                        flag_ini = False
                    elif (dist_in > dist_out):
                        data_out[ind_out,k] += data_in[ind_in,k] * dist_out / \
                                               (bound_out[ind_out + 1] - bound_out[ind_out])
                        ind_out += 1
                        dist = bound_out[ind_out]
                        flag_ini = True
                    else:
                        data_out[ind_out,k] += data_in[ind_in,k] * dist_in / \
                                               (bound_out[ind_out + 1] - bound_out[ind_out])
                        ind_in += 1
                        ind_out += 1
                        dist = bound_in[ind_in,k]
                        flag_ini = True
    else: # sum is conserved
         for k in range(O): # loop over time steps
            ind_max = 0
            for m in range(N):
                if (not isnan(data_in[m,k])):
                    ind_max += 1
                else:
                    break
            if ((ind_max > 0) and (bound_out[P - 2] >= bound_in[0,k])): # short-circuit evaluation
                for m in range(P):
                    if (bound_out[m] >= bound_in[0,k]):
                        ind_out = m
                        break
                for m in range(ind_max + 1):
                    if (bound_in[m,k] <= bound_out[ind_out]):
                        ind_in = m
                    else:
                        break
                dist = bound_out[ind_out]
                flag_ini = True
                while ((dist < bound_out[P - 1]) and (bound_in[ind_max,k] >= bound_out[ind_out + 1])): # short-circuit evaluation
                    if (flag_ini):
                        data_out[ind_out,k] = 0.
                    dist_in = bound_in[(ind_in + 1),k] - dist
                    dist_out = bound_out[ind_out + 1] - dist
                    if (dist_in < dist_out):
                        data_out[ind_out,k] += data_in[ind_in,k] * dist_in / \
                                               (bound_in[(ind_in + 1),k] - bound_in[ind_in,k])
                        ind_in += 1
                        dist = bound_in[ind_in,k]
                        flag_ini = False
                    elif (dist_in > dist_out):
                        data_out[ind_out,k] += data_in[ind_in,k] * dist_out / \
                                               (bound_in[(ind_in + 1),k] - bound_in[ind_in,k])
                        ind_out += 1
                        dist = bound_out[ind_out]
                        flag_ini = True
                    else:
                        data_out[ind_out,k] += data_in[ind_in,k] * dist_in / \
                                               (bound_in[(ind_in + 1),k] - bound_in[ind_in,k])
                        ind_in += 1
                        ind_out += 1
                        dist = bound_in[ind_in,k]
                        flag_ini = True                       
    return data_out

######################################################################################################################## 

cpdef firnproicelay(float[:,::1] data_in, float[:,::1] bound_in, float[::1] bound_out, float thresh_dens):
    """Description: Assigns ice layer to common grid
       Input:
       data_in:     numerical array with ice density (layer, time) (may contain sequential nan-elements at end) [kg m-3]
       bound_in:    numerical array with input boundaries (layer + 1 (monotonously increasing), time)
       bound_out:   numerical array with output boundaries (monotonously increasing)
       thresh_dens: numerical value with density threshold for ice layer [kg m-3]
       Output:
       data_out:   numerical array with ice layer fraction (layer, time) (may contain sequential nan-elements) [-]
       Author: Christian Steger, May 2016"""
       
    cdef int N = data_in.shape[0]
    cdef int O = data_in.shape[1]
    cdef int P = bound_out.shape[0]
    cdef int k, m
    cdef int ind_max, ind_out, ind_in
    cdef float dist, dist_in, dist_out
    cdef bint flag_ini
    cdef float[:, :] data_out = np.empty((P - 1, O), dtype = np.float32)
        
    data_out[:,:] = NAN
    
    # Assign ice layers
    for k in range(O): # loop over time steps
        ind_max = 0
        for m in range(N):
            if (not isnan(data_in[m,k])):
                ind_max += 1
            else:
                break
        if ((ind_max > 0) and (bound_out[P - 2] >= bound_in[0,k])): # short-circuit evaluation
            for m in range(P):
                if (bound_out[m] >= bound_in[0,k]):
                    ind_out = m
                    break
            for m in range(ind_max + 1):
                if (bound_in[m,k] <= bound_out[ind_out]):
                    ind_in = m
                else:
                    break
            dist = bound_out[ind_out]
            flag_ini = True
            while ((dist < bound_out[P - 1]) and (bound_in[ind_max,k] >= bound_out[ind_out + 1])): # short-circuit evaluation
                if (flag_ini):
                    data_out[ind_out,k] = 0.
                dist_in = bound_in[(ind_in + 1),k] - dist
                dist_out = bound_out[ind_out + 1] - dist
                if (dist_in < dist_out):
                    if (data_in[ind_in,k] >= thresh_dens): # input is ice layer
                        data_out[ind_out,k] += dist_in / (bound_out[ind_out + 1] - bound_out[ind_out])
                    ind_in += 1
                    dist = bound_in[ind_in,k]
                    flag_ini = False
                elif (dist_in > dist_out):
                    if (data_in[ind_in,k] >= thresh_dens): # input is ice layer                    
                        data_out[ind_out,k] += dist_out / (bound_out[ind_out + 1] - bound_out[ind_out])
                    ind_out += 1
                    dist = bound_out[ind_out]
                    flag_ini = True
                else:
                    if (data_in[ind_in,k] >= thresh_dens): # input is ice layer   
                        data_out[ind_out,k] += dist_in / (bound_out[ind_out + 1] - bound_out[ind_out])
                    ind_in += 1
                    ind_out += 1
                    dist = bound_in[ind_in,k]
                    flag_ini = True

    return data_out   
         
########################################################################################################################