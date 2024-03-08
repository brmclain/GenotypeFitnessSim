# Written by Brett Mclain. With segments trasnlated from perl code provided 
# by Professor Meisel.


import os
import numpy as np
import pandas as pd
import math
import sys
import time
import warnings
import argparse

from multiprocessing import Pool
from multiprocessing import Process

from datetime import datetime


def recursion_seln_Y1WY2(raw_counts_ref, fitness_ref):
    # Calculates a generation using the probability of producing offspring of a specific genotype and 
    # the fitness value of each genotype
    counts = []
    raw_counts = raw_counts_ref
    
    fitness = fitness_ref

    for i in range(18):
        counts.append(raw_counts[i] * fitness[i])   
    Sum1 = sum(counts)

    counts_norm = []
    for j in range(18):
        counts_norm.append(counts[j] / Sum1)
    
    f01 = counts_norm[0]
    f02 = counts_norm[1]
    f03 = counts_norm[2]
    f04 = counts_norm[3]
    f05 = counts_norm[4]
    f06 = counts_norm[5]
    f07 = counts_norm[6]
    m01 = counts_norm[7]
    m02 = counts_norm[8]
    m03 = counts_norm[9]
    m04 = counts_norm[10]
    m05 = counts_norm[11]
    m06 = counts_norm[12]
    m07 = counts_norm[13]
    m08 = counts_norm[14]
    m09 = counts_norm[15]
    m10 = counts_norm[16]
    m11 = counts_norm[17]
    
    after_mating = [0.5*f01*m01 + 0.25*f02*m01 + 0.125*f04*m01 + 0.25*f01*m02 + 0.125*f02*m02 + 0.0625*f04*m02 + 0.5*f01*m04 + 0.25*f02*m04 + 0.125*f04*m04 + 0.25*f01*m05 + 0.125*f02*m05 + 0.0625*f04*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.03125*f04*m06,
		            0.25*f02*m01 + 0.5*f03*m01 + 0.125*f04*m01 + 0.25*f05*m01 + 0.25*f01*m02 + 0.25*f02*m02 + 0.25*f03*m02 + 0.125*f04*m02 + 0.125*f05*m02 + 0.5*f01*m03 + 0.25*f02*m03 + 0.125*f04*m03 + 0.25*f02*m04 + 0.5*f03*m04 + 0.125*f04*m04 + 0.25*f05*m04 + 0.125*f02*m05 + 0.25*f03*m05 + 0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.0625*f04*m07,
		            0.125*f02*m02 + 0.25*f03*m02 + 0.0625*f04*m02 + 0.125*f05*m02 + 0.25*f02*m03 + 0.5*f03*m03 + 0.125*f04*m03 + 0.25*f05*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.03125*f04*m06 + 0.0625*f05*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.0625*f04*m07 + 0.125*f05*m07,
		            0.125*f04*m01 + 0.25*f05*m01 + 0.25*f06*m01 + 0.5*f07*m01 + 0.125*f04*m02 + 0.125*f05*m02 + 0.25*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f06*m03 + 0.25*f02*m04 + 0.5*f03*m04 + 0.25*f04*m04 + 0.5*f05*m04 + 0.25*f06*m04 + 0.5*f07*m04 + 0.125*f02*m05 + 0.25*f03*m05 + 0.125*f04*m05 + 0.25*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.125*f04*m06 + 0.125*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.125*f04*m07 + 0.125*f06*m07 + 0.5*f02*m08 + 1*f03*m08 + 0.25*f04*m08 + 0.5*f05*m08 + 0.25*f02*m09 + 0.5*f03*m09 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f01*m10 + 0.25*f02*m10 + 0.25*f03*m10 + 0.125*f04*m10 + 0.125*f05*m10 + 0.5*f01*m11 + 0.25*f02*m11 + 0.125*f04*m11,
		            0.0625*f04*m02 + 0.125*f05*m02 + 0.125*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f05*m03 + 0.25*f06*m03 + 0.5*f07*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.125*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.125*f04*m07 + 0.25*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.125*f02*m10 + 0.25*f03*m10 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.25*f02*m11 + 0.5*f03*m11 + 0.125*f04*m11 + 0.25*f05*m11,
		            0.125*f04*m04 + 0.25*f05*m04 + 0.25*f06*m04 + 0.5*f07*m04 + 0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f06*m07 + 0.25*f04*m08 + 0.5*f05*m08 + 0.5*f06*m08 + 1*f07*m08 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f06*m09 + 0.5*f07*m09 + 0.125*f04*m10 + 0.125*f05*m10 + 0.25*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f06*m11,
		            0.03125*f04*m06 + 0.0625*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.125*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f05*m11 + 0.25*f06*m11 + 0.5*f07*m11,
		            0.5*f01*m01 + 0.25*f02*m01 + 0.125*f04*m01 + 0.25*f01*m02 + 0.125*f02*m02 + 0.0625*f04*m02 + 0.25*f01*m05 + 0.125*f02*m05 + 0.0625*f04*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.03125*f04*m06,
		            0.25*f02*m01 + 0.5*f03*m01 + 0.125*f04*m01 + 0.25*f05*m01 + 0.25*f01*m02 + 0.25*f02*m02 + 0.25*f03*m02 + 0.125*f04*m02 + 0.125*f05*m02 + 0.5*f01*m03 + 0.25*f02*m03 + 0.125*f04*m03 + 0.125*f02*m05 + 0.25*f03*m05 + 0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.0625*f04*m07,
		            0.125*f02*m02 + 0.25*f03*m02 + 0.0625*f04*m02 + 0.125*f05*m02 + 0.25*f02*m03 + 0.5*f03*m03 + 0.125*f04*m03 + 0.25*f05*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.03125*f04*m06 + 0.0625*f05*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.0625*f04*m07 + 0.125*f05*m07,
		            0.125*f04*m01 + 0.25*f06*m01 + 0.0625*f04*m02 + 0.125*f06*m02 + 0.5*f01*m04 + 0.25*f02*m04 + 0.25*f04*m04 + 0.25*f06*m04 + 0.25*f01*m05 + 0.125*f02*m05 + 0.125*f04*m05 + 0.125*f06*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.0625*f04*m06 + 0.0625*f06*m06 + 1*f01*m08 + 0.5*f02*m08 + 0.25*f04*m08 + 0.5*f01*m09 + 0.25*f02*m09 + 0.125*f04*m09 + 0.25*f01*m10 + 0.125*f02*m10 + 0.0625*f04*m10,
		            0.125*f04*m01 + 0.25*f06*m01 + 0.0625*f04*m02 + 0.125*f06*m02 + 0.25*f01*m05 + 0.125*f02*m05 + 0.125*f04*m05 + 0.125*f06*m05 + 0.125*f01*m06 + 0.0625*f02*m06 + 0.0625*f04*m06 + 0.0625*f06*m06 + 0.5*f01*m09 + 0.25*f02*m09 + 0.125*f04*m09 + 0.25*f01*m10 + 0.125*f02*m10 + 0.0625*f04*m10,
		            0.125*f04*m01 + 0.25*f05*m01 + 0.25*f06*m01 + 0.5*f07*m01 + 0.125*f04*m02 + 0.125*f05*m02 + 0.25*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f06*m03 + 0.125*f02*m05 + 0.25*f03*m05 + 0.125*f04*m05 + 0.25*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.125*f01*m06 + 0.125*f02*m06 + 0.125*f03*m06 + 0.125*f04*m06 + 0.125*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.25*f01*m07 + 0.125*f02*m07 + 0.125*f04*m07 + 0.125*f06*m07 + 0.25*f02*m09 + 0.5*f03*m09 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f01*m10 + 0.25*f02*m10 + 0.25*f03*m10 + 0.125*f04*m10 + 0.125*f05*m10 + 0.5*f01*m11 + 0.25*f02*m11 + 0.125*f04*m11,
		            0.0625*f04*m02 + 0.125*f05*m02 + 0.125*f06*m02 + 0.25*f07*m02 + 0.125*f04*m03 + 0.25*f05*m03 + 0.25*f06*m03 + 0.5*f07*m03 + 0.0625*f02*m06 + 0.125*f03*m06 + 0.0625*f04*m06 + 0.125*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.125*f02*m07 + 0.25*f03*m07 + 0.125*f04*m07 + 0.25*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.125*f02*m10 + 0.25*f03*m10 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.25*f02*m11 + 0.5*f03*m11 + 0.125*f04*m11 + 0.25*f05*m11,
		            0.125*f04*m04 + 0.25*f06*m04 + 0.0625*f04*m05 + 0.125*f06*m05 + 0.03125*f04*m06 + 0.0625*f06*m06 + 0.25*f04*m08 + 0.5*f06*m08 + 0.125*f04*m09 + 0.25*f06*m09 + 0.0625*f04*m10 + 0.125*f06*m10,
		            0.0625*f04*m05 + 0.125*f06*m05 + 0.03125*f04*m06 + 0.0625*f06*m06 + 0.125*f04*m09 + 0.25*f06*m09 + 0.0625*f04*m10 + 0.125*f06*m10,
		            0.0625*f04*m05 + 0.125*f05*m05 + 0.125*f06*m05 + 0.25*f07*m05 + 0.0625*f04*m06 + 0.0625*f05*m06 + 0.125*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f06*m07 + 0.125*f04*m09 + 0.25*f05*m09 + 0.25*f06*m09 + 0.5*f07*m09 + 0.125*f04*m10 + 0.125*f05*m10 + 0.25*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f06*m11,
		            0.03125*f04*m06 + 0.0625*f05*m06 + 0.0625*f06*m06 + 0.125*f07*m06 + 0.0625*f04*m07 + 0.125*f05*m07 + 0.125*f06*m07 + 0.25*f07*m07 + 0.0625*f04*m10 + 0.125*f05*m10 + 0.125*f06*m10 + 0.25*f07*m10 + 0.125*f04*m11 + 0.25*f05*m11 + 0.25*f06*m11 + 0.5*f07*m11]
    
    Sum2 = sum(after_mating)
    after_mating_std = []
    for k in range(18):
        after_mating_std.append(after_mating[k]/Sum2)
        
    return(after_mating_std)
    
def genotype_fitness_Y1WY2(input_array, dominance):
    
    YM_fem = input_array[0]
    YM_male = input_array[1]
    W_fem = input_array[2]
    W_male = input_array[3]
    Y2_male = input_array[4]

    W_fem_dom = 0
    W_male_dom = 0
    YM_fem_dom = 0
    YM_male_dom = 0

    #========== SECTION GENOTYPE_FITNESS A START ==========#
    # Chooses dom or "h" values based on the dominance type
    if(dominance == 'additive'):
        W_fem_dom = 0.5
        W_male_dom = 0.5
        YM_fem_dom = 0.5
        YM_male_dom = 0.5
        
    if(dominance == 'dominant'):
        if(W_fem <= 0):
            W_fem_dom = 1
        if(W_male <= 0):
            W_male_dom = 1
        if(YM_fem <= 0):
            YM_fem_dom = 1
        if(YM_male <= 0):
            YM_male_dom = 1 
            
    if(dominance == 'recessive'):
        if(W_fem >= 0):
            W_fem_dom = 1
        if(W_male >= 0):
            W_male_dom = 1
        if(YM_fem >= 0):
            YM_fem_dom = 1
        if(YM_male >= 0):
            YM_male_dom = 1 
            
    if dominance == 'random':
        np.random.seed()
        rand_dominance = np.random.uniform(size = 4, low = 0, high = 1)  # array of size 4
        YM_fem_dom = rand_dominance[0]
        YM_male_dom = rand_dominance[1]
        W_fem_dom = rand_dominance[2]
        W_male_dom = rand_dominance[3]
    #========== SECTION GENOTYPE_FITNESS A END ==========#
    
    
    #========== SECTION GENOTYPE_FITNESS B START ==========#
    # Calculates fitness values for each allele combination
    # Y1 female    
    if(YM_fem < 0):
        X1X1_fem = 1
        X1Y1_fem = 1 + YM_fem_dom*YM_fem
        Y1Y1_fem = 1 + YM_fem
    else:
        X1X1_fem = 1 - YM_fem
        X1Y1_fem = 1 - YM_fem_dom*YM_fem
        Y1Y1_fem = 1
        
    # W female
    if(W_fem < 0):
        ZZ_fem = 1
        ZW_fem = 1 + W_fem_dom*W_fem
        WW_fem = 1 + W_fem
    else:
        ZZ_fem = 1 - W_fem
        ZW_fem = 1 - W_fem_dom*W_fem
        WW_fem = 1
    
    # Y1 males
    if(YM_male < 0):
        X1X1_male = 1
        X1Y1_male = 1 + YM_fem_dom*YM_male
        Y1Y1_male = 1 + YM_male
    else:
        X1X1_male = 1 - YM_male
        X1Y1_male = 1 - YM_male_dom*YM_male
        Y1Y1_male = 1
    
    # W males
    if(W_male < 0):
        ZZ_male = 1
        ZW_male = 1 + W_male_dom*W_male
        WW_male = 1 + W_male
    else:
        ZZ_male = 1 - W_male
        ZW_male = 1 - W_male_dom*W_male
        WW_male = 1
    
    # Y2 males
    if(Y2_male < 0):
        X2X2_male = 1
        X2Y2_male = 1 + Y2_male
    else:
        X2X2_male = 1 - Y2_male
        X2Y2_male = 1 
    #========== SECTION GENOTYPE_FITNESS B END ==========#
    
    
    #========== SECTION GENOTYPE_FITNESS C START ==========#
    # Calculated the fitness of each genotype based on allele combinations and normalizes the results.    
    genotype_fit = [X1X1_fem * ZZ_fem ,		            #0  female	f01	X1X1	ZZ	X2X2
		            X1X1_fem * ZW_fem ,					#1  female	f02	X1X1	ZW	X2X2
		            X1X1_fem * WW_fem ,					#2  female	f03	X1X1	WW	X2X2
		            Y1Y1_fem * ZW_fem ,					#3  female	f04	X1Y1	ZW	X2X2
		            X1Y1_fem * WW_fem , 				#4  female	f05	X1Y1	WW	X2X2
		            Y1Y1_fem * ZW_fem ,					#5  female	f06	Y1Y1	ZW	X2X2
		            Y1Y1_fem * WW_fem , 				#6  female	f07	Y1Y1	WW	X2X2
		            X1X1_male * ZZ_male * X2Y2_male ,	#7  male	m01	X1X1	ZZ	X2Y2
		            X1X1_male * ZW_male * X2Y2_male ,	#8  male	m02	X1X1	ZW	X2Y2
		            X1X1_male * WW_male * X2Y2_male ,	#9  male	m03	X1X1	WW	X2Y2
		            Y1Y1_male * ZZ_male * X2X2_male ,	#10 male	m04	X1Y1	ZZ	X2X2
		            X1Y1_male * ZZ_male * X2Y2_male ,	#11 male	m05	X1Y1	ZZ	X2Y2
		            X1Y1_male * ZW_male * X2Y2_male ,	#12 male	m06	X1Y1	ZW	X2Y2
		            X1Y1_male * WW_male * X2Y2_male ,	#13 male	m07	X1Y1	WW	X2Y2
		            Y1Y1_male * ZZ_male * X2X2_male ,	#14 male	m08	Y1Y1	ZZ	X2X2
		            Y1Y1_male * ZZ_male * X2Y2_male ,	#15 male	m09	Y1Y1	ZZ	X2Y2
		            Y1Y1_male * ZW_male * X2Y2_male ,	#16 male	m10	Y1Y1	ZW	X2Y2
		            Y1Y1_male * WW_male * X2Y2_male]	#17 male	m11	Y1Y1	WW	X2Y2]
    
    genotype_norm = genotype_fit
    max_female = max(genotype_norm[0:7])
    max_male = max(genotype_norm[7:18])
    for i in range(0,18):
        if i <= 6:
            genotype_norm[i] = genotype_fit[i] / max_female
        else:
            genotype_norm[i] = genotype_fit[i] / max_male       
    dom_array = [YM_fem_dom, YM_male_dom, W_fem_dom, W_male_dom]
    #========== SECTION GENOTYPE_FITNESS C END ==========#
    
    return(genotype_norm, dom_array)
 
def binary_search(arr, high):
    # This function will search to find where each genotype reaches its point of equilibrium
    # Our high is the number of generations
    low = 0
    high = high - 1
    while(low <= high):
        mid = int((low + high)/2)
        if abs(arr[mid + 1] - arr[mid - 1]) > 1e-6 and abs(arr[mid + 1] - arr[mid - 1]) < 1e-5 : # If the difference between 2 almost consecutive points is 1e-6, then we know we are reaching the point where the asymptote begins
            return mid
        elif abs(arr[mid + 1] - arr[mid - 1]) > 1e-5: 
            low = mid + 1
        elif abs(arr[mid + 1] - arr[mid - 1]) < 1e-6:
            high = mid - 1
    
    return 0 
    
def sex_ratio_Y1WY2(freqs):
    # Calculates the frequencies of sexes based on the frequency of each genotype
    sum_freq = sum(freqs)
    female = (freqs[0] + freqs[1] + freqs[2] + freqs[3] + freqs[4] + freqs[5] + freqs[6] + freqs[7])/sum_freq
    male = 1 - female
    return(female, male)

def allele_freqs_Y1WY2(freqs):
    # Calculates the frequencies of alleles based on the frequency they occur in each genotype.
    sum_freq = sum(freqs)
    
    X_freq = (freqs[0] + freqs[1] + freqs[2] + freqs[7] + freqs[8] + freqs[9]
              + 0.5*(freqs[3] + freqs[4] + freqs[10] + freqs[11] + freqs[12] + freqs[13]))/sum_freq
    Y_freq = 1 - X_freq
    
    III_freq = (freqs[0] + freqs[1] + freqs[2] + freqs[3] + freqs[4] + freqs[5] + freqs[6] + freqs[10] + freqs[14]
                + 0.5*(freqs[7] + freqs[8] + freqs[9] + freqs[11] + freqs[12] + freqs[13] + freqs[15] + freqs[16] + freqs[17]))
    IIIM_freq = 1 - III_freq
    
    tra_freq = (freqs[0] + freqs[7] + freqs[10] + freqs[11] + freqs[14] + freqs[15]
                + 0.5*(freqs[1] + freqs[3] + freqs[5] + freqs[8] + freqs[12] + freqs[16]))
    traD_freq = 1 - tra_freq
    
    return(X_freq, Y_freq, III_freq, IIIM_freq, traD_freq, tra_freq)

def simu_Y1WY2(generations, dominance, initial_freq, allele_range, switches="1111111"):

    # Sets a random distribution for 5 allele fitness values
    np.random.seed()
    allele_fitness = np.random.uniform(size = 5, low = allele_range*(-1), high = allele_range)
    
    # Creates fitness value for each genotype
    fitness_array, dom_result = genotype_fitness_Y1WY2(allele_fitness, dominance)  
    
    # Genotype fitness and frequencies are then passed into a function to simulate breeding probabilities over
    # a specified number of generations. 
    simulation_results = []
    simulation_results.append(initial_freq)
    for a in range(generations):                                
        result_recursion = recursion_seln_Y1WY2(simulation_results[a-1], fitness_array)   # Takes genotoype frequencies and uses same genotype fitness array
        simulation_results.append(result_recursion)  
    
    # Returns the number of generations it took to reach a level of equalibrium.
    geno_sim_results = np.array(simulation_results) # convert simulation results into a numpy array for slicing
    geno_asymp = []
    if switches[6] == '1':                              
        for genotypes in range(len(simulation_results[0])):         # for 1 through 18
            geno_asymp1 = binary_search(geno_sim_results[0:, genotypes], generations)   # binary search takes each genotype's data and the number of generations that were performed
            geno_asymp.append(geno_asymp1)
    # geno_asymp will contain all 18 points of equilibrium (the generation where equilibrium occurs)
    
    # Data is collected and placed into a 2D array to be returned.
    final_sim_result = simulation_results[generations]                          
    allele_freq = allele_freqs_Y1WY2(final_sim_result)
    sex_ratios = sex_ratio_Y1WY2(final_sim_result)
    Desired_data = []
    if switches[0] == '1':
        Desired_data.append(allele_fitness)
    if switches[1] == '1':
        Desired_data.append(fitness_array)
    if switches[2] == '1':
        Desired_data.append(dom_result)
    if switches[3] == '1':
        Desired_data.append(final_sim_result)
    if switches[4] == '1':
        Desired_data.append(allele_freq)
    if switches[5] == '1':
        Desired_data.append(sex_ratios)
    if switches[6] == '1':
        Desired_data.append(geno_asymp)
        
    return Desired_data

def main_Y1WY2(JobID, simulations, generations, processes, dominance, initial_freq, allele_range, switches):
    
    start = time.time()
    cur_path = os.path.abspath(os.getcwd())
    output = "data_" + str(JobID) + ".txt"
    modif_path = os.path.join(cur_path, "OutputArchive")
    try:
        os.mkdir(modif_path)
    except:
        pass
    
    #==========SECTION MAIN A START==========#
    # Creation of dictionaries to organize the resulting data from the simulation into a data frame.
    allele_fitness_dict = {'YM_fem_fitness': [], 'YM_male_fitness': [], 'traD_fitness': [], 'tra_fitness': [],  'IIIM_male_fitness': []}
    fitness_array_dict = {'f1_fit': [], 'f2_fit': [], 'f3_fit': [], 'f4_fit': [], 'f5_fit': [], 'f6_fit': [], 'f7_fit': [],
                        'm1_fit': [], 'm2_fit': [], 'm3_fit': [], 'm4_fit': [], 'm5_fit': [], 'm6_fit': [], 'm7_fit': [], 'm8_fit': [], 'm9_fit': [], 'm10_fit': [], 'm11_fit': []}
    dom_result_dict = {'YM_fem_h': [], 'YM_male_h': [], 'W_fem_h': [], 'W_male_h': []}
    final_sim_result_dict = {'f1': [], 'f2': [], 'f3': [], 'f4': [], 'f5': [], 'f6': [], 'f7': [],
                        'm1': [], 'm2': [], 'm3': [], 'm4': [], 'm5': [], 'm6': [], 'm7': [], 'm8': [], 'm9': [], 'm10': [], 'm11': []}
    allele_freq_dict = {'X_freq': [], 'Y_freq': [], 'III_freq': [], 'IIIM_freq': [], 'traD_freq': [], 'tra_freq': []}
    sex_ratios_dict = {'female' : [], 'male': []}
    genotype_equil_dict = {'f1_equil': [], 'f2_equil': [], 'f3_equil': [], 'f4_equil': [], 'f5_equil': [], 'f6_equil': [], 'f7_equil': [],
                        'm1_equil': [], 'm2_equil': [], 'm3_equil': [], 'm4_equil': [], 'm5_equil': [], 'm6_equil': [], 'm7_equil': [], 'm8_equil': [], 'm9_equil': [], 'm10_equil': [], 'm11_equil': []}
    
    df_allele_fitness = pd.DataFrame(allele_fitness_dict)
    df_fitness_array = pd.DataFrame(fitness_array_dict)
    df_dom_result = pd.DataFrame(dom_result_dict)
    df_final_sim_result = pd.DataFrame(final_sim_result_dict)
    df_allele_freq = pd.DataFrame(allele_freq_dict)
    df_sex_ratios = pd.DataFrame(sex_ratios_dict)
    df_genotype_equil = pd.DataFrame(genotype_equil_dict)
    #==========SECTION MAIN A END==========#
    
    
    #==========SECTION MAIN B START==========#
    # Simulations are run in parallel using starmap. Output data is sorted using "switches" and combined into
    # a 2D array. Then exported into a .txt file.
    
    # The number of simulations that can run in parallel before the results are exported to a file are capped
    # to avoid "data" (The varaible holding the result of starmap) becoming too large. 
    
    # Varaible "section" should be adjusted to be as large as possible without making "data" so large
    # that it cannot be contained in maina memory. Should this occur, processing times will increase dramatically.
    section = 25000
    total_sims = simulations
    p = Pool(processes = def_processes)
    while(simulations>0):
        if(simulations > section):
            data = list(p.starmap(simu_Y1WY2, [(generations, dominance, initial_freq, allele_range, switches) for _ in range (section)]))
        else:
            data = list(p.starmap(simu_Y1WY2, [(generations, dominance, initial_freq, allele_range, switches) for _ in range (simulations)]))
        
        for i in range(len(data)):
            j = 0
            if switches[0] == '1':
                df_allele_fitness.loc[i] = data[i][j]
                j += 1
            if switches[1] == '1':
                df_fitness_array.loc[i] = data[i][j]
                j += 1
            if switches[2] == '1':
                df_dom_result.loc[i] = data[i][j]
                j += 1
            if switches[3] == '1':
                df_final_sim_result.loc[i] = data[i][j]
                j += 1
            if switches[4] == '1':
                df_allele_freq.loc[i] = data[i][j]
                j += 1
            if switches[5] == '1':
                df_sex_ratios.loc[i] = data[i][j]
                j += 1
            if switches[6] == '1':
                df_genotype_equil.loc[i] = data[i][j]
        data_frames = []
        
        if switches[0] == '1':
            data_frames.append(df_allele_fitness)
        if switches[1] == '1':
            data_frames.append(df_fitness_array)
        if switches[2] == '1':
            data_frames.append(df_dom_result)
        if switches[3] == '1':
            data_frames.append(df_final_sim_result)
        if switches[4] == '1':
            data_frames.append(df_allele_freq)
        if switches[5] == '1':
            data_frames.append(df_sex_ratios)
        if switches[6] == '1':
            data_frames.append(df_genotype_equil)
            
        if len(data_frames) > 0:      
            final_data = pd.concat(data_frames, axis = 1)
        else: 
            final_data = pd.DataFrame()
            
        try:
            final_data.to_csv(os.path.join(modif_path, output), sep = '\t', index = None, mode = 'x')
        except:
            final_data.to_csv(os.path.join(modif_path, output), sep = '\t', index = None, mode = 'a', header=False)   
        
        simulations -= section
        if(simulations <= 0):
            break 
    
    p.close()
    p.join()
    #==========SECTION MAIN B END==========#
    
    
    #==========SECTION MAIN C START==========#
    # Parameters used in the simulation are saved along with the jobID and exported into a file containing 
    # a history of previous jobs.
    end = time.time()
    total_time = "{:.2f}".format(end-start)
    print(end - start, "seconds")
    
    with open(os.path.join(modif_path, "ParameterLog.txt"), 'a') as f:
        f.write(str(JobID))
        f.write(" Precessors:")
        f.write(str(def_processes))
        f.write(" Time:")
        f.write(total_time)
        f.write("s Simulations:")
        f.write(str(total_sims))
        f.write(" Generations:")
        f.write(str(generations))
        f.write(" Dominance:")
        f.write(dominance)
        f.write(" Allele Range:")
        f.write(str(allele_range))
        f.write(" Switch Board:")
        f.write(switches)
        f.write(" Initial Freq:[")
        f.write(" ".join(map(str,initial_freq)))
        f.write("] \n")  
    #==========SECTION MAIN C END==========#   


#The only thing here should be command line argument handling and the main function call
if __name__ == "__main__":
    
    #========== Section GLOBAL A START ==========#
    # Placing command line arguments into varaibles
    parser = argparse.ArgumentParser(description = "Enter the JobID first, number of simulations second, and the number of generations third")
    parser.add_argument('SJobID', type=int)
    parser.add_argument('simulations', type = int)
    parser.add_argument('generations', type = int)
    parser.add_argument('def_processes', type = int)
    parser.add_argument('dominance', type = str) 
    parser.add_argument('-o', '--control_file', required=False)
    # Line below defines the upper limit and multiplies it by -1 to create a mirror lower limit for the creation of a uniform distribution in allele_fitness.
    # Defaults to -1 to 1 if not included or given a bad input.
    parser.add_argument('-a', '--allele_range', required=False, type=float)
    #Line below specifies which data catagories to return into the output file with input of 1-'return' or 0-'dont return'.
    #Defaults to all values of "1" if input is bad or not given. EX: "1111111"
    #Order of input: [0] allele_fitness
    #                [1] fitness_array
    #                [2] dom_result
    #                [3] final_sim_result
    #                [4] allele_freq
    #                [5] sex_ratios
    #                [6] genotype_equil
    parser.add_argument('-s', '--switches', required=False, type=str)
    
    args = parser.parse_args()
    JobID = args.SJobID
    simulations = args.simulations
    generations = args.generations
    def_processes = args.def_processes
    dominance = args.dominance
    control_file = args.control_file
    allele_range = args.allele_range
    switches = args.switches
    #========== Section GLOBAL A END ==========#
    

    #========== SECTION GLOBAL B START ==========#
    # Each command line argument is checked if it is a usable. If their is a problem with the value, they are 
    # replaced with default values or stops the program via and error.
    if control_file == None:                # Checks first to see if file is even an argument, if not, default to initial_freq
        initial_freq = [0.05556, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 
        0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556] # PUT FREQUENCIES HERE
    else:                                   # If it is argument, do the following
        try:                    
            control_file = open(control_file, 'r')      # Open the file, if not found in directory, again default to initial_freq
        except FileNotFoundError:                       # Catches the exception if FileNotFound
            initial_freq = [0.05556, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 
        0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556] # PUT FREQUENCIES HERE
            pass
        else:                                           # If found, parse file and load its contents
            new_control = control_file.read()
            new_control = new_control.replace('[', '')
            new_control = new_control.replace(']', '')
            new_control = new_control.split(', ')
            parsed_input = list(map(lambda s: float(s), new_control))# Creates a lambda function that takes each string instance "parsed_input" and converts it into a float.
            # map allows lambda to iterate over indexes in parsed_input to generate a float output. Then a list of all the outputs are placed in a the "parsed_input" list
            initial_freq = parsed_input
            pass

    #assert(simulations > 1), "Simulations must be greater than zero"
    assert(generations > 1), "Generation must be greater than zero"
    if dominance != 'dominant' and dominance != 'recessive' and dominance != 'random' and dominance != 'additive':
        print("WARNING: Since dominance selected is invalid, dominance has defaulted to random")  # prints warning but does not crash program
        dominance = 'random'    # If invalid scenario selected, code will default to random scenario

    if switches is not None:
        if len(switches) == 7:
            switches = list(switches)
            for i in range(7):
                if (switches[i] != '1') and (switches[i] != '0'):
                    switches[i] = '1'
        else:
            switches = "1111111"
    else:
        switches = "1111111"
        
    if type(allele_range).__name__ == "float":
        if (allele_range>0) and (allele_range<=1):
            pass
        else:
            allele_range = 1
    else:
        allele_range = 1
    #========== SECTION GLOBAL B END ==========#
    
    
    main_Y1WY2(JobID, simulations, generations, def_processes, dominance, initial_freq, allele_range, switches)
    