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

# Written in Python by Martin Elias. Perl code taken from Professor Meisel and edited by Martin Elias
# Edited by Brett Mclain Jan 2023
# This specific code simulates populations where the Proto-Y chromosomes have an ADDITIVE fitness effect


#CHANGES MADE
#   Removed command line augment for output file and replaced with dynamic output file naming scheme and Book keeping of all setting appended to master text file.
#   The total number of simulations now executes in batches of size "section" to avoid creating a dataframe too large to fit in memory
#   Command line argument "switches" added to allow specific data to be returned.
#   Small reorganization to increase readability and prevent unwanted excution when imported

# This function should pass the lists by reference (A list is mutable in Python, a tuple is immutable)
def recursion_seln(raw_counts_ref, fitness_ref):    # pass by reference, not copy
    counts = []                                     # declare counts
    raw_counts = raw_counts_ref

    fitness = fitness_ref

    for i in range(18):
        counts.append(raw_counts[i] * fitness[i])   # append to counts since it is empty, this will append in increasing index order (0,1,2...)    
    Sum1 = sum(counts)

    counts_norm = []
    for j in range(18):
        counts_norm.append(counts[j] / Sum1)
    
    #f10 is placed in the right place here (used to be m9)
    f1 = counts_norm[0]
    f2 = counts_norm[1]
    f3 = counts_norm[2]
    f4 = counts_norm[3]
    f5 = counts_norm[4]
    f6 = counts_norm[5]
    f7 = counts_norm[6]
    f8 = counts_norm[7]
    f9 = counts_norm[8]
    f10 = counts_norm[9]
    m1 = counts_norm[10]
    m2 = counts_norm[11]
    m3 = counts_norm[12]
    m4 = counts_norm[13]
    m5 = counts_norm[14]
    m6 = counts_norm[15]
    m7 = counts_norm[16]
    m8 = counts_norm[17]

    # To find the probability values cross the alleles of the mentioned parents (Ex. f9*m1). Once crossed, determine the probability 
    # for the genotype that matches the female or male index the cross is in.
    # (Ex. .25 of genotypes match f1's genotype in a cross between f2 and m1) 

    # This is where the genotype frequency is being calculated
    after_mating = [f1*m1*.5 + f1*m3*.5 + f1*m4*.25 + f2*m1*.25 + f2*m3*.25 + f2*m4*.125 + f3*m1*.25 + f3*m3*.125 + 
    f3*m4*.0625  + f5*m1*.125 + f5*m3*.125 + f5*m4*.0625 + f6*m1*.0625 + f6*m3*.0625 + f6*m4*.03125,
    f2*m1*.25 + f2*m3*.25 + f2*m4*.125 + f3*m1*.125 + f3*m3*.125 + f3*m4*.0625 + f5*m1*.125 +
    f5*m3*.125 + f5*m4*.0625 + f6*m1*.0625 + f6*m3*.0625 + f6*m4*.03125,
    f2*m1*.25 + f2*m2*.5 + f2*m4*.125 + f2*m5*.25 + f3*m1*.25 + f3*m2*.25 + f3*m3*.125 + f3*m4*.125 + 
    f3*m5*.125 + f4*m1*.25 + f4*m3*.25 + f4*m4*.125 + f5*m1*.125 + f5*m2*.25 + f5*m4*.0625 + f5*m5*.125 + f6*m1*.125
    + f6*m2*.125 + f6*m3*.0625 + f6*m4*.0625 + f6*m5*.0625 + f7*m1*.125 + f7*m3*.125 + f7*m4*.0625,
    f3*m1*.125 + f3*m2*.25 + f3*m4*.0625 + f3*m5*.125 + f4*m1*.25 + f4*m2*.5 + f4*m4*.125 + f4*m5*.25 + f6*m1*.0625
    + f6*m2*.125 + f6*m4*.03125 + f6*m5*.0625 + f7*m1*.125 + f7*m2*.25 + f7*m4*.0625 + f7*m5*.125,
    f2*m3*.25 + f2*m4*.125 + f2*m6*.5 + f2*m7*.25 + f3*m3*.125 + f3*m4*.0625 + f3*m6*.25 + f3*m7*.125 + f5*m1*.125
    + f5*m3*.25 + f5*m4*.125 + f5*m6*.25 + f5*m7*.125 + f6*m1*.0625 + f6*m3*.125 + f6*m4*.0625 + f6*m6*.125 + f6*m7*.0625
    + f8*m1*.25 + f8*m3*.25 + f8*m4*.125 + f9*m1*.125 + f9*m3*.125 + f9*m4*.0625,
    f2*m4*.125 + f2*m5*.25 + f2*m7*.25 + f2*m8*.5 + f3*m3*.125 + f3*m4*.125 + f3*m5*.125 + f3*m6*.25 + f3*m7*.25 + 
    f3*m8*.25 + f4*m3*.25 + f4*m4*.125 + f4*m6*.5 + f4*m7*.25 + f5*m1*.125 + f5*m2*.25 + f5*m4*.125 + f5*m5*.25 + 
    f5*m7*.125 + f5*m8*.25 + f6*m1*.125 + f6*m2*.125 + f6*m3*.125 + f6*m4*.125 + f6*m5*.125 + f6*m6*.125 + f6*m7*.125
    + f6*m8*.125 + f7*m1*.125 + f7*m3*.25 + f7*m4*.125 + f7*m6*.25 + f7*m7*.125 + f8*m1*.25 + f8*m2*.5 + f8*m4*.125
    + f8*m5*.25 + f9*m1*.25 + f9*m2*.25 + f9*m3*.125 + f9*m4*.125 + f9*m5*.125 + f10*m1*.25 + f10*m3*.25 + f10*m4*.125,
    f3*m4*.0625 + f3*m5*.125 + f3*m7*.125 + f3*m8*.25 + f4*m4*.125 + f4*m5*.25 + f4*m7*.25 + f4*m8*.5 + f6*m1*.0625
    + f6*m2*.125 + f6*m4*.0625 + f6*m5*.125 + f6*m7*.0625 + f6*m8*.125 + f7*m1*.125 + f7*m2*.25 + f7*m4*.125 + f7*m5*.25
    + f7*m7*.125 + f7*m8*.25 + f9*m1*.125 + f9*m2*.25 + f9*m4*.0625 + f9*m5*.125 + f10*m1*.25 + f10*m2*.5 + f10*m4*.125 + f10*m5*.25,
    f5*m3*.125 + f5*m4*.0625 + f5*m6*.25 + f5*m7*.125 + f6*m3*.0625 + f6*m4*.03125 + f6*m6*.125 +f6*m7*.0625 + f8*m3*.25
    + f8*m4*.125 + f8*m6*.5 + f8*m7*.25 + f9*m3*.125 + f9*m4*.0625 + f9*m6*.25 + f9*m7*.125,
    f5*m4*.0625 + f5*m5*.125 + f5*m7*.125 + f5*m8*.25 + f6*m3*.0625 + f6*m4*.0625 + f6*m5*.0625 + f6*m6*.125 + f6*m7*.125
    + f6*m8*.125 + f7*m3*.125 + f7*m4*.0625 + f7*m6*.25 + f7*m7*.125 + f8*m4*.125 + f8*m5*.25 + f8*m7*.25 + f8*m8*.5 + f9*m3*.125
    + f9*m4*.125 + f9*m5*.125 + f9*m6*.25 + f9*m7*.25 + f9*m8*.25 + f10*m3*.25 + f10*m4*.125 + f10*m6*.5 + f10*m7*.25,
    f1*m1*.5 + f1*m2*1 + f1*m4*.25 + f1*m5*.5 + f2*m1*.25 + f2*m2*.5 + f2*m4*.125 + f2*m5*.25 + f3*m1*.125 + f3*m2*.25
    + f3*m3*.125 + f3*m4*.125 + f3*m5*.125 + f4*m1*.25 + f4*m3*.25 + f4*m4*.125 + f5*m1*.125 + f5*m2*.25 + f5*m4*.0625
    + f5*m5*.125 + f6*m1*.125 + f6*m2*.125 + f6*m3*.0625 + f6*m4*.0625 + f6*m5*.0625 + f7*m1*.125 + f7*m3*.125 + f7*m4*.0625,
    f6*m4*.03125 + f6*m5*.0625 + f6*m7*.0625 + f6*m8*.125 + f7*m4*.0625 + f7*m5*.125 + f7*m7*.125 + f7*m8*.25 + f9*m4*.0625
    + f9*m5*.125 + f9*m7*.125 + f9*m8*.25 + f10*m4*.125 + f10*m5*.25 + f10*m7*.25 + f10*m8*.5,
    f3*m1*.125 + f3*m2*.25 + f3*m4*.0625 + f3*m5*.125 + f4*m1*.25 + f4*m2*.5 + f4*m4*.125 + f4*m5*.25 + f6*m1*.0625
    + f6*m2*.125 + f6*m4*.03125 + f6*m5*.0625 + f7*m1*.125 + f7*m2*.25 + f7*m4*.0625 + f7*m5*.125,
    f1*m3*.5 + f1*m4*.25 + f1*m6*1 + f1*m7*.5 + f2*m3*.25 + f2*m4*.125 + f2*m6*.5 + f2*m7*.25 + f3*m3*.125 + f3*m4*.0625
    + f3*m6*.25 + f3*m7*.125 + f5*m1*.125 + f5*m3*.25 + f5*m4*.125 + f5*m6*.25 + f5*m7*.125 + f6*m1*.0625 + f6*m3*.125
    + f6*m4*.0625 + f6*m6*.125 + f6*m7*.0625 + f8*m1*.25 + f8*m3*.25 + f8*m4*.125 + f9*m1*.125 + f9*m3*.125 + f9*m4*.0625,
    f1*m4*.25 + f1*m5*.5 + f1*m7*.5 + f1*m8*1 + f2*m4*.125 + f2*m5*.25 + f2*m7*.25 + f2*m8*.5 + f3*m3*.125 + f3*m4*.125
    + f3*m5*.125 + f3*m6*.25 + f3*m7*.25 + f3*m8*.25 + f4*m3*.25 + f4*m4*.125 + f4*m6*.5 + f4*m7*.25 + f5*m1*.125
    + f5*m2*.25 + f5*m4*.125 + f5*m5*.25 + f5*m7*.125 + f5*m8*.25 + f6*m1*.125 + f6*m2*.125 + f6*m3*.125 + f6*m4*.125
    + f6*m5*.125 + f6*m6*.125 + f6*m7*.125 + f6*m8*.125 + f7*m1*.125 + f7*m3*.25 + f7*m4*.125 + f7*m6*.25 + f7*m7*.125
    + f8*m1*.25 + f8*m2*.5 + f8*m4*.125 + f8*m5*.25 + f9*m1*.25 + f9*m2*.25 + f9*m3*.125 + f9*m4*.125 + f9*m5*.125
    + f10*m1*.25 + f10*m3*.25 + f10*m4*.125,
    f3*m4*.0625 + f3*m5*.125 + f3*m7*.125 + f3*m8*.25 + f4*m4*.125 + f4*m5*.25 + f4*m7*.25 + f4*m8*.5 + f6*m1*.0625 + f6*m2*.125
    + f6*m4*.0625 + f6*m5*.125 + f6*m7*.0625 + f6*m8*.125 + f7*m1*.125 + f7*m2*.25 + f7*m4*.125 + f7*m5*.25 + f7*m7*.125
    + f7*m8*.25 + f9*m1*.125 + f9*m2*.25 + f9*m4*.0625 + f9*m5*.125 + f10*m1*.25 + f10*m2*.5 + f10*m4*.125 + f10*m5*.25,
    f5*m3*.125 + f5*m4*.0625 + f5*m6*.25 + f5*m7*.125 + f6*m3*.0625 + f6*m4*.03125 + f6*m6*.125 + f6*m7*.0625 + f8*m3*.25
    + f8*m4*.125 + f8*m6*.5 + f8*m7*.25 + f9*m3*.125 + f9*m4*.0625 + f9*m6*.25 + f9*m7*.125,
    f5*m4*.0625 + f5*m5*.125 + f5*m7*.125 + f5*m8*.25 + f6*m3*.0625 + f6*m4*.0625 + f6*m5*.0625 + f6*m6*.125 + f6*m7*.125
    + f6*m8*.125 + f7*m3*.125 + f7*m4*.0625 + f7*m6*.25 + f7*m7*.125 + f8*m4*.125 + f8*m5*.25 + f8*m7*.25 + f8*m8*.5
    + f9*m3*.125 + f9*m4*.125 + f9*m5*.125 + f9*m6*.25 + f9*m7*.25 + f9*m8*.25 + f10*m3*.25 + f10*m4*.125 + f10*m6*.5 + f10*m7*.25,
    f6*m4*.03125 + f6*m5*.0625 + f6*m7*.0625 + f6*m8*.125 + f7*m4*.0625 + f7*m5*.125 + f7*m7*.125 + f7*m8*.25 + f9*m4*.0625
    + f9*m5*.125 + f9*m7*.125 + f9*m8*.25 + f10*m4*.125 + f10*m5*.25 + f10*m7*.25 + f10*m8*.5]

    
    
    Sum2 = sum(after_mating)
    after_mating_std = []       # Declare the "after_mating_std" list
    for k in range(18):
        after_mating_std.append(after_mating[k] / Sum2) 

    # after_mating_std = 18 genotype values
    

    return(after_mating_std)


# Following function determines the genotype fitness from allele fitness values
def genotype_fitness(input_array, dominance): # Input_array will be allele_fitness
    
    IIIM_fem = input_array[0]
    YM_fem = input_array[1]
    traD_fem = input_array[2]
    IIIM_male = input_array[3]
    YM_male = input_array[4]
    
    IIIM_fem_h = 0
    YM_fem_h = 0
    IIIM_male_h = 0
    YM_male_h = 0
    
    # Dominance used here should be the dominance selected by the user
    # Treating dominance for proto Y chromosomes (IIIM and YM)!
    # Treating for Additive dominance scenario
    # h is for the "h" value
    if dominance == 'additive':
        IIIM_fem_h = 0.5
        YM_fem_h = 0.5
        IIIM_male_h = 0.5
        YM_male_h = 0.5
    
    if dominance == 'dominant':
        if IIIM_fem > 0:
            IIIM_fem_h = 0
        else: 
            IIIM_fem_h = 1

        if YM_fem > 0:
            YM_fem_h = 0
        else: 
            YM_fem_h = 1
        
        if IIIM_male > 0:
            IIIM_male_h = 0
        else: 
            IIIM_male_h = 1
        
        if YM_male > 0:
            YM_male_h = 0
        else: 
            YM_male_h = 1

    if dominance == 'recessive':
        if IIIM_fem > 0:
            IIIM_fem_h = 1
        else: 
            IIIM_fem_h = 0

        if YM_fem > 0:
            YM_fem_h = 1
        else: 
            YM_fem_h = 0
        
        if IIIM_male > 0:
            IIIM_male_h = 1
        else: 
            IIIM_male_h = 0
        
        if YM_male > 0:
            YM_male_h = 1
        else: 
            YM_male_h = 0

    if dominance == 'random':
        np.random.seed()
        rand_dominance = np.random.uniform(size = 4, low = 0, high = 1)  # array of size 4
        IIIM_fem_h = rand_dominance[0]
        YM_fem_h = rand_dominance[1]
        IIIM_male_h = rand_dominance[2]
        YM_male_h = rand_dominance[3]

    #Calculate single locus fitness effects of genotypes

    #IIIM Female
    if (IIIM_fem < 0):
        III_fem = 1
        IIIh_fem = 1 + (IIIM_fem_h * IIIM_fem)
        IIIm_fem = 1 + IIIM_fem # New value should be between 1 and 0.9
    else:
        III_fem = 1 - IIIM_fem
        IIIh_fem = 1 - (IIIM_fem_h * IIIM_fem)
        IIIm_fem = 1

    #YM Female
    if (YM_fem < 0):
        XX_fem = 1
        XY_fem = 1 + (YM_fem_h * YM_fem)
        YY_fem = 1 + YM_fem
    else:
        XX_fem = 1 - YM_fem
        XY_fem = 1 - (YM_fem_h * YM_fem)  #0.5 is "h"
        YY_fem = 1

    #traD Female    #would'nt have to touch tra at all, that's why we have 4 variables 
    if (traD_fem < 0):
        tra_fem = 1
        trad_fem = 1 + traD_fem
    else:
        tra_fem = 1 - traD_fem
        trad_fem = 1
    
    #IIIM Male
    if (IIIM_male < 0):
        III_male = 1
        IIIh_male = 1 + (IIIM_male_h * IIIM_male)
        IIIm_male = 1 + IIIM_male
    else:
        III_male = 1 - IIIM_male
        IIIh_male = 1 - (IIIM_male_h * IIIM_male)
        IIIm_male = 1

    #YM Male
    if (YM_male < 0):
        XX_male = 1
        XY_male = 1 + (YM_male_h * YM_male)   #h is only affecting the heterozygote
        YY_male = 1 + YM_male
    else: 
        XX_male = 1 - YM_male
        XY_male = 1 - (YM_male_h * YM_male)
        YY_male = 1

    # List of length 18
    # This array contains all the fitness permutations of the alleles possible in males and females (all 18 genotypes)
    # Note what happens if you remove trailing comma
    genotype_fit = [III_fem * tra_fem * XX_fem,     #   1 female f1 M0 F0 Y0    III/III     tra+/tra+   X/X [0]
                    III_fem * trad_fem * XX_fem,    #   2 female f2 M0 F1 Y0    III/III     traD/tra+   X/X [1]
                    IIIh_fem * trad_fem * XX_fem,   #   3 female f3 M1 F1 Y0    IIIM/III    traD/tra+   X/X [2]
                    IIIm_fem * trad_fem * XX_fem,   #   4 female f4 M2 F0 Y0    IIIM/IIIM   traD/tra+   X/X [3]
                    III_fem * trad_fem * XY_fem,    #   5 female f5 M0 F1 Y1    III/III     traD/tra+   X/Y [4]
                    IIIh_fem * trad_fem * XY_fem,   #   6 female f6 M1 F1 Y1    IIIM/III    traD/tra+   X/Y [5]
                    IIIm_fem * trad_fem * XY_fem,   #   7 female f7 M2 F1 Y1    IIIM/IIIM   traD/tra+   X/Y [6]
                    III_fem * trad_fem * YY_fem,    #   8 female f8 M0 F1 Y2    III/III     traD/tra+   Y/Y [7]
                    IIIh_fem * trad_fem * YY_fem,   #   9 female f9 M1 F1 Y2    IIIM/III    traD/tra+   Y/Y [8]
                    IIIm_fem * trad_fem * YY_fem,   #   10 female f10 M2 F1 Y2  IIIM/IIIM   traD/tra+   Y/Y [9]
                    IIIh_male * XX_male,            #   11 male m1 M1 F0 Y0     IIIM/III    tra+/tra+   X/X [10]
                    IIIm_male * XX_male,            #   12 male m2 M2 F0 Y0     IIIM/IIIM   tra+/tra+   X/X [11]
                    III_male * XY_male,             #   13 male m3 M0 F0 Y1     III/III     tra+/tra+   Y/X [12]
                    IIIh_male * XY_male,            #   14 male m4 M1 F0 Y1     IIIM/III    tra+/tra+   Y/X [13]
                    IIIm_male * XY_male,            #   15 male m5 M2 F0 Y1     IIIM/IIIM   tra+/tra+   Y/X [14]
                    III_male * YY_male,             #   16 male m6 M0 F0 Y2     III/III     tra+/tra+   Y/Y [15]
                    IIIh_male * YY_male,            #   17 male m7 M1 F0 Y2     IIIM/III    tra+/tra+   Y/Y [16]
                    IIIm_male * YY_male]            #   18 male m8 M2 F0 Y2     IIIM/IIIM   tra+/tra+   Y/Y [17]
                    
    genotype_norm = genotype_fit
    
    #Now to normalize genotype_norm

    max_male = max(genotype_norm[10:18])     #start at 10th position but do not include 18th position (all the males)
    max_female = max(genotype_norm[0:10])    #start at 0th position but do not include 10th position (all the females)

    #genotype_norm also returns 18 values
    for i in range(0,18):
        if i <= 9:
            genotype_norm[i] = genotype_fit[i] / max_female
            
        else:
            genotype_norm[i] = genotype_fit[i] / max_male


    dom_array = [IIIM_fem_h, YM_fem_h, IIIM_male_h, YM_male_h]
            
    return(genotype_norm, dom_array)
    

def sex_ratio(freqs):   # Function calculates sex ratio given an array of all possible genotypes
    sum_freq = sum(freqs)
    female = (freqs[0] + freqs[1] + freqs[2] + freqs[3] + freqs[4] + freqs[5] + freqs[6] + freqs[7] + freqs[8] + freqs[9]) / sum_freq
    #male = (freqs[10] + freqs[11] + freqs[12] + freqs[13] + freqs[14] + freqs[15] + freqs[16] + freqs[17])
    #male + female should = 1, so no need for the line above as the line below is faster
    male = 1 - female
    return (female, male)


def allele_freqs(freqs): 
    sum_freq = sum(freqs)

    # X and Y freqs fixed
    X_freq = (freqs[0] + freqs[1] + freqs[2] + freqs[3] + freqs[10] + freqs[11] 
    + 0.5*(freqs[4] + freqs[5] + freqs[6] + freqs[12] + freqs[13] + freqs[14])) / sum_freq
    
    #Y_freq1 = ((freqs[7] + freqs[8] + freqs[9] + freqs[15] + freqs[16] + freqs[17])
    #+ 0.5*(freqs[4] + freqs[5] +freqs[6] + freqs[12] + freqs[13] + freqs[14])) / sum_freq

    Y_freq = 1 - X_freq
   

    # III_freq and IIIM_freq
    III_freq = (freqs[0] + freqs[1] + freqs[4] + freqs[7] + freqs[12] + freqs[15]
	+ 0.5*(freqs[2] + freqs[5] + freqs[8] + freqs[10] + freqs[13] + freqs[16])) / sum_freq
    IIIM_freq = 1 - III_freq

    #tra_freq and traD_freq
    tra_freq = (freqs[0] + freqs[10] + freqs[11] + freqs[12] + freqs[13] + freqs[14] + freqs[15] + freqs[16] + freqs[17]
	+ 0.5*(freqs[1] + freqs[2] + freqs[3] + freqs[4] + freqs[5] + freqs[6] + freqs[7] + freqs[8] + freqs[9])) / sum_freq
    traD_freq = 1 - tra_freq

    return(X_freq, Y_freq, III_freq, IIIM_freq, tra_freq, traD_freq)


# This function will search to find where each genotype reaches its point of equilibrium
# Our high is the number of generations
def binary_search(arr, high):
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


def simu_v2_02():
    simulation_results = []                         # initialize the simulation_results array to having nothing in it before actually using the array
    simulation_results.append(initial_freq)         # Append initial_freq as the first index to simulation_results
             
    np.random.seed()                                # providing no parameter makes the random function completely unpredictable, adding a seed will make it predictable
    #Program now checks if user specifies range to be used for allele_fitness, checks if its usable, and default to [-1, 1] if not.
    
    allele_fitness = np.random.uniform(size = 5,low = allele_range*(-1),high = allele_range)
    
    # You will have a dominance array of size 4 with a low of 0 and a high of 1 (tra-D only exists in heterozygote, so won't be used here)
    # S tells you the difference between the homozygotes and h is the relative position of the heterozygote between the homozygotes
    
    fitness_array, dom_result = genotype_fitness(allele_fitness, dominance)         # Returns fitness array of 18 genotypes, only 1 fitness_array per simulation
    
    #This is where the recursion is happening
    for a in range(generations):                                
    
        result_recursion = recursion_seln(simulation_results[a-1], fitness_array)   # Takes genotoype frequencies and uses same genotype fitness array
        simulation_results.append(result_recursion)                                 # Add each simulation result to the next index


    geno_sim_results = np.array(simulation_results) # convert simulation results into a numpy array for slicing
    geno_asymp = []
    if switches[6] == '1':                              
        for genotypes in range(len(simulation_results[0])):         # for 1 through 18
            geno_asymp1 = binary_search(geno_sim_results[0:, genotypes], generations)   # binary search takes each genotype's data and the number of generations that were performed
            geno_asymp.append(geno_asymp1)
        
    # geno_asymp will contain all 18 points of equilibrium (the generation where equilibrium occurs)
    
    final_sim_result = simulation_results[generations]                          
    allele_freq = allele_freqs(final_sim_result)
    sex_ratios = sex_ratio(final_sim_result)
    
    #TEMPORARY-----
    #Leave uncommented when doing larger number of simulations to remove rows
    #with uninteresting allele frequencies
    # if ((allele_freq[0]<=.01 or allele_freq[0]>=.99) or
    #     (allele_freq[1]<=.01 or allele_freq[1]>=.99) or
    #     (allele_freq[2]<=.01 or allele_freq[2]>=.99) or
    #     (allele_freq[3]<=.01 or allele_freq[3]>=.99) or
    #     (allele_freq[4]<=.76 or allele_freq[4]>=.99) or
    #     (allele_freq[5]<=.01 or allele_freq[5]>=.24)):
    #     return None
    #--------------
    
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
   
    
def main_v2_02():
    
    p = Pool(processes = def_processes)                         # Set equal to 48 for carya cluster
    data = list(p.starmap(simu_v2_02, [() for _ in range (simulations)]))
    p.close()
    p.join()
    
    #TEMPORARY-----
    #Leave uncommented when doing larger number of simulations to remove rows
    #with uninteresting allele frequencies
    data[:] = [x for x in data if x != None]
    #--------------
    
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
    
    
def Debug(out):
    #Section dedicated to trouble shooting program by opening text file to be used as an output for temporary test statements
    #CAUTION!! IF called inside simulation, will be repeated every time
    if type(out).__name__ == "str":
        out = out + "\n"
        with open(os.path.join(cur_path, "Debugfile.txt"), "a") as f:
            f.write(out)


# The below code is needed to run main
# PURPOSE: Stops code from running when file is imported as a module. 
if __name__ == "__main__":
    
       
    #==========SECTION A Start==========#
    # This is the code that takes in arg parameters, simulations first, generations second
    parser = argparse.ArgumentParser(description = "Enter the number of simulations first and the number of generations second")
    parser.add_argument('SJobID', type=int)
    parser.add_argument('simulations', type = int)      # This is a positional argument, must be placed in this sequence in CLI
    parser.add_argument('generations', type = int)      # This is a positional argument, must be placed in this sequence in CLI
    parser.add_argument('def_processes', type = int)
    parser.add_argument('dominance', type = str)      # The '--' makes this argument optional, will default if not provided, must place '--dominance type' in CLI
    # Line below commented out to be replaced by line after to make control file optional
    # parser.add_argument('control_file', type = argparse.FileType('r', encoding='UTF-8'))
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
    
    #----------- File Creation -----------
    # Trust the user inputs correct path here instead
    now = datetime.now()
    cur_path = os.path.abspath(os.getcwd())             # This finds the current directory where the code is being run
    #if JobID != 0:       Commented out as JobID was made into a required parameter
    #    JobID = str(JobID)
    #else:
    #    JobID = now.strftime("%Y%m%d_%H%M%S")
    output = "data_" + str(JobID) + ".txt"
    modif_path = os.path.join(cur_path, "OutputArchive")
    
    #Check to see if the folder "OutputArchive" exist and create it if not
    try:
        os.mkdir(modif_path)
    except:
        pass
    
    #==========SECTION A END==========#


    #==========SECTION B Start==========#
    if control_file == None:                # Checks first to see if file is even an argument, if not, default to initial_freq
        initial_freq = [0.05556, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 
        0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556]
    else:                                   # If it is argument, do the following
        try:                    
            control_file = open(control_file, 'r')      # Open the file, if not found in directory, again default to initial_freq
        except FileNotFoundError:                       # Catches the exception if FileNotFound
            initial_freq = [0.05556, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 
            0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556, 0.05555, 0.05556]
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
    #==========SECTION B END==========#


    #==========SECTION C START==========#
    # Recreation of Data Frames to provide more detailed output.

    #dataout
    allele_fitness_dict = {'IIIM_Fem_fitness': [], 'YM_fem_fitness': [], 'traD_fem_fitness': [], 'IIIM_male_fitness': [], 'YM_male_fitness': []}
    fitness_array_dict = {'f1_fit': [], 'f2_fit': [], 'f3_fit': [], 'f4_fit': [], 'f5_fit': [], 'f6_fit': [], 'f7_fit': [], 'f8_fit': [], 'f9_fit': [], 'f10_fit': [],
                        'm1_fit': [], 'm2_fit': [], 'm3_fit': [], 'm4_fit': [], 'm5_fit': [], 'm6_fit': [], 'm7_fit': [], 'm8_fit': []}
    dom_result_dict = {'IIIM_fem_h': [], 'YM_fem_h': [], 'IIIM_male_h': [], 'YM_male_h': []}

    #dataout2
    final_sim_result_dict = {'f1': [], 'f2': [], 'f3': [], 'f4': [], 'f5': [], 'f6': [], 'f7': [], 'f8': [], 'f9': [], 'f10': [],
                        'm1': [], 'm2': [], 'm3': [], 'm4': [], 'm5': [], 'm6': [], 'm7': [], 'm8': []}
    allele_freq_dict = {'X_freq': [], 'Y_freq': [], 'III_freq': [], 'IIIM_freq': [], 'tra_freq': [], 'traD_freq': []}
    sex_ratios_dict = {'female' : [], 'male': []}

    genotype_equil_dict = {'f1_equil': [], 'f2_equil': [], 'f3_equil': [], 'f4_equil': [], 'f5_equil': [], 'f6_equil': [], 'f7_equil': [], 'f8_equil': [], 'f9_equil': [], 'f10_equil': [],
                        'm1_equil': [], 'm2_equil': [], 'm3_equil': [], 'm4_equil': [], 'm5_equil': [], 'm6_equil': [], 'm7_equil': [], 'm8_equil': []}

    df_allele_fitness = pd.DataFrame(allele_fitness_dict)
    df_fitness_array = pd.DataFrame(fitness_array_dict)
    df_dom_result = pd.DataFrame(dom_result_dict)
    df_final_sim_result = pd.DataFrame(final_sim_result_dict)
    df_allele_freq = pd.DataFrame(allele_freq_dict)
    df_sex_ratios = pd.DataFrame(sex_ratios_dict)
    df_genotype_equil = pd.DataFrame(genotype_equil_dict)
    #==========SECTION C END==========#
    
    
    #==========SECTION D START==========#
    # Run only one "section" at a time
    # Change in varaible names/uses:
    #    Put "simulations" in "total_sims"
    #    Simulation will never be more than the amount specified in "section"
    # Main function will run multiple times, doing one section until the required number of simulations is met
    # "section" can be changed at will to avoid creating a dataframe that is too large to fit in memory before
    # writing to file.
    section = 25000
    start = time.time()
    sim_output = simulations
    total_sims = simulations
    simulations = section
    while total_sims > section:
        main_v2_02()
        total_sims -= section
    simulations = total_sims
    main_v2_02()
    end = time.time()
    total_time = "{:.2f}".format(end-start)
    print(end - start, "seconds") 
    #==========SECTION D END==========#
    
    
    #==========SECTION E START==========#
    # Places all the parameters used into a single file
    
    with open(os.path.join(modif_path, "ParameterLog.txt"), 'a') as f:
        f.write(str(JobID))
        f.write(" Precessors:")
        f.write(str(def_processes))
        f.write(" Time:")
        f.write(total_time)
        f.write("s Simulations:")
        f.write(str(sim_output))
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
    #==========SECTION E END==========#