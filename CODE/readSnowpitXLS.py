# readSnowpitXLS
# hpm 18 May 2019

import pandas as pd
import glob
import numpy as np
import os

def nanmean(d): # function to get mean of density data
    col=list(d) # get a list of columns
    pos=0
    q=[]
    for icol in col:
        d2=pd.to_numeric(d[icol],errors='coerce')
        q.append(d2.mean(skipna=True)) # mean of profile
        pos += 1
    mrho=np.mean(q) # mean of profile means (IMPROVE!)
    return mrho

def readSnowpit(filename):
    print(filename)
    fname=os.path.basename(filename) # get just the filename
    pitID=fname[0:3] # get just the pitID
    xl = pd.ExcelFile(filename) # make an excel object
    d=pd.read_excel(xl,sheet_name=0,usecols='L') # read a column
    pitE=d['Surveyors:'][2] # get easting
    d=pd.read_excel(xl,sheet_name=0,usecols='P') # read a column
    pitN=d['Unnamed: 15'][2] # get northing
    d=pd.read_excel(xl,sheet_name=0,usecols='G') # read a column
    pitHS=d['Unnamed: 6'][4] # get snowdepth
    d=pd.read_excel(xl,sheet_name=0,usecols='E:G') # read a column
    d2=d.iloc[9:,:] # get density data
    mrho=nanmean(d2) # get mean density
    pitMeanRho=float("{0:.1f}".format(mrho)) # format to have one decimal point
    pitSWE=pitMeanRho*pitHS/100 # calculate SWE
    pitSWE=float("{0:.1f}".format(pitSWE)) # format for 1 dec point
    column=['pitID','Easting','Northing','SnowDepth[cm]','MeanDensity[kgpm3]','SWE[mm]'] # define columns of dataframe
    data=[[pitID,pitE,pitN,pitHS,pitMeanRho,pitSWE]] # make list of data
    df = pd.DataFrame(data,columns=column) # create dataframe
    return df # output result

pitfiles=glob.glob('../PITS/*.xlsx') # get all pit files
frames = [ readSnowpit(f) for f in pitfiles ] # iterate over pits with list comprehension
result = pd.concat(frames, ignore_index=True) # concatenate the results
print(result) # print results
# export results to csv file
result.to_csv('../PitSummary.csv', index=None, header=True, na_rep='NaN') # print to csv
