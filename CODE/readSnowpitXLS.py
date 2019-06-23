# readSnowpitXLS
# hpm 18 May 2019

import pandas as pd
import glob
import numpy as np
import os

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
    d=pd.read_excel(xl,sheet_name=0,usecols='E:G') # read 3 columns
    d2=d.iloc[9:,:] # get all density data
    d3=d2.values.flatten() # make into one 1-D array
    d4=pd.to_numeric(d3) # make strings into floats for computation
    mrho=np.nanmean(d4) # calculate the mean, ignoring NaN values
    pitMeanRho=float("{0:.1f}".format(mrho)) # format to have one decimal point
    pitSWE=pitMeanRho*pitHS/100 # calculate SWE
    pitSWE=float("{0:.1f}".format(pitSWE)) # format for 1 dec point
    column=['pitID','Easting','Northing','SnowDepth[cm]','MeanDensity[kgpm3]','SWE[mm]'] # define columns of dataframe
    data=[[pitID,pitE,pitN,pitHS,pitMeanRho,pitSWE]] # make list of data
    df = pd.DataFrame(data,columns=column) # create dataframe
    return df # output result

# !find . -name "~*" -type f -delete  # need to figure out how to do this within python. Issues with " and ~
pitfiles=glob.glob('../PITS/*.xlsx') # get all pit files
frames = [ readSnowpit(f) for f in pitfiles ] # iterate over pits with list comprehension
result = pd.concat(frames, ignore_index=True) # concatenate the results
print(result) # print results
# export results to csv file
result.to_csv('../PitSummary.csv', index=None, header=True, na_rep='NaN') # print to csv
