#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 23:28:11 2026

@author: shawt5
"""
import os, sys
import argparse, datetime
import numpy as np
import pandas as pd
import glob

HiC_A_COL='#chr1'
HiC_B_COL='chr2'
HiC_A_START_COL='x1'
HiC_A_END_COL='x2'
HiC_B_START_COL='y1'
HiC_B_END_COL='y2'
INDEX_COL='INDEX'


def StrinToList(label):
    STRING = label.replace('"', '').replace("'", '').replace("[", '').replace("]", '')
    CELL=[]
    for cell in  STRING.split(','):
        CELL.append(cell.strip())
    return CELL

def combineATACpro(inputFiles, OutputDir):   
    inputFiles= inputFiles.split(' ')
    PATTERN=[]
    for Filename in inputFiles:
        patternFile = Filename.split('_chr')[0]
        PATTERN.append(patternFile)
    PATTERN = list(set(PATTERN).intersection(set(PATTERN)))
    for pattern in PATTERN:
        DATA=[]
        for Filename in inputFiles:
            if Filename.find(pattern)!=-1:
                filepathA = os.path.join(OutputDir, Filename)
                filepathB = Filename
                filepath = filepathA if os.path.exists(filepathA) else filepathB  
                data = pd.read_csv(filepath, sep='\t')
                DATA.append(data)
        DATA = pd.concat(DATA, axis = 0)
        DATA = DATA.drop_duplicates()
        DATA.to_csv(pattern+'.tsv', sep='\t', index=False)

def reIndex(data):
    CHROM_COL = data.columns[data.columns.str.find(HiC_A_COL)!=-1][0]
    X1_COL = data.columns[data.columns.str.find(HiC_A_START_COL)!=-1][0]
    X2_COL = data.columns[data.columns.str.find(HiC_A_END_COL)!=-1][0]
    Y1_COL = data.columns[data.columns.str.find(HiC_B_START_COL)!=-1][0]
    Y2_COL = data.columns[data.columns.str.find(HiC_B_END_COL)!=-1][0]
    INDEX_COL_All = data.columns[data.columns.str.find(INDEX_COL)!=-1]
    data = data.set_index(data[CHROM_COL].astype(str)+':'+data[X1_COL].astype(int).astype(str)+':'+
                          data[X2_COL].astype(int).astype(str)+':'+data[Y1_COL].astype(int).astype(str)+':'+data[Y2_COL].astype(int).astype(str))
    data = data.sort_values(by=[CHROM_COL, X1_COL, X2_COL, Y1_COL, Y2_COL])
    data['CHECK'] = data.index
    i = 1
    for loop_index in data['CHECK'].unique():
        for COL in INDEX_COL_All:
            data.loc[loop_index, COL] = i
        i+=1
    data = data.drop(columns=['CHECK'])
    return data

def combineToBacon(inputFolder, OutputDir, SamplePair):
        dsamplepair = pd.read_csv(SamplePair, sep='\t')
        Target = dsamplepair.loc[0, 'Target']
        Reference = dsamplepair.loc[0, 'Reference']
        File_endstring='ToBacon_'+Target+'_vs_'+Reference.replace(' ', '').replace(',','_')
        ##inputFiles = StrinToList(inputFiles)
        inputFiles= glob.glob(inputFolder+'ToBacon_*')
        DATA=[]
        for Filename in inputFiles:
            Filename = Filename.split('/')[-1]
            if Filename.find(File_endstring)!=-1:
                filepathA = os.path.join(OutputDir, Filename)
                filepathB = Filename
                filepath = filepathA if os.path.exists(filepathA) else filepathB  
                data = pd.read_csv(filepath, sep='\t')
                DATA.append(data)
        DATA = pd.concat(DATA, axis = 0)
        DATA = DATA.drop_duplicates()
        DATA = reIndex(DATA)
        DATA.to_csv(File_endstring+'.tsv', sep='\t', index=False)


def RUN(inputFolder, OutputDir, SamplePair, Pattern):
    if inputFolder==None:
        print ('No files were found')
    else:
        if Pattern==None:
            combineATACpro(inputFolder, OutputDir)
        elif Pattern=='bacon':
            combineToBacon(inputFolder, OutputDir, SamplePair)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFolder", type=str,default=None,
                        help="main directory")    
    parser.add_argument("-o", "--OutputDir", type=str,default='./output/',
                        help="default output folder")
    parser.add_argument("-sp", "--SamplePair", type=str,default=None,
                        help="import Sample Combinations")
    parser.add_argument("-p", "--Pattern", type=str,default=None,
                        help="import Sample Combinations")
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.inputFolder, args.OutputDir, args.SamplePair, args.Pattern)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
