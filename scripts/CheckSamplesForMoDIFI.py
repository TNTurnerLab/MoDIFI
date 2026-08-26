#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:33:32 2025

@author: shawt5
"""

import numpy as np
import pandas as pd
import os, sys
import argparse, datetime
import warnings
import glob


LABEL_HiC='HiC'
LABEL_PRO='Promoter'
LABEL_COL='Label'
FILENAME_COL='Filename'
TARGET_COL='Target'
REF_COL='Reference'
## DEseq
CHROM_COL='CHROM'
START_COL='Start'
END_COL='End'
COORDS_COL='Coords'
Log2FC_COL='log2FoldChange'
FC_COL = 'FoldChange'
FCnorm_COL='NormFoldChange'
GENEID_COL = 'GeneID'
GENENAME_COL='GeneName'
Stat_Score='stat'
Z_Score='Z'
LFCSE_COL='lfcSE'
PValue_COL='pvalue'
Padj_COL='padj'



HiC_A_COL='#chr1'
HiC_B_COL='chr2'
HiC_A_START_COL='x1'
HiC_A_END_COL='x2'
HiC_B_START_COL='y1'
HiC_B_END_COL='y2'
CHECK_COL='Check'
OL_COL='Overlap'
LABEL_PASS='PASS'
LABEL_XY='XY'
LABEL_X='X'
LABEL_Y='Y'
GENEID_X_COL='X_GeneID'
GENEID_Y_COL='Y_GeneID'
INDEX_COL='INDEX'
STRAND_COL='Strand'


def checkColumnHic(f, OutputDir, cell, Nmin=1000):
    files = glob.glob(f)
    HIC=[]
    for file in files:
        dd_hic = pd.read_csv(file, sep='\t')
        if dd_hic.columns[0]!=HiC_A_COL:
            dd_hic = pd.read_csv(file, sep='\t', header=None)
            dd_hic = dd_hic[dd_hic.columns[:6]]
            for i, col in enumerate([HiC_A_COL, HiC_A_START_COL, HiC_A_END_COL, HiC_B_COL, HiC_B_START_COL, HiC_B_END_COL]):
                dd_hic = dd_hic.rename(columns={i:col})
        HIC.append(dd_hic)
    HIC = pd.concat(HIC, axis = 0)
    HIC = HIC.drop_duplicates()
    HIC = HIC.set_index(np.arange(HIC.shape[0]))
    HICPATH=[]
    if HIC.shape[0]<=Nmin: 
        hicpath= cell+'_MergedLoop.tsv'
        HIC.to_csv(hicpath, sep='\t',index=False)
        HICPATH.append(hicpath)
    else:
        N = int(np.ceil(HIC.shape[0]/Nmin))
        for i in range(N):
            if i < (N-1):
                start = i * Nmin 
                end = ((i+1) * Nmin) -1
                partial = HIC.loc[start:end]
            else:
                start = i * Nmin 
                partial = HIC.loc[start:]
            hicpath= cell+'_MergedLoop'+str(i+1)+'.tsv'
            partial.to_csv(hicpath, sep='\t',index=False)
            HICPATH.append(hicpath)
    return HICPATH



def generateSampleInfoFromDEseqOutputs(PromoterPath, ResourcesDir, CELL, FileName, OutputDir, LABEL_RNA, LABEL_ATAC):
    ## RNAseq
    dSampleInfo = pd.DataFrame(columns=[LABEL_COL, FILENAME_COL, TARGET_COL, REF_COL])
    itern=0
    for LABEL in [LABEL_RNA, LABEL_ATAC]:
        files = glob.glob(OutputDir+'*'+LABEL+'.txt')
        if files:
            for f in files:
                dSampleInfo.loc[itern, LABEL_COL] = LABEL
                dSampleInfo.loc[itern, FILENAME_COL] = f
                dSampleInfo.loc[itern, TARGET_COL] = f.split('/')[-1].split('_')[0]
                dSampleInfo.loc[itern, REF_COL] = f.split('/')[-1].split('_')[2]
                itern+=1
    ## HiC
    for i, cell in enumerate(CELL):
        files = FileName[i]
        HICPATH = checkColumnHic(files, OutputDir, cell)
        for hicpath in HICPATH:
            dSampleInfo.loc[itern, LABEL_COL] = LABEL_HiC
            dSampleInfo.loc[itern, FILENAME_COL] = hicpath
            dSampleInfo.loc[itern, TARGET_COL] = cell
            itern+=1
    ## promoters
    files = PromoterPath
    if len(files)==1:
        dSampleInfo.loc[itern, LABEL_COL] = LABEL_PRO
        dSampleInfo.loc[itern, FILENAME_COL] = files
        itern+=1
    else:
        dSampleInfo.loc[itern, LABEL_COL] = LABEL_PRO
        dSampleInfo.loc[itern, FILENAME_COL] = ResourcesDir+'Promoter/Promoter.tsv'
    return dSampleInfo
        

def PathReformat(ResourcesDir):
    if ResourcesDir[-1]=='/':
        ResourcesDir = ResourcesDir
    else:
        ResourcesDir = ResourcesDir+'/'  
    return ResourcesDir
        
    

def importData(SampleInfo, SamplePair, PromoterPath, ResourcesDir, CELL, FileName, OutputDir, LABEL_RNA, LABEL_ATAC, RUN_WITH_WARNINGS=False):
    ## Import data and normalize RNA-seq amd ATAC-seq
    STOP=False
    dSampleInfo = generateSampleInfoFromDEseqOutputs(PromoterPath, ResourcesDir, CELL, FileName, OutputDir,  LABEL_RNA, LABEL_ATAC)
    if SamplePair!=None:
        if os.path.exists(SamplePair):
            dSamplePair = pd.read_csv(SamplePair, sep='\t')
        else:
            dSamplePair=False
    else:
        dSamplePair=False
    RNA={}
    ATAC={}
    HiC={}
    PRO={}
    for label in [LABEL_ATAC, LABEL_RNA, LABEL_HiC, LABEL_PRO]:
        sample = dSampleInfo[dSampleInfo[LABEL_COL]==label]
        if label==LABEL_ATAC:
            for index in sample.index:
                Filename = sample.loc[index, FILENAME_COL]
                Target = sample.loc[index, TARGET_COL]
                REF = sample.loc[index, REF_COL]
                ATAC[Target+'/'+REF]=1
                ## Reverse
        elif label==LABEL_RNA:
            for index in sample.index:
                Filename = sample.loc[index, FILENAME_COL]
                Target = sample.loc[index, TARGET_COL]
                REF = sample.loc[index, REF_COL]
                RNA[Target+'/'+REF] = 1
        elif label==LABEL_HiC:
            for index in sample.index:
                Filename = sample.loc[index, FILENAME_COL]
                Target = sample.loc[index, TARGET_COL]
                HiC[Target] = 1
        elif label==LABEL_PRO:
            for index in sample.index:
                Filename = sample.loc[index, FILENAME_COL]
                Target = sample.loc[index, TARGET_COL]
                if np.isnan(Target)==True:
                    PRO = pd.read_csv(Filename, sep='\t')
                else:
                    PRO[Target] = pd.read_csv(Filename, sep='\t')
    ## Check the imported data are enough to calculate MoDIFI
    dSamplePair, STOP = CheckImportedData(RNA, ATAC, HiC, PRO, dSampleInfo, dSamplePair, STOP)
    
    if SamplePair==None or SamplePair=='None':
        SamplePairOut= 'SamplePair.tsv'
        dSamplePair.to_csv(SamplePairOut, sep='\t', index=False)
    else:
        SamplePairOut = SamplePair.split('/')[-1]
        dSamplePair.to_csv(SamplePair, sep='\t', index=False)
    for i in range(dSamplePair.shape[0]):
        dSamplePair.loc[[i]].to_csv('piece_'+str(i)+'.tsv', sep='\t', index=False)
    if SampleInfo==None or SampleInfo=='None':
        #dSamplePair.to_csv(ResourcesDir+SamplePairOut, sep='\t', index=False)
        SampleInfoOut='SampleInfo.tsv'
        dSampleInfo.to_csv(SampleInfoOut, sep='\t', index=False)
    else:
        SampleInfoOut = SampleInfo.split('/')[-1]
        dSampleInfo.to_csv(SampleInfoOut, sep='\t', index=False)
    if STOP==True and RUN_WITH_WARNINGS==False:
        print ('!!!Warning!!! Script Stopped, missing some inputs!!')
        print ('Missing input data shown above or details in this file %s%s. Please fill the necessary data information in %s' %(ResourcesDir, SamplePairOut,SampleInfoOut))
        print ('Or if you would like the statistics generated for combinations with fully available data you can turn on RUN_WITH_WARNINGS (--RUN_WITH_WARNINGS)')
        print ("")
        sys.exit(1)
    return RNA, ATAC, HiC, PRO, dSamplePair


def CheckImportedData(RNA, ATAC, HiC, PRO, dSampleInfo, dSamplePair, STOP):
    ## Check whether the input data are sufficient to calculate MoDIFI
    if isinstance(dSamplePair, pd.DataFrame):
        for index in dSamplePair.index:
            Target = dSamplePair.loc[index, TARGET_COL]
            REF = dSamplePair.loc[index, REF_COL].split(',')
            cleaned_values = [value.strip() for value in REF]
            print ('Target:%s, Ref:%s' %(Target,  dSamplePair.loc[index, REF_COL]))
            lackString=''
            for ref in cleaned_values:
                name =Target+'/'+ref
                if name  not in list(RNA.keys()):
                    print ('Missing RNA-seq for %s/%s' %(Target, ref))
                    lackString = lackString+ 'Missing RNA-seq for %s/%s' %(Target, ref)+'\n'
                    STOP=True
                if name not in list(ATAC.keys()):
                    print ('Missing ATAC-seq for %s/%s' %(Target, ref))
                    lackString = lackString+ 'Missing ATAC-seq for %s/%s' %(Target, ref)+'\n'
                    STOP=True
            if Target not in list(HiC.keys()):
                print ('Missing HiC for %s' %(Target))
                STOP=True
                lackString = lackString+ 'Missing HiC for %s' %(Target)+'\n'
            if lackString!='':
                dSamplePair.loc[index, CHECK_COL]=lackString[:-1]
            else:
                dSamplePair.loc[index, CHECK_COL]=LABEL_PASS
                print ('PASS')
            print ('')
    else: ### if SamplePair=None (dSamplePair=False)
        CellT = dSampleInfo[dSampleInfo[TARGET_COL].notna()][TARGET_COL].unique()
        CellR = dSampleInfo[dSampleInfo[REF_COL].notna()][REF_COL].unique()
        CELLS = list(set(list(CellT)+list(CellR)))
        ddSamplePair = pd.DataFrame(columns=[TARGET_COL, REF_COL, CHECK_COL])
        itern = 0
        for cell in CELLS:
            Target = cell
            ddSamplePair.loc[itern, TARGET_COL] = cell
            cleaned_values = CELLS.copy()
            cleaned_values.remove(cell)
            ref_string=''
            for string in cleaned_values:
                ref_string = ref_string+ string + ', '
            ref_string=ref_string[:-2]
            ddSamplePair.loc[itern, REF_COL] = ref_string
            print ('Target:%s, Ref:%s' %(cell, ref_string))
            lackString=''
            for ref in cleaned_values:
                name =Target+'/'+ref
                if name  not in list(RNA.keys()):
                    print ('Missing RNA-seq for %s/%s' %(Target, ref))
                    STOP=True
                    lackString = lackString+ 'Missing RNA-seq for %s/%s' %(Target, ref)+'\n'
                if name not in list(ATAC.keys()):
                    print ('Missing ATAC-seq for %s/%s' %(Target, ref))
                    STOP=True
                    lackString = lackString+ 'Missing ATAC-seq for %s/%s' %(Target, ref)+'\n'
            if Target not in list(HiC.keys()):
                print ('Missing HiC for %s' %(Target))
                STOP=True
                lackString = lackString+ 'Missing HiC for %s' %(Target)+'\n'
            if lackString!='':
                ddSamplePair.loc[itern, CHECK_COL]=lackString[:-1]
            else:
                ddSamplePair.loc[itern, CHECK_COL]=LABEL_PASS
                print ('PASS')
            print ('')
            itern+=1
        dSamplePair = ddSamplePair.copy()
    return dSamplePair, STOP

def StrinToList(label):
    STRING = label.replace('"', '').replace("'", '').replace("[", '').replace("]", '')
    CELL=[]
    for cell in  STRING.split(','):
        CELL.append(cell.strip())
    return CELL

def RUN(SampleInfo, SamplePair, PromoterPath, ResourcesDir, label, HiCFileName, OutputDir, LABEL_RNA, LABEL_ATAC, RUN_WITH_WARNINGS=False):
    warnings.filterwarnings("ignore")
    ResourcesDir = PathReformat(ResourcesDir)
    OutputDir = PathReformat(OutputDir)
    CELL=StrinToList(label)
    FileName=StrinToList(HiCFileName)
    RNA, ATAC, HiC, PRO, dSamplePair = importData(SampleInfo, SamplePair, PromoterPath, ResourcesDir, CELL, FileName,  OutputDir,LABEL_RNA, LABEL_ATAC, RUN_WITH_WARNINGS)
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ResourcesDir", type=str,default='./resources',
                        help="main directory")
    parser.add_argument("-p", "--PromoterPath", type=str,default='./resources/Promoter/Promoter.tsv',
                        help="Promoter File Path")
    parser.add_argument("-si", "--SampleInfo", type=str,default=None,
                        help="import RNA-seq (DESeq formate), ATAC-seq, HiC and GeneList")
    parser.add_argument("-sp", "--SamplePair", type=str,default=None,
                        help="import Sample Combinations")
    parser.add_argument("-l", "--label", type=str, default='["HEK293", "IMR32"]',
                    help="List of labels for cell lines")
    parser.add_argument("-c", "--HiCFileName", type=str,default='[./resources/hek293/hic/*, ./resources/imr32/hic/*]',
                        help="path")
    parser.add_argument("-o", "--OutputDir", type=str,default='./output',
                        help="default output folder") 
    parser.add_argument("--LABEL_RNA", type=str,default='RNAseq',
                        help="labels for RNAseq from DEseq")    
    parser.add_argument("--LABEL_ATAC", type=str,default='ATACseq',
                        help="labels for RNAseq from DEseq")    
    parser.add_argument('--RUN_WITH_WARNINGS', action='store_true',
                        help="Run the script while allowing warnings")
    parser.set_defaults(RUN_WITH_WARNINGS=False)
    
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    RUN(args.SampleInfo, args.SamplePair, args.PromoterPath, args.ResourcesDir, args.label, args.HiCFileName, args.OutputDir, args.LABEL_RNA, args.LABEL_ATAC, args.RUN_WITH_WARNINGS)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()
