import os, sys, subprocess
from shutil import rmtree
from stat import *
import pandas as pd

def baseNameFormat(df, cols, delim = "-"):
    # Input: dataframe, cols = colnames in order of basename
    # delim = string delimiter for basename column values
   
    # Checks:
    #
    # -- cols are members of df:
    nonExistColnames = [c for c in cols if c not in list(df)]
    if len(nonExistColnames) != 0:
        errorString = "Some columns in 'cols' don't exist:\n" + " ".join(nonExistColnames)
        raise ValueError(errorString)

    ###
    # make format string:
    baseNameFormat = delim.join(["{}" for c in cols])
    return(baseNameFormat)

def baseNameFromCols(df, cols, formatString, outputColName = "baseName"):
    # Appends basename column to dataframe
    # Input:
    df[outputColName] = df[cols].apply(lambda x : formatString.format(*x), axis = 1)
    return(df)

def addBaseName(df, cols, delim = "-", outputColName = "baseName"):
    # Takes sampleSheet & list of columns as input, 
    # returns sampleSheet with 'basename' column appended
    formatString = baseNameFormat(df, cols, delim)
    df = baseNameFromCols(df, cols, formatString, outputColName = outputColName)
    return(df)

#def addExt(df, ext, baseNameCol = 'basename'):
#    # input sampleSheet, return extension of file

def makeSampleSheets(path, idcols, delim, baseNameColumn = "baseName", fileDelimiter = "\t"):
    df = pd.read_table(path, delimiter = fileDelimiter)
    df = addBaseName(df, idcols, delim)
    keep_cols = idcols.copy()
    keep_cols.append(baseNameColumn)
    pool_df = df[keep_cols].drop_duplicates()
    return(df, pool_df)

def move_fastq(read1, read2, baseNames):
    # rename fastq read1 and read2 files to basename_R1.fq.gz basename_R2.fq.gz
    # vectorized rename
    # requires .fastq.gz filetype
    # TODO: make generic & use input type as output type
    ## TODO:
    # Make better solution to this:
    # -- add column to sampleSheet with old fastq names so they can be un-renamed,
    # -- change fastq_r1 and fastq_r2 to NEW filename in sampleSheet
    fastq_dir = 'Fastq/'
    fastqs = [fastq_dir + f for f in os.listdir(fastq_dir)]
    fq_permissions = [os.stat(f) for f in fastqs]
    
    # Add write permission to all bits so script can remove the tree
    [os.chmod(f, S_IREAD + S_IWRITE + S_IWGRP + S_IWOTH) for f in fastqs]
    
    # Remove fastq_dir recursively
    if os.path.exists(fastq_dir):
        #os.rmdir(fastq_dir)
        rmtree(fastq_dir)

    os.mkdir(fastq_dir)

    for r1, r2, baseName in zip(read1, read2, baseNames):
        #shutil.copyfile(r1, fastq_dir + baseName + "_R1.fastq.gz")
        #shutil.copyfile(r2, fastq_dir + baseName + "_R2.fastq.gz")
        copy_r1 = lambda fq : ["cp", fq, fastq_dir + baseName + "_R1.fastq.gz"]
        copy_r2 = lambda fq : ["cp", fq, fastq_dir + baseName + "_R2.fastq.gz"]
        subprocess.run(copy_r1(r1))
        subprocess.run(copy_r2(r2))


def main(path, idcols, delim = '-'):
    # format & write sampleSheet/poolSampleSheet
    sampleSheet, poolSampleSheet = makeSampleSheets(path, idcols, delim)
    
    out = {"sampleSheet": sampleSheet, "poolSampleSheet": poolSampleSheet}
    
    for df, name in zip(out.values(), out.keys()):
        filename = name + ".tsv"
        df.to_csv(filename, sep = "\t", index = False)
    # Copy & rename fastqs
    move_fastq(sampleSheet.fastq_r1, sampleSheet.fastq_r2, sampleSheet.baseName)
    

if __name__ == "__main__":
    # Drop call from argv
    sys.argv.pop(0)
    # Path to configFile is first argument
    path = sys.argv.pop(0)
    # TODO: check that path exists ?
    # Remaining function calls are the id variables for the basename
    idcols = sys.argv

    main(path, idcols, '-')
