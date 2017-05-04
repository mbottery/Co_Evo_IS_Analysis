#!/usr/local/bin/python3.5

# Script used to identity possible insertion sequences using read template length. The 
# script takes aligned reads from an indexed bed file and calculates reads with a 
# template length within the top 99th percentile for that reference. Read frequencies are 
# binned to create a histogram, with peaks representing areas within the genome that have 
# an over representation of reads with long template length. These were subsequently 
# confirmed or rejected as IS sequences through visualisation using Integrative Genomics 
# Viewer (http://software.broadinstitute.org/software/igv/)

# import modules
from optparse import OptionParser
from optparse import OptionGroup
import subprocess
import pysam
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    parser = OptionParser()

    # get command line options
    REQUIRED = OptionGroup(parser,"REQUIRED")
    REQUIRED.add_option("-b",action="store",type="string",dest="IN",help="Input bam file")
    OPTIONAL = OptionGroup(parser,"OPTIONAL")
    OPTIONAL.add_option("-o",action="store",type="string",dest="OUT",default="",help="Out file name")
    OPTIONAL.add_option("-w",action="store",type="int",dest="BIN",default="5000",help="Number of bins per template")
    OPTIONAL.add_option("-p",action="store",choices=['y','n'],dest="PLOT",default="y",help="Produce Plot, y/n")
    parser.add_option_group(REQUIRED)
    parser.add_option_group(OPTIONAL)
    (options, args) = parser.parse_args()

    # Verify arguments
    if options.IN == None:
        parser.error("Must provide input bam file with -b")
    return(options)

class OptionsInsertSize():
    def __init__(self,b,o=""):
        self.IN=b
        self.OUT=o
        self.BIN=w
        self.PLOT=p
    
def plotLongReads(longReads,bin,sample,ref_name,ref_length,plot):
    (n, bins, patches) = plt.hist(longReads[:,1],bins = bin)
    plt.xlabel("Reference Position")
    plt.ylabel("Long read frequency")
    axis = plt.gca()
    axis.set_xlim([0,ref_length])
    plt.title("{0} - {1}".format(sample,ref_name))
    if plot == 'y':
        plt.show()
    return n, bins

def InsertSize(options):
    # Open bam file
    bamfile = pysam.AlignmentFile(options.IN, "rb")
    # Open out file
    outfn = "{0}_large_insert.txt".format(options.IN.split('/')[-1].split('.rmdup')[0])
    out = open(outfn,'w')
    
    # required indexed bam file
    assert bamfile.header['HD']['SO']=='coordinate', 'bam file must be indexed'
    print bamfile.header['SQ']
        
    # Find reference sequence names from reference sequence dictionary
    for SQ in bamfile.header['SQ']:
        SN=(SQ['SN'])
        print SN
        reads=[]
        # Get observed template lengths, and read start and end
        for read in bamfile.fetch(SN):
            if read.is_paired:
                reads.append([abs(read.template_length),read.reference_start,read.reference_end])
        reads=np.array(reads)
        # Find 99th percentile of observed template lengths
        percentile99 = np.percentile(reads,99,axis=0)[0]
        # Filter reads based on top 99th % template length
        reads = reads[reads[:,0]>percentile99]
        # reads[:,0] = read template length
        # reads[:,1] = read start
        # reads[:,2] = read end
        print "{0} reads with observed template lenth in top 99th percentile".format(len(reads))
        
        # plot density
        (n, bins) = plotLongReads(reads,options.BIN,options.IN,SN,SQ['LN'],options.PLOT)
        
        # write hist data to out file
        out.write('Ref\tFrequency\tBin_Start\n')
        for i in range(len(n)):
            out.write('{0}\t{1}\t{2}\n'.format(SN,str(n[i]),str(bins[i])))
    out.close()
    
if __name__=="__main__":
    InsertSize(parse_args())