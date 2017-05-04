#InsertSize.py
Script used to identity possible insertion sequences using read template length. The script takes aligned reads from an indexed bed file and calculates reads with a template length within the top 99th percentile for that reference. Read frequencies are binned to create a histogram, with peaks representing areas within the genome that have an over representation of reads with long template length. These were subsequently confirmed or rejected as IS sequences through visualisation using Integrative Genomics Viewer (http://software.broadinstitute.org/software/igv/).

Python 2.7.10

Requires:
numpy 1.11.1 
matplotlib 1.5.1
pysam 0.8.4
