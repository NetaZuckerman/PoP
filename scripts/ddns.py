#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 12:13:32 2023

@author: hagar
"""


MIN_READ_LEN = 1000
MAX_READ_LEN = 1300

import pandas as pd 
from pathlib import Path
import csv
import subprocess
from utils.utils import create_dirs, mafft, change_header
import os
import pysam
from io import StringIO
from statistics import mean
import gzip
from scripts.mutations.signatures import run_ddns
from Bio import SeqIO




def get_barcode_name(barcode_csv):
    
    with open(barcode_csv, mode='r') as infile:
        reader = csv.reader(infile)
        barcodes = dict((rows[0],rows[1]) for rows in reader)
    
    return barcodes
    



def filter_reads_by_length(dir_in,reads_out):
    
    for barcode in os.listdir(dir_in):
        if not barcode == "unclassified" and not barcode == reads_out:
       
            if not os.path.exists(reads_out + barcode):
                os.mkdir(reads_out + barcode)
                
            for f in os.listdir(dir_in + barcode):
                fastq_records = []
                with open(reads_out + barcode + "/" + f,"w") as fw:
                    if f.endswith(".gz"):
                            with gzip.open(os.path.join(dir_in + barcode, f), "rt") as handle:
                                for record in SeqIO.parse(handle, "fastq"):
                                    length = len(record)
                                    if length > MIN_READ_LEN and length < MAX_READ_LEN:
                                        fastq_records.append(record)
                                        
        
                    elif f.endswith(".fastq") or f.endswith(".fq"):
                        for record in SeqIO.parse(os.path.join(dir_in + barcode,f),"fastq"):
                            length = len(record)
                            if length > MIN_READ_LEN and length < MAX_READ_LEN:
                                fastq_records.append(record)
        
                    SeqIO.write(fastq_records,fw, "fastq")






def minimap(fastq_path, reference, barcodes, output):
    minimap = "minimap2 -ax map-ont %(ref)s %(fastq_dir)s/*.fastq* | \
                samtools sort > %(output)s.bam"
    split = "bamtools split -in %(bam)s.bam -reference"
    
    for fastq_dir in os.listdir(fastq_path):
        if os.path.isdir(fastq_path+fastq_dir) and fastq_dir in barcodes.keys():
            sample = barcodes[fastq_dir]
            subprocess.call(minimap % dict(ref=reference, fastq_dir=fastq_path+fastq_dir,\
                                           output=output+"BAM/"+sample), shell=True)
            subprocess.call(split % dict(bam=output+"BAM/"+sample), shell=True)
            
            
def depth(output):
    samtools_index = "samtools index %(bam)s"
    depth = "samtools depth -a %(bam)s > %(output)s.txt"
    for bam in os.listdir(output+"BAM/"):
        if "bai" not in bam:
            subprocess.call(samtools_index % dict(bam=output+"BAM/"+bam), shell=True)
            if "REF" in bam:
                subprocess.call(depth % dict(bam=output+"BAM/"+bam, \
                                             output=output+"depth/"+bam.split(".bam")[0]), shell=True) 

def cns(output):
    for bam in os.listdir(output+"BAM/"):
        if "REF" in bam and not "bai" in bam and not "unmapped" in bam:
            sample = bam.split(".bam")[0]
            ivar_cns = "samtools mpileup -A %(bam)s.bam | ivar consensus -t 0.6 -m 1 -p %(cns)s.fa"
            ivar_cns5 = "samtools mpileup -A %(bam)s.bam | ivar consensus -t 0.6 -m 5 -p %(cns)s.fa"
            subprocess.call(ivar_cns % dict(bam=output+"BAM/"+sample, cns=output+"CNS/"+sample), shell=True) 
            subprocess.call(ivar_cns5 % dict(bam=output+"BAM/"+sample, cns=output+"CNS5/"+sample), shell=True) 
            os.remove(output+"CNS/"+sample+".qual.txt")
            os.remove(output+"CNS5/"+sample+".qual.txt")
    change_header(output+"CNS5/")
    change_header(output+"CNS/")
    

def align(s1_ref, s2_ref, s3_ref, output):
    subprocess.call("cat " + output+"CNS5/*REF_Sabin1* > " + output + "alignment/sabin1_not_aligned.fasta", shell=True)
    subprocess.call("cat " + output+"CNS5/*REF_Sabin2* > " + output + "alignment/sabin2_not_aligned.fasta", shell=True)
    subprocess.call("cat " + output+"CNS5/*REF_Sabin3* > " + output + "alignment/sabin3_not_aligned.fasta", shell=True)
    mafft(output+"alignment/sabin1_not_aligned.fasta", s1_ref, output+"alignment/sabin1_aligned.fasta")     
    mafft(output+"alignment/sabin2_not_aligned.fasta", s2_ref, output+"alignment/sabin2_aligned.fasta")
    mafft(output+"alignment/sabin3_not_aligned.fasta", s3_ref, output+"alignment/sabin3_aligned.fasta")
    
  #  subprocess.call("cat " + output+"CNS5/*REF_Poliovirus1* > " + output + "alignment/'wt1_not_aligned.fasta", shell=True)
   # subprocess.call("cat " + output+"CNS5/*REF_Poliovirus2-wt* > " + output + "alignment/wt2_not_aligned.fasta", shell=True)
   # subprocess.call("cat " + output+"CNS5/*REF_Poliovirus3-wt* > " + output + "alignment/wt3_not_aligned.fasta", shell=True)
   # mafft(output+"alignment/wt1_not_aligned.fasta", s1_ref, output+"alignment/wt1_aligned.fasta")     
   # mafft(output+"alignment/wt2_not_aligned.fasta", s2_ref, output+"alignment/wt2_aligned.fasta")
   # mafft(output+"alignment/wt3_not_aligned.fasta", s3_ref, output+"alignment/wt3_aligned.fasta")

def get_stat(depth_file, sample, stat, ref):
    
    mapped_reads = int(stat.loc[ref]["numreads"])
    covered_bases = int(stat.loc[ref]["covbases"])
    coverage = int(stat.loc[ref]["coverage"])
    
    depths = [int(x.split('\t')[2]) for x in open(depth_file).readlines()] if Path(depth_file).is_file() else ''
    mean_depth= str(round(mean(depths),2)) if depths else ''
    min_depth = min(depths) if depths else ''
    max_depth = max(depths) if depths else ''
    genome_size = len(depths)
    
    breadth_cns5 = len([i for i in depths if i > 5])
    coverage_cns5 = round(breadth_cns5 / genome_size * 100,2)  if genome_size else ''
    breadth_cns1 = len([i for i in depths if i > 1])
    coverage_cns1 = round(breadth_cns1 / genome_size * 100,2)  if genome_size else ''
    
    return mapped_reads, covered_bases, coverage, mean_depth, min_depth, max_depth, coverage_cns1, coverage_cns5


def QC_reports(barcodes, output):
    general_qc = pd.DataFrame(columns=["total_reads", "mapped_reads", "%mapped", "Sabin1","Sabin2","Sabin3","WT1","WT2","WT3", "%Sabin1 (from mapped_reads)", "%Sabin2 (from mapped_reads)", "%Sabin3 (from mapped_reads)","%wt1 (from mapped_reads)", "%wt2 (from mapped_reads)","%wt3 (from mapped_reads)"])
    sabin1_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    sabin2_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    sabin3_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    wt1_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    wt2_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    wt3_df = pd.DataFrame(columns=["mapped_reads", "covered_bases", "coverage", "mean_depth", "min_depth", "max_depth", "coverage", "coverage_cns5"])
    
    for sample in barcodes.values():
        
        #general information
        total_reads = pysam.AlignmentFile(output + "BAM/" + sample + ".bam").count(until_eof=True)  if Path(output + "BAM/" + sample + ".bam").is_file() else 0
        if total_reads == 0:
            continue
        if Path(output + "BAM/" + sample + ".REF_unmapped.bam").is_file():
            mapped_reads = total_reads - pysam.AlignmentFile(output + "BAM/" + sample + ".REF_unmapped.bam").count(until_eof=True)
        else:
            mapped_reads = total_reads
        mapped_percentage = round(mapped_reads / total_reads * 100, 2) 
        stat = pd.read_csv(StringIO(pysam.coverage(output + "BAM/" + sample + ".bam")), sep='\t').set_index("#rname")
        sabin1_reads = (stat.loc["Sabin1"]["numreads"])
        sabin1_reads_pre = round(sabin1_reads/mapped_reads*100,2)
        sabin2_reads = (stat.loc["Sabin2"]["numreads"])
        sabin2_reads_pre = round(sabin2_reads/mapped_reads*100,2)
        sabin3_reads = (stat.loc["Sabin3"]["numreads"])
        sabin3_reads_pre = round(sabin3_reads/mapped_reads*100,2)
       
        #merge all wt
        wt1_rows = [row for row in stat.index if "Poliovirus1-wt" in row]
        wt1 = stat.loc[wt1_rows].sum()
        wt1_reads = wt1["numreads"]
        wt1_reads_pre = round(wt1_reads / mapped_reads*100,2)
        wt1 = wt1.rename({"numreads": "mapped_reads", "covbases": "covered_bases","meandepth":"mean_depth"})
        
        
        wt2_rows = [row for row in stat.index if "Poliovirus2-wt" in row]
        wt2 = stat.loc[wt2_rows].sum()
        wt2_reads = wt2["numreads"]
        wt2_reads_pre = round(wt2_reads/mapped_reads*100,2)
        wt2 = wt2.rename({"numreads": "mapped_reads", "covbases": "covered_bases","meandepth":"mean_depth"})
        
        wt3_rows = [row for row in stat.index if "Poliovirus3-wt" in row]
        wt3 = stat.loc[wt3_rows].sum()
        wt3_reads = wt3["numreads"]
        wt3_reads_pre = round(wt3_reads / mapped_reads*100,2)
        wt3 = wt3.rename({"numreads": "mapped_reads", "covbases": "covered_bases","meandepth":"mean_depth"})
        
        general_qc.loc[sample] = (total_reads, mapped_reads, mapped_percentage,
                                  sabin1_reads, sabin2_reads, sabin3_reads, 
                                  wt1_reads, wt2_reads, wt3_reads,
                                  sabin1_reads_pre,
                                  sabin2_reads_pre, sabin3_reads_pre, wt1_reads_pre,
                                  wt2_reads_pre, wt3_reads_pre)
        
        # info by reference
        depth_file = output + "depth/" + sample +".REF_Sabin1.txt"
        sabin1_df.loc[sample] = get_stat(depth_file, sample, stat, "Sabin1" ) 
        depth_file = output + "depth/" + sample +".REF_Sabin2.txt"
        sabin2_df.loc[sample] = get_stat(depth_file, sample, stat, "Sabin2" )
        depth_file = output + "depth/" + sample +".REF_Sabin3.txt"
        sabin3_df.loc[sample] = get_stat(depth_file, sample, stat, "Sabin3" )
        wt1_df.loc[sample] = wt1
        wt2_df.loc[sample] = wt2
        wt3_df.loc[sample] = wt3
        
        
    with pd.ExcelWriter(output + "reports/" + 'reports.xlsx') as writer:
        general_qc.to_excel(writer, sheet_name='General_QC')
        sabin1_df.to_excel(writer, sheet_name='Sabin1')
        sabin2_df.to_excel(writer, sheet_name='Sabin2')
        sabin3_df.to_excel(writer, sheet_name='Sabin3')
        wt1_df.to_excel(writer, sheet_name='wt1')
        wt2_df.to_excel(writer, sheet_name='wt2')
        wt3_df.to_excel(writer, sheet_name='wt3')
    

def run(fastq_path, barcode_csv, output):
    par_dir = os.path.dirname(os.path.abspath(os.path.dirname(__file__ )))
    reference = par_dir + "/refs/VP1_sabins_WT.fasta"
    s1_ref = par_dir + "/refs/VP1_sabin1.fasta"
    s2_ref = par_dir + "/refs/VP1_sabin2.fasta"
    s3_ref = par_dir + "/refs/VP1_sabin3.fasta"
  #  create_dirs([fastq_path + "min_900/", output+"BAM/", output+"depth/", output+"CNS/", output+"CNS5/", output+"alignment/", output + "reports/"])
    barcodes = get_barcode_name(barcode_csv)
  #  filter_reads_by_length(fastq_path,fastq_path + "min_900/", )
    fastq_path = fastq_path + "min_900/"
   # minimap(fastq_path, reference, barcodes, output)
   # depth(output)
   # cns(output)
    align(s1_ref, s2_ref, s3_ref, output)
   # QC_reports(barcodes, output)
    
    run_ddns(output + "alignment/sabin1_aligned.fasta" ,par_dir + "/refs/sabin1.csv", output + "reports/sabin1_mutations.xlsx")
    run_ddns(output + "alignment/sabin2_aligned.fasta" ,par_dir + "/refs/sabin2.csv", output + "reports/sabin2_mutations.xlsx")
    run_ddns(output + "alignment/sabin3_aligned.fasta" ,par_dir + "/refs/sabin3.csv", output + "reports/sabin3_mutations.xlsx")
    
#    run_ddns(output + "alignment/wt1_aligned.fasta" ,par_dir + "/refs/sabin1.csv", output + "reports/wt1_mutations.xlsx")
 #   run_ddns(output + "alignment/wt2_aligned.fasta" ,par_dir + "/refs/sabin2.csv", output + "reports/wt2_mutations.xlsx")
  #  run_ddns(output + "alignment/wt3_aligned.fasta" ,par_dir + "/refs/sabin3.csv", output + "reports/wt3_mutations.xlsx")    

