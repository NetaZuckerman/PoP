#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:03:13 2022

@author: hagar
"""
import os
import sh
import subprocess
import multiprocessing as mp
from Bio import SeqIO
import pandas as pd
import numpy as np

#mafft commands
MAFFT = os.path.dirname(__file__)+"/MAFFT.sh  %(not_aligned)s %(reference)s %(aligned)s"
RM_SPADES = "find spades -type f ! -name 'transcripts.fasta' -delete"
#split bam file
SPLIT = "bamtools split -in %(bam)s -reference"


ambiguous_nucleotides = ["W", "Y", "R", "S", "D","K","M","V","H","B","X"]

#return list of touple (R1,R2) file names
def get_r1r2_list(fastq_path):
    r1r2_list = []
    skip_files=["R2", "Undetermined", "unpaired", "singletons"]
    for r1 in os.listdir(fastq_path):
        if any(skip_file in r1 for skip_file in skip_files) or "fast" not in r1:
            continue
        r2 = r1.replace("R1","R2")
        sample = r1.split("_")[0].split(".fastq")[0] #sample short name
        r1r2_list.append([sample,r1,r2])
    return r1r2_list
 
#split all bam files by segments
def split_bam(dir):
    for bam_file in os.listdir(dir):        
        if "sorted" in bam_file and "bai" not in bam_file:
            subprocess.call(SPLIT % dict(bam=dir+bam_file), shell=True)
            os.remove(dir + bam_file)
  
#removes dir is exists and recreates it
def create_dirs(dirs):
        for dir in dirs:
            if os.path.exists(dir):
                continue
            os.makedirs(dir)
            
def remove_from_name(dir, to_remove):
    for file in os.listdir(dir):
        file = dir + file
        os.rename(file, file.replace(to_remove, ''))

#change first line in all files in a directory to ">sample_name"
def change_header(dir):
    for file in os.listdir(dir):
        if ".fa" in file:
            new_header = ">" + file.split(".fa")[0].split(".REF")[0]
            file = dir + file        
            sh.sed("-i", "1s/.*/" + new_header + "/", file)
        

def mafft(not_aligned, reference, aligned):
    subprocess.call(MAFFT % dict(not_aligned=not_aligned, reference=reference, aligned=aligned), shell=True)

def rm_spades():
    subprocess.call(RM_SPADES, shell=True)

def get_sequences(alignment_file):
    sequences = {}
    alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
    for sample, record in alignment.items():
        sequences[sample] = str(record.seq).upper()

        sequences[sample.replace("Consensus_", "").split("threshold")[0]] = sequences.pop(sample)

    return sequences




#get mutations positions list by comparing all sequences to each other.
#"-" and "N" are excluded.
def mutations_positions(sequences, no_n = 0):
    mutations_positons = []
    seq_length = len(next(iter(sequences.values())))
    for pos in range(seq_length-1):
        temp = ""
        for sample, record in sequences.items():
            if not temp:
                temp = record[pos]
            if no_n and record[pos] in ["N", "-"]:
                break
            if record[pos] in ambiguous_nucleotides:
                break
            if not temp == record[pos] :
                mutations_positons.append(pos)
                break
        continue
    return mutations_positons

def run_mp(threads, func, arg):
    with mp.Pool(threads) as pool:
        pool.map(func,arg)
        pool.close()
        pool.join()
 
def hamming_distance(seq1, seq2):
 df = pd.DataFrame()
 df["seq1"] = pd.Series(seq1)
 df["seq2"] = pd.Series(seq2)
 df['difference'] = np.where((df["seq1"] == df["seq2"]) | (df["seq1"] == "N") | (df["seq2"] == "N") | (df["seq1"] == "-") | (df["seq2"] == "-"), 0, 1)
 
 return (df["difference"].sum())


translate_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}
