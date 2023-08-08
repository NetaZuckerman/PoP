#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 14:06:34 2022

@author: hagar

----------------------------------------------------------------------------------------------
This code produces mutations report
input:
    argv[1] = alignment file
    argv[2] = regions file (output of parse gb file script)
    argv[3] = output file (path+name)
    
The code will find the mutations of the sequences and present it ordered by genes.
----------------------------------------------------------------------------------------------

"""
import sys
sys.path.insert(1, '/home/hagar/PoP/')
from sys import argv
from math import floor
import pandas as pd
from Bio.Seq import Seq
from utils.utils import translate_table, get_sequences, mutations_positions
from utils.format_xl import save_format_xl
ambiguous_nucleotides = ["W", "Y", "R", "S", "D","K","M","V","H","B","X"]


def mutations_by_sample(mutations_position,sequences):
    '''
    Parameters
    ----------
    mutations_position : list
        all mutations positions in all sequences.
    sequences : dict
        {sample : fasta(str)}.

    Returns
    -------
    mutations_by_sample : dict
        {sample : nucleotide list}. the dict value is a list of the nucletide in the positions of mutations_position list

    '''
    mutations_by_sample = {}
    for sample, record in sequences.items():
        mutations = []
        for pos in mutations_position:
            mutations.append(record[pos])
        mutations_by_sample[sample] = mutations
    return mutations_by_sample


def get_regions(regions_csv):
    '''

    Parameters
    ----------
    regions_csv : str
        path to regions file (output of parse gb file script).

    Returns
    -------
    regions : dict
        {gene : (start, end)}.

    '''
    regions = {}
    f = open(regions_csv)
    for line in f.readlines():
        if line.upper().startswith("GENE"):
            continue
        line = line.split(",")
        regions[line[0]] = (int(line[1]),int(line[2]),line[3])
    f.close()
    return regions


def get_gene(mutations_positions_nt, regions):
    '''
    iterate all mutations positions and get the gene and the position on the gene (amino acid and nucleotide)
    Parameters
    ----------
    mutations_positions_nt : list
        mutations positions on nucleotides.
    regions : dict
        {gene : (start, end)}.

    Returns
    -------
    gene_names : list
        gene names list in the size of the mutations number
    position_on_gene_nt : list
    position_on_gene_aa : list

    '''
    gene_names = []
    position_on_gene_nt = []
    position_on_gene_aa = []
    for mut in mutations_positions_nt:
        gene_name = "UTR"
        pos_gene_nt = ""
        pos_gene_aa = ""
        for gene, pos in regions.items():
            start = pos[0]
            end = pos[1]
            if mut in range(start,end):
                gene_name = gene
                pos_gene_nt = mut-start + 2
                pos_gene_aa = floor(pos_gene_nt/3) if pos_gene_nt%3 == 0 else  floor(pos_gene_nt/3) + 1
        gene_names.append(gene_name)
        position_on_gene_nt.append(pos_gene_nt)
        position_on_gene_aa.append(pos_gene_aa)

    return gene_names,position_on_gene_nt,position_on_gene_aa 

def aa_sum(df, sequences):
    '''
    summarize the amino acid properties.
    find out if a mutation is is a replacement(R) or seilent(S) by iteratind the amino acids.
    if it is a replacement find the amino acids properties and compare it with each other. 

    Parameters
    ----------
    df : DataFrame
        the mutations report.
    sequences : dict
        {sample : fasta}.

    Returns
    -------
    None. it updates the results in the DataFrame

    '''
    aa_groups = pd.read_csv("/home/hagar/UPv/mutations/AAproperties.txt", sep = '\t')
    for index, row in df.iloc[:,-(len(sequences)):].iterrows():
        row.reset_index(drop=True,inplace=True)
        no_x_aa = list(dict.fromkeys([x for x in row.to_list() if x != "X"]))
        if len(no_x_aa) > 1:
            df.at[index,"R/S"] =  "R"
            if len(no_x_aa) == 2:
                groups = aa_groups.loc[aa_groups['Abbv1'] == row.iloc[0], 'properties'].iloc[0] + "(" + row.iloc[0] +"), " 
                no_x_aa.remove(row.iloc[0])
                groups += aa_groups.loc[aa_groups['Abbv1'] == no_x_aa[0], 'properties'].iloc[0] + "(" + no_x_aa[0] +")"
                df.at[index, "aa_group"] = groups
            if len(no_x_aa) > 2:
                groups = aa_groups.loc[aa_groups['Abbv1'] == row.iloc[0], 'properties'].iloc[0] + "(" + row.iloc[0] +"), " 
                no_x_aa.remove(row.iloc[0])
                groups += aa_groups.loc[aa_groups['Abbv1'] == no_x_aa[0], 'properties'].iloc[0] + "(" + no_x_aa[0] +")"
                no_x_aa.remove(no_x_aa[0])
                groups += aa_groups.loc[aa_groups['Abbv1'] == no_x_aa[0], 'properties'].iloc[0] + "(" + no_x_aa[0] +")"
                df.at[index, "aa_group"] = groups
        else:
            df.at[index,"R/S"] =  "S"


def get_single_aa(seq, position, region):
    '''
    find the translation of a codon by the reading frame of the gene.

    Parameters
    ----------
    seq : string
        fasta sequence.
    position : int
        the nucleotide position that needs to be translated.
    start : int
        the gene start - indecates the reading frame.

    Returns
    -------
    aa : chr
        the translated amino acid.

    '''
    start = region[0]
    end = region[1]
    strand= region[2]
    
    

    #get codon
    if strand == "-":
        pos_on_gene = end - position + 1
        mod = pos_on_gene % 3
        if mod == 0:
            codon_pos = (position + 2, position + 1, position)
        if mod == 1:
            codon_pos = (position, position - 1, position - 2)
        if mod == 2:
            codon_pos = (position + 1, position, position - 1)
            
        codon = seq[codon_pos[0]-1] + seq[codon_pos[1]-1] + seq[codon_pos[2]-1]
        codon = str(Seq(codon).complement())
    else:
        
        pos_on_gene = position - start + 1
        
        mod = pos_on_gene % 3
        if mod == 0:  # third nuc on the codon
            codon_pos = (position -2 , position - 1, position)
        if mod == 1:  # first nuc on the codon
            codon_pos = (position , position + 1, position + 2)
        if mod == 2:  # second nuc on the codon
            codon_pos = (position - 1, position, position + 1)
            
        codon = seq[codon_pos[0]] + seq[codon_pos[1]] + seq[codon_pos[2]]
    
    if codon in translate_table:
        aa = translate_table[codon] if not '-' in codon and not 'N' in codon else 'X'
    else:
        aa = 'X'
    return aa


def get_all_aa(mutations_positions_nt, sequences, gene_names, regions):
    '''
    
    Parameters
    ----------
    mutations_positions_nt : list
        nucleotide mut position.
    sequences : dict
    gene_names : list
    regions : dict

    Returns
    -------
    mutations_by_sample_aa : dict
        {sample : list of the amino acids in mutations positions}.

    '''    
    mutations_by_sample_aa = {}
    for sample, seq in sequences.items():
        sample_aa = []
        for i in range(len(mutations_positions_nt)):
            if gene_names[i].endswith('UTR'):
                sample_aa.append('X')
            else:
                pos = mutations_positions_nt[i]
                region = regions[gene_names[i]]
                aa = get_single_aa(seq, pos, region)
                sample_aa.append(aa)
        mutations_by_sample_aa[sample] = sample_aa

    return mutations_by_sample_aa
    


def run_ddns(alignment_file,regions_csv,output):
        
    
    '''
    run all functions.

    '''
    sheet_name = regions_csv.split("/")[-1].split(".csv")[0]
    sequences = get_sequences(alignment_file)
    mutations_positions_nt = mutations_positions(sequences,1)
    mutations_by_sample_nt = mutations_by_sample(mutations_positions_nt,sequences)
    regions = get_regions(regions_csv)
    gene_names, position_on_gene_nt, position_on_gene_aa = get_gene(mutations_positions_nt, regions)
    get_all_aa(mutations_positions_nt, sequences, gene_names, regions)
    
    mutations_by_sample_aa = get_all_aa(mutations_positions_nt, sequences, gene_names, regions)
    
    
    df = pd.DataFrame()
    df["nt_position_on_gene"] = position_on_gene_nt
    df.insert(0, "gene_name", "VP1")
    df["nt_position_on_genome"] = mutations_positions_nt 
    df["nt_position_on_genome"] =df["nt_position_on_genome"] + 1 + regions["VP1"][0]
    for sample, mut in mutations_by_sample_nt.items():
        df[sample+"_NT"] = mut
    df["aa_position_on_gene"] = position_on_gene_aa
    for sample, mut in mutations_by_sample_aa.items():
        df[sample+"_AA"] = mut    
    
    aa_sum(df, sequences)
    save_format_xl(df, len(sequences)-1, output, sheet_name)
    
    
def run(alignment_file,regions_csv,output, noN = 0):
            
    
    '''
    run all functions.
    Parameters:
        noN: ignore N's and gups in mutation analysis if noN=1. 
        default = 0
            
        analysing N's and gaps option is made for nOPV2.
        why? nOPV2 5'UTR is longer than the other polioviruses
        if you see N's in this region you should suspect its nOPV2
    '''
    sheet_name = regions_csv.split("/")[-1].split(".csv")[0]
    sequences = get_sequences(alignment_file)
    mutations_positions_nt = mutations_positions(sequences,noN)
    mutations_by_sample_nt = mutations_by_sample(mutations_positions_nt,sequences)
    regions = get_regions(regions_csv)
    gene_names, position_on_gene_nt, position_on_gene_aa = get_gene(mutations_positions_nt, regions)
    get_all_aa(mutations_positions_nt, sequences, gene_names, regions)
    
    mutations_by_sample_aa = get_all_aa(mutations_positions_nt, sequences, gene_names, regions)
    
    
    df = pd.DataFrame()
    df["gene_name"] = gene_names
    df["nt_position_on_gene"] = position_on_gene_nt
    df["nt_position_on_genome"] = mutations_positions_nt
    df["nt_position_on_genome"] += 1
    for sample, mut in mutations_by_sample_nt.items():
        df[sample+"_NT"] = mut
    df["aa_position_on_gene"] = position_on_gene_aa
    for sample, mut in mutations_by_sample_aa.items():
        df[sample+"_AA"] = mut    
    
    aa_sum(df, sequences)
    save_format_xl(df, len(sequences)-1, output, sheet_name)

    
    
if __name__ == "__main__":
    run(argv[1], argv[2], argv[3])
    