#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 08:59:48 2023

@author: hagar
"""
import argparse
from os import getcwd

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='fastq folder')
    parser.add_argument('-r','--reference' ,help='reference')
    parser.add_argument('-o','--output_dir' ,help='output_dir', default=getcwd() + "/PoP_output/")
    parser.add_argument('-b','--barcodes' ,help='barcodes.csv file for nanopore output')
    parser.add_argument('--ddns', action='store_true', help="PolioVirus analysis") #store_true will store the argument as true
    parser.add_argument('--nopv', action='store_true', help="get nOPV attenuaions and mutations table. input should be fasta file")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    run()