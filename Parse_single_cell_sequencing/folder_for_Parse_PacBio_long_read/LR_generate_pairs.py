import os
import re
import gzip
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from fuzzysearch import find_near_matches

parser = argparse.ArgumentParser(description="Detect and change reverse linker seq")
parser.add_argument("--fastq",help="fastq input")
parser.add_argument("--new_fname",help="New filename")
parser.add_argument("--out_dir",help="Output directory")
parser.add_argument("--l1dist",type=int,help="Edit distance for first linker")
parser.add_argument("--l2dist",type=int,help="Edit distance for second linker")
parser.add_argument("--chemistry",help="Kit chemistry version (e.g. v1, v2, v3)")
parser.add_argument("--multiple_fq",action="store_true",help="Multiple fastq files")
parser.add_argument("--max_deletions",default=1,type=int,help="Max deletions for edit distance")
args = parser.parse_args()

chem = args.chemistry
out_dir = args.out_dir
filename = args.fastq
l1dist = args.l1dist
l2dist = args.l2dist
new_fname = args.new_fname
multiple = args.multiple_fq
max_d = args.max_deletions

if chem == "v1" or chem == "v2":
    L1F = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
    L2F = "ATCCACGTGCTTGAGACTGTGG"
    L1RC = "AGTCGTACGCCGATGCGAAACATCGGCCAC"
    L2RC = "CCACAGTCTCAAGCACGTGGAT"
elif chem == "v3":
    L1F = "ATGAGGGGTCAG"
    L2F = "TCCAACCACCTC"
    L1RC = "CTGACCCCTCAT"
    L2RC = "GAGGTGGTTGGA"
else:
    print("chemistry not set")  # Assuming you want to print this message

r1fwd_list=[]
r2fwd_list=[]
r1rev_list=[]
r2rev_list=[]

if multiple:
    set_num = re.findall("_[0-9]*\.", filename)[0]
else:
    set_num = "."

cnts = {"n_reads":0, "rc_match":0,"fwd_match":0, "no_match":0,
    "no_bc3":0, "rev_bad_pos":0, "fwd_bad_pos":0}

def clean_header(record):
    record.description = record.description.split(" ")[0]

"""
Max deletions has large impact on processing time.
It's generally better to not go above 1.
"""

def fnm(L1RC, max_l_dist):
    record_seq = record.seq
    return find_near_matches(L1RC, record_seq, max_l_dist=max_l_dist, max_deletions=max_d)

if not os.path.exists(filename):
    print(f"File: {e} doesn't exist, skipping")
else:
    try:
        with gzip.open(filename, "rt") as handle:
            records = SeqIO.parse(handle, "fastq")
            for record in records:
                
                l1rc_match = fnm(L1RC,l1dist)
                l2rc_match = fnm(L2RC,l2dist)
                l1f_match = fnm(L1F,l1dist)
                l2f_match = fnm(L2F,l2dist)
                
                if len(l1rc_match) and len(l2rc_match) > 0:
                    if l1rc_match[0].start - l2rc_match[0].end == 8:
                        clean_header(record)
                        end = l1rc_match[0].end + 18
                        r1rev = record[1:(end - 85)]
                        r2rev = record[(end - 86):end]
                        if len(r2rev) < 86:
                            cnts["rev_bad_pos"] += 1
                        elif end <= len(record):
                            r2rev.seq = r2rev.seq.reverse_complement()
                            r1rev_list.append(r1rev)
                            r2rev_list.append(r2rev)
                            cnts["rc_match"] += 1
                        else:
                            cnts["no_bc3"] += 1
                        
                elif len(l1f_match) and len(l2f_match) > 0:
                    if l2f_match[0].start - l1f_match[0].end == 8:
                        clean_header(record)
                        start = l1f_match[0].start - 18
                        r1fwd = record[(start + 86):len(record)]
                        r2fwd = record[start:(start + 86)]
                        
                        if len(r2fwd) < 86:
                            cnts["fwd_bad_pos"] += 1
                        elif start > 0:
                            r1fwd.seq = r1fwd.seq.reverse_complement()
                            r1fwd_list.append(r1fwd)
                            r2fwd_list.append(r2fwd)
                            if len(r2fwd) == 0:
                                break
                            cnts["fwd_match"] += 1
                        else:
                            cnts["no_bc3"] += 1
                else:
                    cnts["no_match"] += 1
                    SeqIO.write(record, out_dir + "unmatched.fastq", "fastq")
                
                cnts["n_reads"] += 1
                if cnts["n_reads"] % 1000 == 0:
                    now = datetime.now()
                    curr_time = now.strftime("%H:%M:%S")
                    print(cnts["n_reads"], "reads processed", curr_time)
            
            r1_complete = r1fwd_list + r1rev_list
            r2_complete = r2fwd_list + r2rev_list
            out_r1 = out_dir + f"{new_fname}_R1" + set_num + "fastq"
            out_r2 = out_dir + f"{new_fname}_R2" + set_num + "fastq"
            
            SeqIO.write(r1_complete, out_r1, "fastq")
            SeqIO.write(r2_complete, out_r2, "fastq")
            
            print(filename, out_r1, out_r2)
            print(cnts)
            print("Fraction of reads retained",(cnts["rc_match"]+cnts["fwd_match"])/cnts["n_reads"], "\n")
    except Exception as exptn:
        print(f"An error occurred while processing the file: {exptn}")
