#! /usr/bin/env python

import os
import sys
import glob
import random
import math
path = 'partial_bpseq_rfam-12.0/'
jp = os.path.join

def readBpseq(bpseq_fn):
    """
    Reads sequence information. 
    """
    content = open(bpseq_fn).readlines()
    seq = [-1] * len(content)
    struct = [-1] * len(content)
    pairing = [-1] * len(content)
    for i, entry in enumerate(content):
        pos, base, pair = entry.strip().split()
        seq[i] = base
        p = int(pair)
        struct[i] = [1, p][p == 0]
        pairing[i] = p
    return "".join(seq), struct, pairing


kReadLength = 60
kReadLengthVariability = 0
kMinReadLength = 20
kMinCoverage = 7
kMaxCoverage = 12
MODE = 'V1'

def generateReads(seq, struct, s1_fh, v1_fh):
    """
    Read generator.
    """
    print('lenseq', len(seq))
    for i, s in enumerate(struct):
        if s == 1:
            dsym = 'v'
            fh = v1_fh
        else:
            fh = s1_fh
            dsym = 's'
        for r in range(random.randint(kMinCoverage, kMaxCoverage)):
            b = i + 1
            e = min(b + kReadLength + random.randint(-kReadLengthVariability, kReadLengthVariability), len(seq))
            if (e-b) >= kMinReadLength:
                read = seq[b:e]
                fh.write("%s\n" % read)
                
def ReadCoverageCalculation(ReadsS1_fn, Sequence):
    """
    Read Coverage for S1. 
    """
    ReadCoverageS1 = [0] * len(Sequence)
    contentS1 = open(ReadsS1_fn).readlines()
    for i, read in enumerate(contentS1):
        read2 = read.strip()
        if read2 in Sequence:
            ReadCoverageS1[Sequence.find(read2)+1] += 1
    ReadCoverageN = NormaliseReadCoverage(ReadCoverageS1)
    return ReadCoverageN;

def NormaliseReadCoverage(coverage):
    scaled_coverage = []
    for i in range(len(coverage)):
        if min(coverage) != max(coverage):
            scaled_coverage.append((coverage[i]-min(coverage))/(max(coverage)-min(coverage)))
        else:
            scaled_coverage.append((coverage[i]-min(coverage))/(max(coverage)-min(coverage)+0.001))
    return scaled_coverage;

# Generate input for run_qrisp.sh file.
def NucTransform(nuc):
    if nuc=='A':
        return 0
    if nuc=='C':
        return 1
    if nuc=='U':
        return 2
    if nuc=='G':
        return 3

allfiles = glob.glob(path + 'RF*.bpseq')
output_dir = ''
training_data = open(jp(output_dir, "Rfam-trainingdata.txt"), "w")
ReadCoverage_list = []
StructTags = []
training_test_set = ["RF00012_B", "RF00026_B", "RF00045_A", "RF00050_B", "RF00059_B", "RF00065_A",
                     "RF00083_B", "RF00102_B", "RF00114_B", "RF00128_B", "RF00140_A", "RF00168_B",
                     "RF00182_B", "RF00389_A", "", "RF00433_B", "RF00461_A", "RF00484_A",
                     "RF00511_B", "RF02532_A",
                     "RF02530_B", "RF02531_A", "RF02536_A", "RF02537_A", "RF02538_B", "RF01794_A",
                     "RF01796_B", "RF01518_B",
                     "RF02359_A", "RF02358_B", "RF02534_B", "RF01518_B"]
v1_fh = open(jp(output_dir, "reads-V1.txt"), 'w+')
s1_fh = open(jp(output_dir, "reads-S1.txt"), 'w+')

for i in range(len(allfiles)):
    tag = allfiles[i][len(allfiles[i])-15:len(allfiles[i])-6]
    if tag in training_test_set:
        seq, struct, pairing = readBpseq(allfiles[i])
        generateReads(seq, struct, s1_fh, v1_fh)
        ReadCoverage1 = ReadCoverageCalculation(output_dir+"reads-S1.txt", seq)
        ReadCoverage_list.append(ReadCoverage1)
        training_data.write('item { \n')
        training_data.write("\tid: '"+tag+"'\n")
        StructTags.append(tag)
        training_data.write('\tstructure { \n')
        for j in range(len(seq)):
            data = ('\t\trows { pos: '+str(j+1)+' pair: '+str(pairing[j])+' base: '+str(NucTransform(seq[j]))+' score: '+str(ReadCoverage1[j])+'}')
            training_data.write(data+' \n')
        training_data.write('\t} \n')
        training_data.write('} \n')
training_data.close()
