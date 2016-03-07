#! /usr/bin/env python

import sys
import os
import os.path
import random

jp = os.path.join

# Enable if you want to see some debug output.
DEBUG = False

# V1 in paired domains.
# S1 in unpaired areas.
kReadLength = 60
kReadLengthVariability = 0
kMinReadLength = 50

kMinCoverage = 7
kMaxCoverage = 12

MODE = 'V1'

def readBpseq(bpseq_fn):
  """
  Read in a single file in BPSEQ format.
  """
  content = open(bpseq_fn).readlines()
  seq = [-1] * len(content)
  struct = [-1] * len(content)
  for i, entry in enumerate(content):
    pos, base, pair = entry.strip().split()
    seq[i] = base
    p = int(pair)
    struct[i] = [1, p][p == 0]
  return "".join(seq), struct


def generateReads(seq, struct, v1_fh, s1_fh):
  """
  """
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
      if (e - b) >= kMinReadLength:
        read = seq[b:e]
        if DEBUG:
          print "".join([['|', ' '][p==0] for p in struct])
          print dsym * (b-1) + read
        fh.write("%s\n" % read)


def processData(data_dir, output_dir):
  """
  """
  # V1 cuts at single-stranded positions
  v1_fh = open(jp(output_dir, "reads-V1.txt"), 'w+')
  # S1 cuts at single-stranded positions
  s1_fh = open(jp(output_dir, "reads-S1.txt"), 'w+')

  for bpseq_fn in os.listdir(data_dir):
    path = jp(data_dir, bpseq_fn)
    seq, struct = readBpseq(path)
    generateReads(seq, struct, v1_fh, s1_fh)


if __name__ == '__main__':
  processData(sys.argv[1], sys.argv[2])
