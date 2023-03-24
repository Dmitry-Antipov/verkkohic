#!/usr/bin/python3

import sys
import os
import random
from os import listdir
from os.path import isfile, join

def evaluate_set(contig_set, lengths, colors):
    m_len = 0
    p_len = 0
    u_len = 0
    for contig in contig_set:
        if contig in lengths:
            l = lengths[contig]
            if not (contig in colors):
                u_len += l
            else:
                if colors[contig] == 'm':
                    m_len += l
                elif colors[contig] == 'p':
                    p_len += l
    total_l = m_len + p_len + u_len
    if total_l > 1000000:
        print(f'{total_l}    {m_len / total_l:.3f}/{p_len / total_l:.3f}/{u_len / total_l:.3f}')
        if m_len > 0 and p_len > 0:
            print('BAD')
#pat_from_utig4-2246     utig4-835+,utig4-2245+,utig4-2246+,utig4-2520+,utig4-2521+      PATERNAL

def get_phased_edges(phasedfile):
    phased = set()
    for line in open(phasedfile, 'r'):
        arr = line.split()
        name = arr[0]
        color = "a"
        if arr[4] == "#8888FF" or arr[4] == "#FF8888":
            phased.add(name)
    return phased

#Evaluate paths using only phased edges. If empty set passed - all edges are used.
def evaluate_rukki(rukkifile, gfafile, triofile, phased_edges, out_f):
    lengths = {}
    for line in open(gfafile, 'r'):
        arr = line.split()
        if arr[0] == "S":
            lengths[arr[1]] = int(arr[3].split(':')[-1])
    colors = {}
    for line in open(triofile, 'r'):
        arr = line.split()
        name = arr[0]
        color = "a"
        if arr[4] == "#8888FF":
            color = "m"
        elif arr[4] == "#FF8888":
            color = "p"
        colors[name] = color

    if len(phased_edges) == 0:
        for c in colors.keys():
            phased_edges.add(c)
    unassigned = 0
    assigned = 0
    errors = 0
    error_l = 0
    total_l = 0
    for line in open(rukkifile):
        strpath = line.split()[1]
        path = strpath.split(',')
        state = "0"
        prev_contig = ""
        hap_l = [0,0]
        for sp in path:
            p = sp[:-1]
            if p in colors:
                if colors[p] == "a" or not p in phased_edges:
                    unassigned += 1
                    continue
                else:
                    if colors[p] != state and state != "0":
                        out_f.write(f"Discordant colors between {prev_contig} {p} !!!\n")
                        out_f.write(strpath +"\n")
                        errors += 1
                    assigned += 1            
                    prev_contig = p
                    state = colors[p]
                    if state == 'p':
                        hap_l[0] += lengths[p]
                    else:
                        hap_l[1] += lengths[p]
        error_l += min(hap_l[0], hap_l[1])
        total_l += hap_l[0] + hap_l[1]
    out_f.write(f"Among contigs in paths {rukkifile}, using uncolored/colored {unassigned}/{assigned} edges, we see {errors} errors\n")
    out_f.write(f"Hamming error rate estimate for {rukkifile}: {error_l/total_l:.4f}\n")
if __name__ == "__main__":                
    if len(sys.argv) < 4:
        print(f'Usage: {sys.argv[0]} <rukkifile.tsv> <gfa file> <trio.csv> [binned edges csv]')
        exit()
    if len(sys.argv) == 4:
        evaluate_rukki(sys.argv[1], sys.argv[2], sys.argv[3], set(), sys.stdout)
    else:
        classified = set()
        for line in open (sys.argv[4],'r'):
            classified.add(line.split()[0])
        evaluate_rukki(sys.argv[1], sys.argv[2], sys.argv[3], classified, sys.stdout)

