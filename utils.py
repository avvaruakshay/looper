#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
import gzip
from itertools import takewhile, repeat

import sys

kmer_names = { 
    1: 'Monomer', 2: 'Dimer', 3: 'Trimer', 4: 'Tetramer', 5: 'Pentamer',
    6: 'Hexamer', 7: 'Heptamer', 8: 'Octamer', 9: 'Nonamer', 10: 'Decamer',
    11: 'Undecamer', 12: 'Dodecamer', 13: 'Tridecamer', 14: 'Tetradecamer',
    15: 'Pentadecamer', 16: 'Hexadecamer', 17: 'Heptadecamer', 18: 'Octadecamer',
    19: 'Nonadecamer', 20: 'Icosamer', 21: 'Uncosamer', 22: 'Docosamer',
    23: 'Tricosamer', 24: 'Tetracosamer', 25: 'Pentacosamer', 26: 'Hexacosamer',
    27: 'Heptacosamer', 28: 'Octacosamer', 29: 'Nonacosamer', 30: 'Triacontamer',
    31: 'Untriacontamer', 32: 'Dotriacontamer', 33: 'Tritriacontamer',
    34: 'Tetratriacontamer', 35: 'Pentatriacontamer', 36: 'Hexatriacontamer',
    37: 'Heptatriacontamer', 38: 'Octatriacontamer', 39: 'Nonatriacontamer',
    40: 'Tetracontamer', 41: 'Untetracontamer', 42: 'Dotetracontamer',
    43: 'Tritetracontamer', 44: 'Tetratetracontamer', 45: 'Pentatetracontamer',
    46: 'Hexatetracontamer', 47: 'Heptatetracontamer', 48: 'Octatetracontamer', 
    49: 'Nonatetracontamer', 50: 'Pentacontamer',
}


def rawcharCount(filename, char):
    if filename.endswith('gz'):
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(char.encode('ASCII')) for buf in bufgen if buf )


def get_cycles(string):
    cycles = set()
    for i in range(len(string)):
        cycles.add(string[i:] + string[:i])
    cycles = sorted(list(cycles))
    return cycles


def build_cycVariations(string):
    cycles = get_cycles(string)
    rev_cycles = get_cycles(rev_comp(string))
    for r in rev_cycles:
        if r not in cycles: cycles.append(r)
    return cycles


def getGC(basesCounter):
    totalBases = sum(basesCounter.values())
    try:
        GC = (float(basesCounter['G'] + basesCounter['C'])/(totalBases-basesCounter['N']))*100
    except KeyError:
        GC = (float(basesCounter['G'] + basesCounter['C'])/totalBases)*100
    return GC


def rev_comp(string):
    """Outputs reverse complement of a nucleotide sequence"""
    if sys.version_info.major == 2:
        import string as st
        complement = string.translate(st.maketrans('ACGT', 'TGCA'))
    else:
        complement = string.translate(str.maketrans('ACGT', 'TGCA'))
    return complement[::-1]
