#!/usr/bin/env python
import os, argparse
from os.path import splitext

from analyse import generate_defaultInfo

"""Edits the ebeta.cpp template to compile user inputs as constants"""


def motif_factors(m, M):
    factors = []
    for i in range(M, m-1, -1):
        check = 0
        for j in factors:
            if j % i == 0: check = 1; break
        if check == 0: factors.append(i)
    
    return factors

def getArgs():
    """
    Parses command line arguments and returns them to the caller
    """
    __version__ = 'v0.1.0'
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='<FILE>', help='Input sequence file.')
    
    optional = parser.add_argument_group('Optional arguments')
    
    #Basic options
    optional.add_argument('-o', '--output', metavar='<FILE>', help='Output file name. Default: Input file name + _perf.tsv')
    # optional.add_argument('--format', metavar='<STR>', default='fasta', help='Input file format. Default: fasta, Permissible: fasta, fastq')
    # optional.add_argument('--version', action='version', version='looper ' + __version__)
        
    #Selections options based on motif size and seq lengths
    optional.add_argument('-m', '--min-motif-size', type=int, metavar='<INT>', default=1, help='Minimum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-M', '--max-motif-size', type=int, metavar='<INT>', default=6, help='Maximum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-l', '--min-length', type=int, metavar='<INT>', help='Minimum length cutoff of repeat')
    
    # Analysis options
    optional.add_argument('-a', '--analyse', action='store_true', default=False, help='Generate a summary HTML report.')

    args = parser.parse_args()

    if args.min_length is None:
        args.min_length = 2 * args.max_motif_size
    
    if args.output is None:
        args.output = splitext(args.input)[0] + '_looper.tsv'

    return args

def main():
    """Main function of looper"""


    args = getArgs()

    m = args.min_motif_size
    M = args.max_motif_size
    cutoff = args.min_length

    motif_checks = motif_factors(m, M)
    N = len(motif_checks)
    motif_checks_string = '{'
    for a in motif_checks: motif_checks_string += str(a) + ','
    motif_checks_string = motif_checks_string[:-1] + '}'

    divisor = []
    rem_shift = []
    for i, a in enumerate(motif_checks): 
        d = cutoff // motif_checks[i]
        r = cutoff % motif_checks[i]
        D = ('0'*((2*a)-1) + '1')*d + '0'*(2*r)
        divisor.append(int(D, 2))
        rem_shift.append( int(2*(cutoff - (cutoff % motif_checks[i]))) )

    div_string = f'uint64_t divisor[{N}] = ' + '{'
    for a in divisor: div_string += str(a) + ','
    div_string = div_string[:-1] + '}'

    rem_shift_string = f'uint rem_shift[{N}] = ' + '{'
    for a in rem_shift: rem_shift_string += str(a) + ','
    rem_shift_string = rem_shift_string[:-1] + '}'

    script = open('./pylooper.cpp', 'w')
    with open('./pylooper_template.cpp') as fh:
        for line in fh:
            line = line.rstrip()
            if '$' in line:
                prefix = line[:line.find('$')]
                line = f'{prefix}uint cutoff = {cutoff};\
                    \n{prefix}uint m = {m};\
                    \n{prefix}uint M = {M};\
                    \n{prefix}uint motif_checks[{N}] = {motif_checks_string};\
                    \n{prefix}uint N = {N};\
                    \n{prefix}{div_string};\
                    \n{prefix}{rem_shift_string};'
            print(line, file=script)
    script.close()
    
    os.system('g++ ./pylooper.cpp -O3 -o pylooper')
    os.system(f'./pylooper  {args.input} {args.output}')

    if args.analyse: generate_defaultInfo(args)

if __name__ == "__main__":
    main()
