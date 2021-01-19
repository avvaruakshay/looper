#!/usr/bin/env python

"""
@author: Avvaru Akshay
@filename: looper.py

Looper python script.

This code calculates and inserts the desired parameters in C++ code as constants.
Thus the C++ code doesn't lag on calculating the parameters during the runtime
and also optimises the divisions as the divisors are known at the compile time.
Compiles the generated c++ code and runs repeat indentification using the
compiled c++ code on input sequence file.

"""

# Future
from __future__ import print_function, division

# Generic/Built-in
import os, argparse, sys
from os.path import splitext

# Owned
from analyse import analyse_fasta, analyse_fastq
from annotation import annotate_repeats

if sys.version_info[0] == 2:
    pass


def motif_factors(m, M):
    """
    Generates the set of non redundant motif sizes to be checked for division
    rule.

    Parameters
    ----------
    m : <int> minimum motif size

    M : <int> maximum motif size

    Returns
    -------
    factors : <list> sorted(descending) list of non-redundant motifs.
    """
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
    optional.add_argument('--format', metavar='<STR>', default='fasta', help='Input file format. Default: fasta, Permissible: fasta, fastq')
    optional.add_argument('-v-', '--version', action='version', version='looper ' + __version__)
        
    #Selections options based on motif size and seq lengths
    optional.add_argument('-m', '--min-motif-size', type=int, metavar='<INT>', default=1, help='Minimum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-M', '--max-motif-size', type=int, metavar='<INT>', default=6, help='Maximum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-l', '--min-length', type=int, metavar='<INT>', help='Minimum length cutoff of repeat')
    optional.add_argument('--filter-reads', action='store_true', default=False, help='Seperate out the reads with repeats in a new file.')
    
    # Analysis options
    optional.add_argument('-a', '--analyse', action='store_true', default=False, help='Generate a summary HTML report.')


    # Annotation options
    annotation = parser.add_argument_group('Annotation arguments')
    annotation.add_argument('-g', '--annotate', metavar='<FILE>', help='Genic annotation input file for annotation, Both GFF and GTF can be processed. Use --anno-format to specify format.')
    annotation.add_argument('--anno-format', metavar='<STR>',default='GFF', type=str, help='Format of genic annotation file. Valid inputs: GFF, GTF. Default: GFF')
    annotation.add_argument('--gene-key', metavar='<STR>', default='gene', type=str, help='Attribute key for geneId. The default identifier is "gene". Please check the annotation file and pick a robust gene identifier from the attribute column.')
    annotation.add_argument('--up-promoter', metavar='<INT>', type=int, default=1000, help='Upstream distance(bp) from TSS to be considered as promoter region. Default 1000')
    annotation.add_argument('--down-promoter', metavar='<INT>', type=int, default=1000, help='Downstream distance(bp) from TSS to be considered as promoter region. Default 1000')

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

    div_string = 'uint64_t divisor[{N}] = '.format(N=N) + '{'
    for a in divisor: div_string += str(a) + ','
    div_string = div_string[:-1] + '}'

    rem_shift_string = 'uint rem_shift[{N}] = '.format(N=N) + '{'
    for a in rem_shift: rem_shift_string += str(a) + ','
    rem_shift_string = rem_shift_string[:-1] + '}'

    template_file = ''
    filtered_file = ''
    if args.format == 'fasta':
        template_file = './pylooper_fasta_template.cpp'
    elif args.format == 'fastq':
        if args.input.endswith('.gz'):
            template_file = './pylooper_fastq_gzip_template.cpp'
        else:
            template_file = './pylooper_fastq_template.cpp'
        if args.filter_reads: 
            filtered_file = splitext(args.input)[0] + '_looper.filtered.fastq'

    script = open('./pylooper.cpp', 'w')
    with open(template_file) as fh:
        for line in fh:
            line = line.rstrip()
            if '$' in line:
                prefix = line[:line.find('$')]
                line = '{prefix}uint cutoff = {cutoff};'.format(prefix=prefix, cutoff=cutoff)
                line += '\n{prefix}uint m = {m};'.format(prefix=prefix, m=m)
                line += '\n{prefix}uint M = {M};'.format(prefix=prefix, M=M)
                line += '\n{prefix}uint motif_checks[{N}] = {motif_checks_string};'.format(prefix=prefix, N=N, motif_checks_string=motif_checks_string)
                line += '\n{prefix}uint N = {N};'.format(prefix=prefix, N=N)
                line += '\n{prefix}{div_string};'.format(prefix=prefix, div_string=div_string)
                line += '\n{prefix}{rem_shift_string};'.format(prefix=prefix, rem_shift_string=rem_shift_string)
            print(line, file=script)
    script.close()
    status = os.system('g++ ./pylooper.cpp -O3 -o pylooper')
    if sys.platform.startswith('linux') and status != 0:
        print('\n\nError compiling the cpp auxiliary code. Please check for proper\
               installation of g++. ')
        sys.exit()
    
    if args.format == 'fastq' and args.filter_reads:
        os.system('./pylooper  {input} {output} {filtered_file}'.format(input=args.input, output=args.output, filtered_file=filtered_file))
    if args.format == 'fastq' and args.input.endswith('.gz'):
        os.system('zcat {input} | ./pylooper {output} '.format(input=args.input, output=args.output))
    else:
        os.system('./pylooper  {input} {output}'.format(input=args.input, output=args.output))

    if args.annotate: annotate_repeats(args)

    if args.analyse:
        if args.format == 'fasta': analyse_fasta(args)
        elif args.format == 'fastq': analyse_fastq(args)

if __name__ == "__main__":
    main()
