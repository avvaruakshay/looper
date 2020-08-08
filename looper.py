#!/usr/bin/env python
import os, sys
import math
from datetime import datetime

"""Edits the ebeta.cpp template to compile user inputs as constants"""

def motif_factors(m, M):
    factors = []
    for i in range(M, m-1, -1):
        check = 0
        for j in factors:
            if j % i == 0: check = 1; break
        if check == 0: factors.append(i)
    
    return factors

start = datetime.now()
m = 1
M = int(sys.argv[3])
cutoff = 2*M


f = motif_factors(m, M)
n = len(f)
motifs = '{'
for a in f: motifs += str(a) + ','
motifs = motifs[:-1] + '}'

norm = str(int('1'*(2*cutoff), 2)) + 'ull'
divisor = []
rem_shift = []
for i, a in enumerate(f): 
    d = cutoff // f[i]
    r = cutoff % f[i]
    D = ('0'*((2*a)-1) + '1')*d + '0'*(2*r)
    divisor.append(int(D, 2))
    rem_shift.append( int(2*(cutoff - (cutoff % f[i]))) )

norm_string = f'uint64_t norm = {norm}'

div_string = f'uint64_t divisor[{n}] = ' + '{'
for a in divisor: div_string += str(a) + ','
div_string = div_string[:-1] + '}'

rem_shift_string = f'uint rem_shift[{n}] = ' + '{'
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
                \n{prefix}uint motif_checks[{n}] = {motifs};\
                \n{prefix}uint N = {n};\
                \n{prefix}{norm_string};\
                \n{prefix}{div_string};\
                \n{prefix}{rem_shift_string};'
        print(line, file=script)
script.close()

# print('Script generated - ', datetime.now()-start)
os.system('g++ ./pylooper.cpp -O3 -o pylooper')
# print('Script compiled - ', datetime.now()-start)
os.system(f'./pylooper  {sys.argv[1]} {sys.argv[2]}')
# print('Complete - ', datetime.now()-start)
