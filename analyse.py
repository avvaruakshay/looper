#! /usr/bin/env python

from __future__ import print_function, division
import sys, os, json
from collections import Counter, defaultdict
import numpy as np
from pprint import pprint

from utils import kmers, get_cycles, build_cycVariations, rev_comp

def generate_defaultInfo(args):
    repeats_file = args.output
    input_file = args.input

    html_report = os.path.splitext(repeats_file)[0] + '.html'
    print("\nGenerating HTML report. This may take a while..", end="\n\n")

    inf = float('inf')

    repeat_classes = []
    kmer_classes = defaultdict(list)
    cyclical_variations = dict()
    total_repeat_bases = 0
    total_repeat_freq = 0
    longest_lengths = [['seq', 'start', 'stop', 'repeat_class', 0, '+', 0, 'actualrep']]*100
    most_units = [['seq', 'start', 'stop', 'repeat_class', 0, '+', 0, 'actualrep']]*100
    min_length = inf
    min_units = inf

    plot_data = {}
    defaultInfo = {}
    defaultInfo['info'] = { 'SeqInfo': {}, 'RepInfo': {} }

    with open(repeats_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#'):
                fields = line[1:].split(': ')
                defaultInfo['info']['SeqInfo'][fields[0]] = fields[1]
            else:
                fields = line.split('\t')
                fields = line.split('\t')
                fields[1] = int(fields[1])
                fields[2] = int(fields[2])
                fields[4] = int(fields[4])
                fields[6] = int(fields[6])

                seq_name = fields[0]
                repeat_start = fields[1]
                repeat_end = fields[2]
                repeat_class = fields[3]
                repeat_length = fields[4]
                repeat_ori = fields[5]
                repeat_units = fields[6]
                repeat_actual = fields[7]

                if repeat_class not in repeat_classes:
                    repeat_classes.append(repeat_class)
                    cyclical_variations[repeat_class] = build_cycVariations(repeat_class)

                    plot_data[repeat_class] = dict()
                    plot_data[repeat_class][repeat_length] = [0]*len(cyclical_variations[repeat_class])

                total_repeat_bases += repeat_length
                total_repeat_freq += 1

                if repeat_length not in plot_data[repeat_class]: 
                    plot_data[repeat_class][repeat_length] = \
                        [0]*len(cyclical_variations[repeat_class])

                plot_data[repeat_class][repeat_length] \
                    [cyclical_variations[repeat_class].index(repeat_actual)] += 1

                if min_units > repeat_units: min_units = repeat_units
                if min_length > repeat_length: min_length = repeat_length

                if (longest_lengths[-1][4] < repeat_length) or \
                (longest_lengths[-1][4] == repeat_length and \
                repeat_class < longest_lengths[-1][3]):
                    longest_lengths[-1] = fields
                    longest_lengths.sort(key=lambda x: x[4])
                    longest_lengths.reverse()
                
                if (most_units[-1][6] < repeat_units) or \
                (most_units[-1][6] == repeat_units and \
                repeat_class < longest_lengths[-1][3]):
                    most_units[-1] = fields
                    most_units.sort(key=lambda x: x[6])
                    most_units.reverse()
    
    total_bases = int(defaultInfo['info']['SeqInfo']['GenomeSize'])
    defaultInfo['info']['RepInfo']['PlotData'] = plot_data
    defaultInfo['info']['RepInfo']['NumRepClasses'] = len(plot_data.keys())
    defaultInfo['info']['RepInfo']['TotalRepBases'] = total_repeat_bases
    defaultInfo['info']['RepInfo']['TotalRepFreq'] = total_repeat_freq
    defaultInfo['info']['RepInfo']['PercentGenomeCovered'] = \
      str(round((total_repeat_bases/total_bases)*100, 2)) + "%"
    defaultInfo['info']['RepInfo']['RepDensityByFreq'] = \
      round((total_repeat_freq/total_bases)*1000000, 2)
    defaultInfo['info']['RepInfo']['RepDensityByBases'] = \
      round((total_repeat_bases/total_bases)*1000000, 2)
    defaultInfo['info']['RepInfo']['minLength'] = min_length
    defaultInfo['info']['RepInfo']['minUnits'] = min_units
    defaultInfo['info']['RepInfo']['LongestRepeats'] = []
    defaultInfo['info']['RepInfo']['MostRepeatUnits'] = []
    for a in longest_lengths:
        testDict = {
            'Seq': a[0], 'Start': a[1], 'End': a[2], 'RepClass': a[3], 
            'RepLength': a[4], 'RepOri': a[5], 'RepUnit': a[6], 'ActualRep': a[7]
        }
        defaultInfo['info']['RepInfo']['LongestRepeats'].append(testDict)
    for a in most_units:
        testDict = {
            'Seq': a[0], 'Start': a[1], 'End': a[2], 'RepClass': a[3], 
            'RepLength': a[4], 'RepOri': a[5], 'RepUnit': a[6], 'ActualRep': a[7]
        }
        defaultInfo['info']['RepInfo']['MostRepeatUnits'].append(testDict)
    defaultInfo = 'const data =' + json.dumps(defaultInfo)
    print(defaultInfo)
    # writetoHTML(html_report, defaultInfo, repeat_options, 'fasta')


def writetoHTML(html_file, defaultInfo, repeat_options, input_format):
    html_handle = open(html_file, 'w')
    current_dir = os.path.dirname(__file__)

    template = open('%s/lib/template_%s.html' %(current_dir, input_format), 'r').read()

    fontawesome_js = open('%s/lib/src/all.js' %(current_dir), 'r').read()
    semantic_css = open('%s/lib/styles/semantic.min.css' %(current_dir), 'r').read()
    multiselect_css = open('%s/lib/styles/multi-select.min.css' %(current_dir), 'r').read()
    apexcharts_css = open('%s/lib/styles/apexcharts.min.css' %(current_dir), 'r').read()
    main_css = open('%s/lib/styles/main.css' %(current_dir), 'r').read()

    jquery_js = open("%s/lib/src/jquery-3.5.0.min.js" %(current_dir), "r").read()
    semantic_js = open("%s/lib/src/semantic.min.js" %(current_dir), "r").read()
    multiselect_js = open('%s/lib/src/jquery.multi-select.min.js' %(current_dir), 'r').read()
    apexcharts_js = open('%s/lib/src/apexcharts.min.js' %(current_dir), 'r').read()
    lodash_js = open('%s/lib/src/lodash.min.js' %(current_dir), 'r').read()
    main_js = open('%s/lib/src/main_%s.js' %(current_dir, input_format), 'r').read()
    tables_js = open('%s/lib/src/tables_%s.js' %(current_dir, input_format), 'r').read()
    annocharts_js = ''
    if input_format == 'fasta':
        annocharts_js = open('%s/lib/src/anno_charts.js' %(current_dir), 'r').read()

    template = template.format(
        fontawesome_js = fontawesome_js, 
        semantic_css = semantic_css, 
        multiselect_css = multiselect_css, 
        apexcharts_css = apexcharts_css, 
        main_css = main_css, 
        jquery_js = jquery_js, 
        semantic_js = semantic_js, 
        multiselect_js = multiselect_js, 
        apexcharts_js = apexcharts_js, 
        lodash_js = lodash_js, 
        analyse_data_js = defaultInfo, 
        main_js = main_js, 
        tables_js = tables_js, 
        annocharts_js = annocharts_js,
        repeat_options = repeat_options,
    )

    print(template, file=html_handle)
    html_handle.close()
    print("HTML report successfully saved to " + html_file)


def get_parameters(args):
    runCommand = 'PERF' + ' '.join(sys.argv)

