#!usr/bin/python
from __future__ import print_function, division
from operator import itemgetter
from collections import defaultdict
from tqdm import tqdm
import sys, gzip
from os.path import splitext
from pprint import pprint

from utils import rawcharCount


"""

    CAUTION: Works currently for only sorted bed files.

    Preferential order of assigning annotation:
        Promoter >> Overlapping >> Intergenic

    Defaults:
        > Promoter distance is 1kb upstream and downstream of TSS.
        > Gene id considered is "gene".

    Built by checking on GFF3 file.

    # Sample GTF
    1 transcribed_unprocessed_pseudogene  gene        11869 14409 . + . gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; 
    1 processed_transcript                transcript  11869 14409 . + . gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_sourc e "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";

    # Sample GFF
    X	Ensembl	Repeat	2419108	2419128	42	.	.	hid=trf; hstart=1; hend=21
    X	Ensembl	Repeat	2419108	2419410	2502	-	.	hid=AluSx; hstart=1; hend=303
    X	Ensembl	Repeat	2419108	2419128	0	.	.	hid=dust; hstart=2419108; hend=2419128
"""

def select_anno(List):
    """Function to assign the hierarchically right choice of annotation"""
    if 'Exon' in List:
        return 'Exon'
    elif 'Intron' in List:
        return 'Intron'
    elif 'Genic' in List:
        return 'Genic'
    elif 'Intergenic' in List:
        return 'Intergenic'


def promoter(check):
    if check == 1:
        return 'Promoter'
    else:
        return 'Non-Promoter'


# Need to be updated for better parsing of the attributes
def process_attrs(attribute, annotype):
    """Processes the attribute field to build a dictionary with detailed
     profiling of the feature"""
    
    attr_obj = {}
    attributes = attribute.split(";")
    subdelim = " " if annotype=='GTF' else "="
    for a in attributes:
        attr = a.split(subdelim)
        attrName = attr[0].strip()
        attr_obj[attrName] = attr[1].strip()
    return attr_obj


def process_annofile(annofile, annotype, gene_id):

    """
        Processes the input annotation file and builds an object inclusive of 
        all the available genic features.
        Order of columns for a typical GFF or GTF file:
        - seqname  source  feature start   end score   strand  frame   attribute

        The output is an object is a constituent of two dictionaries:
            - An object for all the gene features. (key: chromosome_name, 
              value: list of features in the chromosome).
            - An object for all the sub gene (exon, cds, etc.) features. 
                key: chromosome_name,
                value:
                    key: parent_geneid, 
                    value: list of features in the chromosome
       
        The features for each chromosome are sorted based on their starts 
        for easier downstream processing.
    """

    gene_obj = defaultdict(list)
    subgene_obj = defaultdict()
    if (annofile.endswith('gz')):
        annohandle = gzip.open(annofile, 'rt')
    else:
        annohandle = open(annofile, 'r')
    for line in annohandle:
        line = line.strip()
        if line.startswith('#'):
            pass
        else:
            fields = line.split('\t')
            seqname, source, feature = fields[:3]
            start = int(fields[3])
            end = int(fields[4])
            score, strand, frame, attribute = fields[5:9]
            attr_obj = process_attrs(attribute, annotype)

            if feature in set(['gene', 'exon']):
                try:
                    gene_name = attr_obj[gene_id]
                except KeyError:
                    print('\nGeneKeyError:')
                    print('The attribute "%s" is not among the attributes for \
                        gene. Please select a different one.' %(gene_id))
                    print('The available ones are [' + 
                        (", ".join(list(attr_obj.keys()))) +']', end='\n\n')
                    sys.exit(1)

            if feature == 'gene':
                gene_obj[seqname].append([gene_name, start, end, strand])
            elif feature == 'exon':
                try:
                    subgene_obj[gene_name][feature].append([start, end, strand])
                except KeyError:
                    subgene_obj[gene_name] = {feature: [[start, end, strand]]}
    for i in gene_obj:
        #sorting based on the start of the feature
        gene_obj[i] = sorted(gene_obj[i], key=itemgetter(1)) 
    for a in subgene_obj:
        for b in subgene_obj[a]:
            #sorting based on the start of the feature
            subgene_obj[a][b] = sorted(subgene_obj[a][b], key=itemgetter(0)) 

    return {'gene': gene_obj, 'subgene': subgene_obj}


def annotate_repeats(args):
    """
        Main function which iterates over the given input bedfile(perf_output)
        Annotates each repeat location based on the close genic features.

        Simple outline:
            - Works with the assumption that perf_output is sorted based on co-ordinates.
            - For each repeat 
                * The features on its chromosome are retrived
                *
    """

    rep_file = args.output
    anno_file = args.annotate
    annotype = args.anno_format
    output_file = open(splitext(rep_file)[0] + '_annotation.tsv', 'w')
    gene_id = args.gene_key

    promUp = args.up_promoter
    promDown = args.down_promoter

    features = process_annofile(anno_file, annotype, gene_id)
    gene_obj = features['gene']
    subgene_obj = features['subgene']
    # pprint(features)

    print('', end='\n')
    print('Generating annotations for identified repeats..')
    print('', end='\n')
    total_regions = rawcharCount(rep_file, '\n')
    # Counting the number of lines in bed -------------------------------------
    with open(rep_file) as bed:
        prev_seqname = "Initialise" # Initialise for checking the prev_seqname
        min_start_index = 0
        for line in tqdm(bed, total=total_regions):
            # Object for the output entries to be appended --------------------
            Annotations = {'Genic': [], 'Exon': [], 'Intron': []}
            line = line.strip()
            if line.startswith('#'):
                print(line.strip(), file = output_file)
            else:
                fields = line.split('\t')
                seqname = fields[0]
                # If the seqname is not same the previous seq name the check
                # starts from the first gene on the sequence.
                if seqname != prev_seqname:
                    min_start_index = 0
                prev_seqname = seqname
                S1 = int(fields[1])     # start of the region
                E1 = int(fields[2])     # end of the region
                least_dist = float('inf')
                break_check = 0
                promoter_check = 0
                try:
                    for i, a in enumerate(gene_obj[seqname][min_start_index:]):
                        annotation = ''
                        geneName = a[0]

                        # checking for sub gene elements of the gene
                        try: subgeneElements = subgene_obj[geneName]
                        except KeyError: subgeneElements = {}
                        S2 = a[1]   # start of the feature
                        E2 = a[2]   # stop of the feature
                        Ori = a[3]  # strand of the feature
                        # Transcription Start site
                        TSS = S2 if Ori == '+' else E2
                        
                        """
                            The below conditions make an optimal choice of a 
                            feature from where distance relation comparisons for
                             this paricular repeat can be initiated.
                            
                            Basic assumption-

                            S2-------Feature-------E2 |--- > maximum promoter distance---|                      |--- > maximum promoter distance---| S2-------Feature-------E2
                                                      |--- > maximum promoter distance---| S1-----Repeat-----E1 |--- > maximum promoter distance---|
                            
                            With the condition of choosing the closest feature 
                            which is at least at a distance length of the 
                            promoter, we can omit comparisons with features 
                            which are much farther away.

                            The point where the comparisons stop is the closest 
                            feature which is greater than the distance of 
                            promoter from the end of the repeat.
                        """
                        if i == 0: # if it's the first gene on the chromosome
                            least_start = S1 - E2
                            min_index = i
                        else:
                            # if the region goes beyond the promoter distance
                            # from the end of the feature change the min_index
                            if S1 - E2 > max([promUp, promDown]):
                                if least_start > (S1 - E2):
                                    least_start = S1 - E2
                                    min_index = i
                        if break_check == 1: break
                        #point to break the comparisons done
                        if (S2 - E1 > max([promUp, promDown])): 
                            break_check = 1
                        
                        
                        # Checking if region comes in promoter 
                        # For positive strand orientation 
                        if Ori == '-' and (TSS-promDown <= S1 <= TSS+promUp or TSS-promDown <= E1 <= TSS+promUp):
                            promoter_check = 1
                        elif Ori == '+' and (TSS-promUp <= S1 <= TSS+promDown or TSS-promUp <= E1 <= TSS+promDown):
                            promoter_check = 1
                        # If no Promoter found 
                        # Checking if it overlaps 
                        if (E2 - S1 >=0 and S2 - S1 <=0) or (E2 - E1 >=0 and S2 - E1 <= 0):
                            annotation = 'Genic'
                            # Removes the Intergenic entries cause a Genic overlap is found 
                            Annotations['Intergenic'] = []
                            TSS = S2
                            diffSS = S2 - S1
                            diffES = E2 - S1
                            diffSE = S2 - E1
                            diffEE = E2 - E1
                            if abs(diffSS) < abs(diffSE): tss_dist = diffSS
                            else: tss_dist = diffSE
                            distance = tss_dist
                            # Checking overlap with subgene elements
                            if 'exon' in subgeneElements:
                                for site in subgeneElements['exon']:
                                    S3 = site[0]
                                    E3 = site[1]
                                    if (E3 - S1 >=0 and S3 - S1 <=0) or (E3 - E1 >=0 and S3 - E1 <= 0):
                                        annotation = "Exon"
                                        break
                                    else:
                                        annotation = "Intron"
                        elif len(Annotations['Exon']) == 0 and len(Annotations['Intron']) == 0 and len(Annotations['Genic']) == 0:
                            # If no Genic annotations are found closest distace from the closest gene is calculated
                            TSS = S2
                            diffSS = S2 - S1
                            diffES = E2 - S1
                            diffSE = S2 - E1
                            diffEE = E2 - E1
                            if abs(diffSS) < abs(diffSE):
                                tss_dist = diffSS
                            else:
                                tss_dist = diffSE
                            min_dist = min([abs(diffSS), abs(diffEE), abs(diffSE), abs(diffES)])
                            annotation = 'Intergenic'
                            distance = tss_dist
                            if min_dist < least_dist:
                                least_dist = min_dist
                                Annotations[annotation] = [line + '\t' + '\t'.join(str(b) for b in a) + '\t' + annotation + '\t' + promoter(promoter_check) + '\t' + str(distance)]

                        if (annotation == "Exon" or annotation == "Intron" or annotation == "Genic"):
                            Annotations[annotation].append(line + '\t' + '\t'.join(str(b) for b in a) + '\t' + annotation + '\t' + promoter(promoter_check) + '\t' + str(distance))

                    min_start_index += min_index
                    # print(min_start_index)

                    if min_start_index > 0:
                        # Cautious assignment the closest genic feature to start comparisons from
                        min_start_index -= 1
                # If sequence is not found, reports as annotation not available
                except KeyError:
                    Annotations = {'Genic': [], 'Exon': [], 'Intron': []}
                    print(line + '\t' + '\t'.join(['-']*7), file = output_file)
                for anno in list(Annotations.keys()):
                    if len(Annotations[anno]) == 0: del Annotations[anno]
                for anno in Annotations:
                    feature_leastdist = float('inf')
                    closest_entry = ""
                    for entry in Annotations[anno]:
                        feature_dist = int(entry.split('\t')[-1])
                        if feature_dist < feature_leastdist:
                            feature_leastdist = feature_dist
                            closest_entry = entry
                    if closest_entry != "":
                        Annotations[anno] = closest_entry
                if len(Annotations) > 1:
                    anno_selected = select_anno(list(Annotations.keys()))
                    print(Annotations[anno_selected], file = output_file)
                else:
                    for anno in Annotations:
                        print(Annotations[anno], file = output_file)
    output_file.close()

if __name__ == "__main__":
    anno_file = sys.argv[1]
    anno_format = sys.argv[2]
    gene_id = 'gene' if (anno_format == 'GFF') else 'gene_id'
    annotation_obj = process_annofile(anno_file, anno_format, gene_id)
    pprint(annotation_obj['subgene'])