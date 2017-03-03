"""A script that generates the SRA submission files for our project."""
# -*- coding: utf-8 -*-
import pandas as pd
import argparse as arg
import os
import numpy

parser = arg.ArgumentParser(description='Generate the SRA metadata files')
parser.add_argument('biosample', type=str,
                    help='path to the biosample file')
parser.add_argument('library_map', type=str,
                    help='The library ID to folder name mapping')
parser.add_argument('dir', type=str, help='Path to the read folders')
parser.add_argument('output', type=str, help='filename to write to')

args = parser.parse_args()

biosample_info = pd.read_csv(args.biosample, sep='\t')
library_map = pd.read_csv(args.library_map, sep='\t', comment='#')
# biosample_info = pd.read_csv('sra_biosamples.txt', sep='\t')
# library_map = pd.read_csv('library_to_ID.txt', sep='\t', comment='#')

column_names = ['library_strategy', 'library_source', 'library_selection',
                'library_layout', 'platform', 'instrument_model',
                'design_description', 'filetype', 'filename']

# biosample_info =biosample_info[['bioproject_accession','biosample_accession',
#                                  'library_ID', 'title']]

proj_to_file = {}
for project in os.listdir(args.dir):
    if ~(library_map.project_name == project).any():
        continue
    if project == '.DS_Store':
        continue
    for fname in os.listdir(args.dir + project):
        if '.fastq' in fname:
            if project in proj_to_file.keys():
                proj_to_file[project] += [fname]
            else:
                proj_to_file[project] = [fname]

with open(args.output, 'w') as file:
    columns = ['bioproject_accession', 'biosample_accession', 'library_ID',
               'title'] + column_names
    columns = '\t'.join(columns)
    file.write(columns + '\n')
    for genotype in library_map.genotype.unique():
        line = ''
        ind = (library_map.genotype == genotype)
        project_ids = library_map[ind].project_name.unique()
        for project in project_ids:
            for f in proj_to_file[project]:
                line += str(f) + '\t'
        line = line[:-1]
        ind = biosample_info.library_ID == genotype
        x = biosample_info[ind].values[0].tolist()
        x = [str(i) for i in x]
        x = '\t'.join(x) + '\t' + line + '\n'
        file.write(x)
