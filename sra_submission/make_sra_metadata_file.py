"""
A script that generates the SRA submission files for our project.
"""
# -*- coding: utf-8 -*-
import pandas as pd
import os
import numpy as np

# list_of_terms = ['RNA-Seq', 'TRANSCRIPTOMIC', 'PCR', 'single', 'ILLUMINA',
#                  'Illumina', 'HiSeq 2500', 'fastq.gz']
#
column_names = ['library_strategy', 'library_source', 'library_selection',
                'library_layout', 'platform', 'instrument_model',
                'design_description', 'filetype', 'filename']

biosample_info = pd.read_csv('sra_biosamples.txt', sep='\t')
library_map = pd.read_csv('library_to_ID.txt', sep='\t', comment='#')

# biosample_info = biosample_info[['bioproject_accession', 'biosample_accession',
#                                  'library_ID', 'title']]

proj_to_file = {}
for project in os.listdir('../input/rawseq'):
    if ~(library_map.project_name == project).any():
        continue
    if project == '.DS_Store':
        continue
    for fname in os.listdir('../input/rawseq/' + project):
        if '.fastq' in fname:
            if project in proj_to_file.keys():
                proj_to_file[project] += [fname]
            else:
                proj_to_file[project] = [fname]

with open('sra_metadata.tsv', 'w') as file:
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
