#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  29 20:05:00 2020

@author: catherinecho
"""
import sys
import json
import pandas as pd
import io
import os
import requests
import json

new_data = []
new_line =''
inFile = sys.argv[1]
with open(inFile,'r') as f:
    fhir = json.load(f)

#Extrating dictionary with key 'contained'
main = fhir['contained']

def api_SPDI_generator(b37):
    '''
    Parameters
    ----------
    b37 : b37 contextual spdi from JSON files

    Returns
    -------
    appropriate format of b37 to request API on NCBI

    '''  
    splitted_b37 = b37.split(":")
    combined_b37 = "%3A".join(splitted_b37)
    return combined_b37

def indel_api_request(translated_b37):
    ''' 
    Parameters
    ----------
    b37 : b37 contextual spdi from JSON files

    Returns
    -------
    appropriate format of b37 to request API on NCBI

    ''' 
    url = "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/"+ translated_b37 + "/canonical_representative"
    response = requests.get(url)
    data = response.text
    data_dict = json.loads(data)
    for dict_val in data_dict.values():
        api_generated_spdi = dict_val['seq_id']+":"+ str(dict_val['position'])+":"+ dict_val['deleted_sequence']+":"+ dict_val['inserted_sequence']
    return api_generated_spdi
    return api_generated_spdi



def parse_dict(lis, comp):
    '''
    Parameters
    ----------
    lis : list of dictionary.
    comp : key of dictionary like to extract.

    Returns
    -------
    comp_lis : list of dictionary with key.

    '''    
    comp_lis = []
    for item in lis:
        for i in item.keys():
            if i in comp:
                comp_lis.append(item[i])
    return comp_lis

comp_lis = parse_dict(main, 'component')

comp = ['valueCodeableConcept', 'valueRange', 'valueString']
new_comp = []
for i in comp_lis:
    value = parse_dict(i, comp)
    new_comp.append(value)

# Parsing reference sequence ref_seq
ref_seq = ''
for ls in new_comp[0]:
    if 'coding' in ls.keys():
        for lis in ls['coding']: 
            if 'code' in lis.keys():
                ref_seq = lis['code']
    
# Creating position array 
pos = [] 
for lis in new_comp:
    for dt in lis:
        if type(dt) != str:
            for key in dt.keys():
                if key == 'low':
                    pos.append(dt[key]['value'])
# First one is just overall info, need to pop
pos.pop(0)

# Creating dictionary of REF and ALT allele for each position
pos_dict = {}
for i in range (1, len(new_comp)):
    allele = []
    for lis in new_comp[i]:
        if type(lis) == str:
            allele.append(lis)
            pos_dict[pos[i-1]] = allele
            
# Creating dictionary of zygosity for each position
pos_zygosity = {}
for i in range (1, len(new_comp)):
    for ls in new_comp[i]:
        if type(ls) == dict:
            if 'coding' in ls.keys():
                for lis in ls['coding']: 
                    if 'display' in lis.keys():
                        if lis['display'] == 'homozygous':
                            pos_zygosity[pos[i-1]] = '1/1:.'
                        elif lis['display'] == 'heterozygous':
                            pos_zygosity[pos[i-1]] = '0/1:.'

# Creating dataframe for output
column_names = ["#CHROM", "POS", "ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE", "ZYGOSITY"]
fhir_parsing = pd.DataFrame(columns = column_names)
# Pre-filled columns with fixed values
fhir_parsing['POS'] = pos
fhir_parsing['#CHROM'] = 10
fhir_parsing['FORMAT'] = 'GT'
# Take out the value from dictionary "pos_zygosity" to fill Column "ZYGOSITY"
zygosity_list = pos_zygosity.values()
fhir_parsing['ZYGOSITY'] = zygosity_list
# Fill Columns "REF" and "ALT" variant allele
for index, row in fhir_parsing.iterrows():
    for key in pos_dict.keys():
        if row['POS'] == key:
            fhir_parsing.at[index,'REF'] = pos_dict[key][0]
            fhir_parsing.at[index,'ALT'] = pos_dict[key][1]

# Create b37spdi
for index, row in fhir_parsing.iterrows():
    row_ls = [ref_seq, str(row['POS']-1),row['REF'], row['ALT']]
    fhir_parsing.at[index, 'b37spdi'] = ':'.join(row_ls)
    
# Filll the rest of Null value with "."
fhir_parsing.fillna('.')

# Import PharmCAT conversion table
pharmcat = pd.read_csv('PharmCAT Conversion Table - PharmCAT 0.7 Conversions.csv')

#Parsing b37 info from patient varients file
b37 = fhir_parsing[["b37spdi","ZYGOSITY", "REF","ALT"]].to_numpy()

#Parsing b37 and b38 info from Pharmcat conversion table
pharmcatspdi = pharmcat[['B37SPDI','B38SPDI','B38CHROM','B38POS','B38REF','B38ALT']].to_numpy()

#Creating dictionary of matching variants from patients varient file with Pharmcat conversion table
variant_lookup = dict()
for spdi in b37:
    if (len(spdi[2]) == 1) and (len(spdi[3]) == 1): ## checking length of REF and ALT (Single Nucleotide variation)
        for pharmspdi in pharmcatspdi:
            if spdi[0] == pharmspdi[0]:
                zygosity = spdi[1]
                variant_lookup[spdi[0]] = [pharmspdi[2],pharmspdi[3],pharmspdi[4],pharmspdi[5], zygosity]
    else:
        b37_spdi = api_SPDI_generator(spdi[0])
        api_gen_spdi = indel_api_request(b37_spdi)
        for pharmspdi in pharmcatspdi:
            if api_gen_spdi == pharmspdi[1]:
                zygosity = spdi[1]
                variant_lookup[api_gen_spdi] = [pharmspdi[2],pharmspdi[3],pharmspdi[4],pharmspdi[5], zygosity]

#Convert vcf file to pandas dataframe 
def read_vcf(path):
    '''

    Parameters
    ----------
    path : name of vcf file to convert.

    Returns
    -------
    pandas dataframe.

    '''
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, 'PharmCAT': str},
        sep='\t'
    )

#Calling read_vcf
    
pharmcat_vcf_template = read_vcf('pharmcat.v0.7.0.template.vcf')


#Searching matching varients and update zygosity for the matched varients
for lis in variant_lookup.values():
    for index, row in pharmcat_vcf_template.iterrows():
        if(lis[0] == row['#CHROM']) and (lis[1] == row['POS']) and (lis [2] == row ['REF']):
            pharmcat_vcf_template.at[index,'ALT'] = lis[3]
            pharmcat_vcf_template.at[index,'Sample'] = lis[4]
        

#Metadata for vcf file
header_vcf = "##fileformat=VCFv4.1\n##fileDate=2015-08-04\n##source=Electronic, version: hg38_2.0.1\n##reference=hg38\n##INFO=<ID=PX,Number=.,Type=String,Description=\"PGX\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FILTER=<ID=PASS,Description=\"All filters passed\"> \n"

#Convert dataframe to vcf file again
pharmcat_vcf_template.to_csv(inFile +".vcf", sep ='\t' ,index = False)

#Re-open pharmacat vcf file and write metadata and export as final vcf file
with open(inFile +".vcf", 'w') as f:
    f.write(header_vcf)
    pharmcat_vcf_template.to_csv(f, sep = '\t', index = False)


    