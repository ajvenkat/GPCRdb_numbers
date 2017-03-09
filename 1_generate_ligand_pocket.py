##### Program to create structure-based alignments of the binding pockets of class A GPCRs
from urllib2 import urlopen
import json
import re
from Bio.SeqUtils import seq1
import sys

##### Build class-wide and subgroup-based alignments of the binding pocket
def Map_gpcr_family_to_uniprot():
	family_to_uniprot_dict = {}
	url = 'http://gpcrdb.org:80/services/proteinfamily/proteins/001/'
	response = urlopen(url)
	gpcr_family_data = json.loads(response.read().decode('utf-8'))
	for entry in gpcr_family_data:
		# if entry['species'] == 'Homo sapiens':
		(family, uniprot) = (entry['family'], entry['entry_name'])
		if family in family_to_uniprot_dict.keys():
			family_to_uniprot_dict[family][uniprot] = 1
		else:
			family_to_uniprot_dict[family] = {}
			family_to_uniprot_dict[family][uniprot] = 1
	return family_to_uniprot_dict

##### Get list of class A GPCR ids from all species
#### This can be obtained from GPCRdb downloaded master file containing the GPCRdb numbers
""" function for obtaining all GPCRs (uniprot ids) in GPCRdb """
def Get_GPCRdb_uniprot_ids():
	gpcrdb_uniprot_dict = {}
	gpcrdb_generic_numbers_file = "/home/ajvenkat/gpcr-bw/All_species_gpcrdb_numbers_revised_17May2016.txt"
	with open(gpcrdb_generic_numbers_file) as GENERIC:
		for line in GENERIC:
			(uniprot, aaNum, aaName, TM, generic_num) = line.rstrip().split("\t")
			if uniprot not in gpcrdb_uniprot_dict.keys():
				gpcrdb_uniprot_dict[uniprot] = 1;
	return gpcrdb_uniprot_dict


""" function for mapping all GPCRdb numbers based on uniprot id """
def Get_GPCRdb_Numbers():
	generic_numbers_dict = {}
	gpcrdb_generic_numbers_file = "/home/ajvenkat/gpcr-bw/All_species_gpcrdb_numbers_revised_17May2016.txt"
	with open(gpcrdb_generic_numbers_file) as GENERIC:
		for line in GENERIC:
			(uniprot, aaNum, aaName, TM, generic_num) = line.rstrip().split("\t")
			generic_num = re.sub("\.\d+", "", generic_num)
			if uniprot in generic_numbers_dict.keys():
				generic_numbers_dict[uniprot][generic_num] = aaName
			else:
				generic_numbers_dict[uniprot] = {}
				generic_numbers_dict[uniprot][generic_num] = aaName
	return generic_numbers_dict

##### Get subgroup information for GPCRs from GPCRdb
generic_numbers_dict = Get_GPCRdb_Numbers()
family_to_uniprot_dict = Map_gpcr_family_to_uniprot()

"""
Binding pocket residues were obtained from the analysis of ~150 GPCR structures
"""
super_binding_pocket_combined = "1x35 1x39 2x53 2x57 2x60 2x61 2x63 2x64 23x48 23x50 3x21 3x25 3x28 3x29 3x30 3x32 3x33 3x34 3x36 3x37 3x40 3x44 3x47 3x51 4x56 4x57 4x61 4x62 45x50 45x51 45x52 5x36 5x39 5x40 5x43 5x44 5x46 5x461 5x47 5x48 5x50 5x51 5x54 5x57 6x42 6x44 6x46 6x48 6x49 6x51 6x52 6x53 6x54 6x55 6x58 6x59 6x61 6x62 7x27 7x30 7x31 7x34 7x35 7x37 7x38 7x39 7x41 7x42"
super_binding_pocket_gpcrdb_list = super_binding_pocket_combined.split(" ")

"""
For class A GPCRs (in family_to_uniprot_dict), obtain the uniprot identifier. 
"""
binding_pocket_allSpecies_outfile = "/home/ajvenkat/projects/customized-pocket/results/classA-binding-pockets-allSpecies.fasta"
binding_pocket_allSpecies_FH = open(binding_pocket_allSpecies_outfile, 'w')

binding_pocket_human_outfile = "/home/ajvenkat/projects/customized-pocket/results/classA-binding-pockets-human.fasta"
binding_pocket_human_FH = open(binding_pocket_human_outfile, 'w')

for family_id in sorted (family_to_uniprot_dict.keys()):
	for uniprot_id in sorted (family_to_uniprot_dict[family_id].keys()):
		binding_pocket_res_list = []
		uniprot_id = uniprot_id.upper()

		if uniprot_id in generic_numbers_dict.keys():
			for gpcrdb_number in super_binding_pocket_gpcrdb_list:
				if gpcrdb_number in generic_numbers_dict[uniprot_id].keys():
					aminoacid = seq1(generic_numbers_dict[uniprot_id][gpcrdb_number])
					binding_pocket_res_list.append(aminoacid)
				else:
					binding_pocket_res_list.append('-')

			binding_pocket_res_combined = ''.join(binding_pocket_res_list)

			# Writing out in fasta format
			fasta_header = family_id + ':' + uniprot_id
			print ">%s\n%s" % (fasta_header, binding_pocket_res_combined)
			binding_pocket_allSpecies_FH.write (">%s\n%s\n" % (fasta_header, binding_pocket_res_combined))

			match = re.search('HUMAN', uniprot_id)
			if match:
				print ">%s\n%s" % (fasta_header, binding_pocket_res_combined)
				binding_pocket_human_FH.write (">%s\n%s\n" % (fasta_header, binding_pocket_res_combined))


print binding_pocket_allSpecies_outfile
print binding_pocket_human_outfile