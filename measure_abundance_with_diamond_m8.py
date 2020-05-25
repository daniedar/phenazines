import sys,os,re
import subprocess
import numpy as np
import pandas as pd
import configparser
import copy
import map_and_count_functions


##############################################################################################################################################
### Description							
# Extract total phz+ and phzdA/podA+ abundance (estimted % of all bacteria in sample)
# - The script takes a diamond result file in .m8 format (see example file and reference diamond fasta and databse)
# - A configuration file can be used to change parameters (e.g., mapping identity threholds) and to adjust paths (see example run_conf.ini).
# - There are several files used by the script:
#	- "all_genes_with_categories_table.txt" = gene groups defined into phenazine, biodegradation or reference (total-bacteria)
#	- "sra_sample_table.txt" = a list of sra samples with their annotation for the final output table
#	- "taxonomy_db.txt" = full taxonomic data for all 181 phenazine producer genomes (

### How to run
#   python run_pipeline_local.py [configuration file] [Section]
# 	 - ex: python measure_abundance_with_diamond_m8.py run_conf.ini example_run
#
#	Script makes a new dir with 3 result files:
#	 - .diamond.m8.normalized_read_counts.txt: gene_name, total read counts, average gene size (kb), RPK normalized counts
#	 - ._abundances_results.txt = abundance results and relative abundance of different producer clades + normalized raw counts
#	 - .relative_genome_id_for_tax.txt = the relative abundance of each genome out of all phenazine producers in sample
#
##############################################################################################################################################

### Collect run parameters from conf file
if len(sys.argv) < 3:
	print('missing conf file and/or section...\n')
	exit(0)

config = configparser.ConfigParser()
config.read(sys.argv[1])
section = sys.argv[2]

min_pindet = float(config.get('parameters','min_pindet'))
reference_gene_db_path = config.get('parameters','reference_gene_db_path')
output_dir = config.get('parameters','output_dir')
taxonomic_output_level = config.get('parameters','taxonomic_output_level') # report the results according to this taxonomic level. e.g., "Family"
ref_sra_table = config.get('parameters','ref_sra_table')

#make result dir
map_and_count_output_dir = '%s/%s/'%(output_dir,section)
subprocess.Popen('mkdir -p %s/%s'%(output_dir,section), stdout=subprocess.PIPE, shell=True)


# collect all sra accessions and their annotations 
sra_ref_db = {}
with open(ref_sra_table,'r') as f:
	f.readline()
	for line in f:
		run_acc,sample_type,sample_subtype,organism,gold = line.rstrip().split('\t')
		sra_ref_db[run_acc] = {
			'run_acc':run_acc,
			'sample_type':sample_type,
			'sample_subtype':sample_subtype,
			'organism':organism,
			'gold':gold,
		}

#get genes and their categories for output stage: e.g., reference, phenazines, degradation
gene_categories = {}
with open(config.get('parameters','genes_and_categories')) as f:
	f.readline()
	for line in f:
		gene_name,category = line.rstrip().split('\t')
		gene_categories[gene_name] = category

### taxonomy db preparation
org_taxonomy_db = {}
class_names = {}
tax_ref = {}
tax_to_class = {}
with open(config.get('parameters','taxonomy_db')) as f:
	f.readline() 
	for line in f: #IMG genome ID	Phylum	Class	Order	Family	Genus	Species
		
		genome_id,Phylum,Class,Order,Family,Genus,Species = line.rstrip().split('\t')
		org_taxonomy_db[genome_id] = {
			'Phylum':Phylum,
			'Class':Class,
			'Order':Order,
			'Family':Family,
			'Genus':Genus,
			'Species':Species
		}
		for tax_level,specific_tax in org_taxonomy_db[genome_id].items():
			if tax_level == taxonomic_output_level:
				tax_to_class[specific_tax] = org_taxonomy_db[genome_id]['Class']
			if not tax_level in tax_ref:
				tax_ref[tax_level] = {
					specific_tax : 0
				}
			else:
				tax_ref[tax_level][specific_tax] = 0
		
		
sorted_by_class_tax_values = [] 
for x in sorted(tax_to_class.items(), key=lambda kv: kv[1]):
	sorted_by_class_tax_values.append(x[0])


#get files from result dir (m8, etc)
path_to_results_dir = config.get(section,'path_to_sra_results')
all_files_in_results_dir = os.listdir(path_to_results_dir)


#run over each sample
for sample in all_files_in_results_dir:
	if re.search('\.m8',sample):
		m = re.search('(\w+)\.diamond.m8',sample)
		run_acc = m.group(1)
		
		file_path = '%s/%s'%(path_to_results_dir,sample)
		print(file_path)
		
		run_db = {}
		run_taxonomy_db = {}
		tax_freq_db = copy.deepcopy(tax_ref)

		
		reference_gene_count_db,genome_id_abundance = map_and_count_functions.map_and_count_m8(sample,file_path,map_and_count_output_dir,reference_gene_db_path,min_pindet)
		
		# collect data for analysis 
		run_data_and_results = {
			'run_acc' : run_acc,
			'sample_type':sra_ref_db[run_acc]['sample_type'],
			'sample_subtype':sra_ref_db[run_acc]['sample_subtype'],
			'organism':sra_ref_db[run_acc]['organism'],
			'gold':sra_ref_db[run_acc]['gold'],
			'ref_median_value' : 'na',
			'phz_median_value' : 'na',
			'median_based_freq_estimate' : 'na'
		}
		
		#prep gene count dbs per group
		ref_genes_counts = {}
		phz_genes_counts = {}
		deg_genes_counts = {}
		for gene_name,category in gene_categories.items():
			if category == 'reference': ref_genes_counts[gene_name] = 0
			if category == 'phenazines': phz_genes_counts[gene_name] = 0
			if category == 'degradation': deg_genes_counts[gene_name] = 0
		
		#read results 
		for gene_name in reference_gene_count_db:
			gene_rpk = reference_gene_count_db[gene_name]['rpk']
			if gene_name in ref_genes_counts:
				ref_genes_counts[gene_name] += float(gene_rpk)
			elif gene_name in phz_genes_counts:
				phz_genes_counts[gene_name] += float(gene_rpk)
			elif gene_name in deg_genes_counts:
				deg_genes_counts[gene_name] += float(gene_rpk)
			else:
				print(gene_name)
		
		phz_genes_counts['phzAB'] = phz_genes_counts['phzA'] + phz_genes_counts['phzB'] #merge phzA and phzB as they are very similar sequence wise
		total_phz_signal = sum([phz_genes_counts['phzAB'], phz_genes_counts['phzD'],phz_genes_counts['phzE'],phz_genes_counts['phzF'],phz_genes_counts['phzG']])
		genome_id_relative_abundance = {}
		for gene_name in ['phzA','phzB','phzD','phzE','phzF','phzG']:
			for genome_id in genome_id_abundance[gene_name]:
				if not genome_id in genome_id_relative_abundance:
					genome_id_relative_abundance[genome_id] = genome_id_abundance[gene_name][genome_id]
				else:
					genome_id_relative_abundance[genome_id] += genome_id_abundance[gene_name][genome_id]
	
		### calculate the median scores
		ref_gene_median = np.median( list(ref_genes_counts.values()) ) 
		phz_gene_median = np.median( [phz_genes_counts['phzAB'], phz_genes_counts['phzD'],phz_genes_counts['phzE'],phz_genes_counts['phzF'],phz_genes_counts['phzG']] )
		run_data_and_results['ref_median_value'] = round(ref_gene_median,2)
		run_data_and_results['phz_median_value'] = round(phz_gene_median,2)
		
		### calculate the abundance estimate for phz+ bacteria
		median_based_freq_estimate = (phz_gene_median / ref_gene_median) * 100
		run_data_and_results['median_based_freq_estimate'] = round(median_based_freq_estimate,2)		
		
		podA_fraction = deg_genes_counts['podA'] / run_data_and_results['ref_median_value'] * 100 
		phdA_fraction = deg_genes_counts['phdA'] / run_data_and_results['ref_median_value'] * 100
	
		run_data_and_results['podA_freq'] = podA_fraction #single gene
		run_data_and_results['phdA_freq'] = phdA_fraction #single gene
		
		run_data_and_results = {**run_data_and_results, **ref_genes_counts}
		run_data_and_results = {**run_data_and_results, **phz_genes_counts}
		run_data_and_results = {**run_data_and_results, **deg_genes_counts}
		
		
		#get taxonomic data
		for genome_id in genome_id_relative_abundance:
			if total_phz_signal > 0:
				genome_id_relative_abundance[genome_id] = genome_id_relative_abundance[genome_id] / float(total_phz_signal) * 100
				for tax_level,specific_tax in org_taxonomy_db[genome_id].items():
					tax_freq_db[tax_level][specific_tax] += genome_id_relative_abundance[genome_id]	
		
		run_data_and_results = {**run_data_and_results, **tax_freq_db[taxonomic_output_level]}
		
		run_db[sample] = run_data_and_results
		run_taxonomy_db[sample] = genome_id_relative_abundance


		# Write output files 
		output_list = []
		for gene_name in phz_genes_counts: output_list.append(gene_name)
		for gene_name in deg_genes_counts: output_list.append(gene_name)
		for gene_name in ref_genes_counts: output_list.append(gene_name)
	
		column_out_order = ['run_acc','sample_type','organism','gold','sample_subtype','median_based_freq_estimate','phz_median_value','ref_median_value','podA_freq','phdA_freq'] +  sorted_by_class_tax_values + output_list #list(phz_genes_counts.keys())
		
		run_db_df = pd.DataFrame(run_db).transpose()
		run_db_df.to_csv('%s/%s/%s_abundances_results.txt'%(output_dir,section,run_acc),sep='\t',index=False,columns=column_out_order)
		
		run_taxonomy_db_df = pd.DataFrame(run_taxonomy_db)
		run_taxonomy_db_df.index.name = "genome_id"
		run_taxonomy_db_df.to_csv('%s/%s/%s_relative_genome_id_for_tax.txt'%(output_dir,section,run_acc),sep='\t')
