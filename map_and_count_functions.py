import sys,os,re
import subprocess
import configparser
import numpy as np


def main():
	if __name__ == "__main__":
		main()


def map_and_count_m8(name,diamond_output_file_path,output_dir,reference_gene_db_path,min_pindet):
	
	#go_over_reference_set_and_collect genes
	reference_gene_count_db = {}
	gene_full_name_to_length = {}
	genome_id_abundance ={}
	
	# read ref genes 
	with open(reference_gene_db_path,'r') as f:
		while True:
			h = f.readline()
			s = f.readline()
			h = re.sub('>','',h)
			if not s: break
			gene_name = h.split('_')[0]#re.match('^([^_]+)_([^_]+),([^_]+)',h).group(1)
			if not gene_name in reference_gene_count_db:
				reference_gene_count_db[gene_name] = {
					'count' : 0,
					'gene_lengths' : [],
					'genome_ids': []
					}
			gene_full_name_to_length[h.rstrip()] = len(s) * 3 / 1000 # length to kb 
			
			#for relative abundance per taxon id consider only phenazine genes
			if re.search('phz',gene_name):
				genome_id = h.split('_')[2]
				if not gene_name in genome_id_abundance:
						genome_id_abundance[gene_name] = {
							genome_id : 0
						}
				else:
					genome_id_abundance[gene_name][genome_id] = 0
	
	#collect_data and count reads
	with open(diamond_output_file_path,'r') as diamond_output:
		for line in diamond_output:
			query,hit,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line.rstrip().split('\t')

			hit_sp = hit.split('_')
			hit_gene_name = hit_sp[0]
			hit_gene_id = hit_sp[1]
			hit_genome_id = hit_sp[2]

			# if reads pass threshold then count
			if float(pident) >= min_pindet:
				hit_gene_length = gene_full_name_to_length[hit]
				reference_gene_count_db[hit_gene_name]['count'] += 1
				reference_gene_count_db[hit_gene_name]['gene_lengths'].append( hit_gene_length )
				
				reference_gene_count_db[hit_gene_name]['genome_ids'].append(hit_genome_id)
				
				if hit_gene_name in genome_id_abundance:
					genome_id_abundance[hit_gene_name][hit_genome_id] += 1		
			
	
	### normalize and output data
	normalized_read_count_output_file_path = '%s/%s.normalized_read_counts.txt'%(output_dir,name)
	with open(normalized_read_count_output_file_path,'w') as output:
		for gene_name in reference_gene_count_db:
			gene_read_count = reference_gene_count_db[gene_name]['count']
			
			gene_median_length_kb,gene_rpk = 0,0
			if gene_read_count > 0:
				gene_median_length_kb = np.median( reference_gene_count_db[gene_name]['gene_lengths'] )#median gene length to be used for normalization per kb
				gene_rpk = gene_read_count / gene_median_length_kb #correct counts
			
				#normalized abundance per genome_id (for taxonomic calculations)
				if gene_name in genome_id_abundance:
					for hit_genome_id in genome_id_abundance[gene_name]:
						genome_norm_abundance = genome_id_abundance[gene_name][hit_genome_id] / gene_median_length_kb			
						genome_id_abundance[gene_name][hit_genome_id] = genome_norm_abundance
			
			reference_gene_count_db[gene_name]['rpk'] = gene_rpk
			output.write( '%s\t%s\t%s\t%s\n'%(gene_name,gene_read_count,gene_median_length_kb,gene_rpk) )
	
	return reference_gene_count_db,genome_id_abundance

	
