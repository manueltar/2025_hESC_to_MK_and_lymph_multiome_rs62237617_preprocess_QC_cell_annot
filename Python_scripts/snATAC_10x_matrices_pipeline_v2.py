#!/usr/bin/env python3

import os
import sys
import gzip
import argparse
import logging
import time
import subprocess
import pysam
import numpy as np
import pandas as pd
import scipy.sparse
from multiprocessing import Pool

def convert_10X_bam(args):
	logging.info(f"convert_10X_bam started at {time.strftime('%H:%M:%S')}")
	bf = pysam.AlignmentFile(args.input_bam, 'rb')
	cf = pysam.AlignmentFile(args.output_prefix + '.compiled.filt.bam', 'wb', template=bf)

	for read in bf:
		try:
                    barcode = read.get_tag('CB').split(':')[0]
		except:
			pass
			continue
		read.query_name = barcode + '_' + read.query_name
		cf.write(read)
	bf.close()
	cf.close()
	logging.info(f"convert_10X_bam finished at {time.strftime('%H:%M:%S')}")
	return

def remove_duplicate_reads(args):
	logging.info(f"def remove_duplicate_reads started at {time.strftime('%H:%M:%S')}")
	# Renamed for clarity as input to this function to avoid confusion with filt_temp_bam
	filt_bam_input = args.output_prefix + '.compiled.filt.bam'
	markdup_bam = args.output_prefix + '.filt.md.bam'
	rmdup_bam = args.output_prefix + '.filt.rmdup.bam'

	# Define intermediate temporary BAM names
	# These will be created sequentially and deleted after their next step
	filt_temp_bam = args.output_prefix + '.filtered.temp.bam'
	sortname_temp_bam = args.output_prefix + '.sortname.temp.bam'
	fixmate_temp_bam = args.output_prefix + '.fixmate.temp.bam'

	# List to keep track of temporary files for robust cleanup in case of error
	temp_files_to_clean = []

	try:
		# Step 1: Filter BAM reads (quality, unmapped, secondary, supplementary)
		# Output: filt_temp_bam
		filt_cmd = [
			'samtools', 'view', '-bu', '-h', '-q', str(args.map_quality),
			'-F', '256', # Skip secondary alignments
			'-F', '512', # Skip supplementary alignments
			'-F', '2048', # Skip unmapped reads
			filt_bam_input, # Input BAM
			'-o', filt_temp_bam # Output to file
		]
		logging.info(f"running filtering command started at {time.strftime('%H:%M:%S')}")
		result = subprocess.run(filt_cmd, capture_output=True, text=True, check=False)
		if result.returncode != 0:
			logging.error(f"filt_cmd failed with code {result.returncode}. Stderr: {result.stderr.strip()}")
			raise RuntimeError(f"filt_cmd failed: {result.stderr.strip()}")
		temp_files_to_clean.append(filt_temp_bam) # Add to cleanup list

		# Step 2: Sort by read name (required for samtools fixmate)
		# Input: filt_temp_bam, Output: sortname_temp_bam
		sortname_cmd = [
			'samtools', 'sort', '-n', # Sort by read name
			'-m', str(args.memory)+'G', # Max memory per thread for sorting
			'-@', str(args.threads), # Number of threads for sorting
			filt_temp_bam, # Input BAM
			'-o', sortname_temp_bam # Output to file
		]
		logging.info(f"running sortname command started at {time.strftime('%H:%M:%S')}")
		result = subprocess.run(sortname_cmd, capture_output=True, text=True, check=False)
		if result.returncode != 0:
			logging.error(f"sortname_cmd failed with code {result.returncode}. Stderr: {result.stderr.strip()}")
			raise RuntimeError(f"sortname_cmd failed: {result.stderr.strip()}")
		temp_files_to_clean.append(sortname_temp_bam)
		
		# Clean up the previous intermediate file to free disk space immediately
		if os.path.exists(filt_temp_bam):
			os.remove(filt_temp_bam)
			logging.info(f"Removed intermediate file: {filt_temp_bam}")
			# Remove from the list of files to clean up in finally block if it's already gone
			if filt_temp_bam in temp_files_to_clean:
				temp_files_to_clean.remove(filt_temp_bam)

		# Step 3: Fixmate (corrects mate-pair information)
		# Input: sortname_temp_bam, Output: fixmate_temp_bam
		fixmate_cmd = [
			'samtools', 'fixmate', '-r', # Add mate score tags
			sortname_temp_bam, # Input BAM
			fixmate_temp_bam # Output BAM
		]
		logging.info(f"running fixmate command started at {time.strftime('%H:%M:%S')}")
		result = subprocess.run(fixmate_cmd, capture_output=True, text=True, check=False)
		if result.returncode != 0:
			logging.error(f"fixmate_cmd failed with code {result.returncode}. Stderr: {result.stderr.strip()}")
			raise RuntimeError(f"fixmate_cmd failed: {result.stderr.strip()}")
		temp_files_to_clean.append(fixmate_temp_bam)

		# Clean up previous intermediate file
		if os.path.exists(sortname_temp_bam):
			os.remove(sortname_temp_bam)
			logging.info(f"Removed intermediate file: {sortname_temp_bam}")
			if sortname_temp_bam in temp_files_to_clean:
				temp_files_to_clean.remove(sortname_temp_bam)

		# Step 4: Sort by position to prepare for marking/removing duplicates
		sortpos_cmd = [
			'samtools', 'sort',
			'-m', str(args.memory)+'G', # Max memory per thread for sorting
			'-@', str(args.threads), # Number of threads for sorting
			fixmate_temp_bam, # Input BAM
			'-o', markdup_bam # Output to final file
		]
		logging.info(f"running sortpos command started at {time.strftime('%H:%M:%S')}")
		result = subprocess.run(sortpos_cmd, capture_output=True, text=True, check=False)
		if result.returncode != 0:
			logging.error(f"sortpos_cmd failed with code {result.returncode}. Stderr: {result.stderr.strip()}")
			raise RuntimeError(f"Pipeline failed at sortpos_cmd: {result.stderr.strip()}")
		# Clean up previous intermediate file
		if os.path.exists(fixmate_temp_bam):
			os.remove(fixmate_temp_bam)
			logging.info(f"Removed intermediate file: {fixmate_temp_bam}")
			if fixmate_temp_bam in temp_files_to_clean:
				temp_files_to_clean.remove(fixmate_temp_bam)

		# Step 5: Index the deduplicated BAM file
		index_cmd = ['samtools', 'index', markdup_bam]
		logging.info(f"Running index_cmd: {' '.join(index_cmd)}")
		result_index1 = subprocess.run(index_cmd, capture_output=True, text=True, check=False)
		if result_index1.returncode != 0:
			logging.error(f"index_cmd failed with code {result_index1.returncode}. Stderr: {result_index1.stderr.strip()}")
			raise RuntimeError(f"index_cmd failed: {result_index1.stderr.strip()}")
		
		# Step 6: Filter out duplicate reads (if marked) and specific chromosomes (e.g., mitochondrial)
		# This 'if' block should be at the same indentation level as the 'index_cmd' block above it.
		if os.path.isfile(markdup_bam):
			rmdup_cmd = [
				'samtools', 'view', '-@', str(args.threads),
				'-b', # Output BAM
				'-f', '3', # Reads must be properly paired (0x3)
				'-F', '1024', # Filter out duplicate reads (0x400)
				markdup_bam # Input BAM
			]
			rmdup_cmd.extend(['chr{}'.format(c) for c in list(map(str, range(1,23))) + ['X','Y']])

			with open(rmdup_bam, 'wb') as bam_out: # Use 'wb' for binary write mode
				logging.info(f"Running rmdup_cmd: {' '.join(rmdup_cmd)}")
				result_rmdup = subprocess.run(rmdup_cmd, stdout=bam_out, stderr=subprocess.PIPE, check=False)
				if result_rmdup.returncode != 0:
					logging.error(f"rmdup_cmd failed with code {result_rmdup.returncode}. Stderr: {result_rmdup.stderr.decode().strip()}")
					raise RuntimeError(f"rmdup_cmd failed: {result_rmdup.stderr.decode().strip()}")

			# Index the final filtered BAM file
			index_cmd2 = ['samtools', 'index', rmdup_bam]
			logging.info(f"Running index_cmd2: {' '.join(index_cmd2)}")
			result_index2 = subprocess.run(index_cmd2, capture_output=True, text=True, check=False)
			if result_index2.returncode != 0:
				logging.error(f"index_cmd2 failed with code {result_index2.returncode}. Stderr: {result_index2.stderr.strip()}")
				raise RuntimeError(f"index_cmd2 failed: {result_index2.stderr.strip()}")
		else: # This 'else' correctly pairs with `if os.path.isfile(markdup_bam):`
			raise FileNotFoundError('{} not found!'.format(markdup_bam))
		
	except Exception as e:
		logging.critical(f"An error occurred during remove_duplicate_reads: {e}")
		# Re-raise the exception to terminate the script and mark the SLURM job as failed
		raise 
	finally:
		# This block ensures that all temporary files are cleaned up,
		# even if an error occurred during processing.
		for f in temp_files_to_clean:
			if os.path.exists(f):
				os.remove(f)
				logging.info(f"Cleaned up temporary file: {f}")

	logging.info(f"def remove_duplicate_reads finished at {time.strftime('%H:%M:%S')}")
	return

				





def qc_metrics(args):
	logging.info(f"qc_metrics started at {time.strftime('%H:%M:%S')}")
	md_bam = args.output_prefix + '.filt.md.bam'
	rmdup_bam = args.output_prefix + '.filt.rmdup.bam'
	tagalign_file = args.output_prefix + '.filt.rmdup.tagAlign.gz'
	
	if os.path.isfile(rmdup_bam) and not os.path.isfile(tagalign_file):
		with gzip.open(tagalign_file, 'wt') as f:
			bamfile = pysam.AlignmentFile(rmdup_bam, 'rb')
			genome_size = {item['SN']:item['LN'] for item in bamfile.header['SQ']}
			for read in bamfile:
				if not read.is_proper_pair:
					continue
				barcode = read.query_name.split('_')[0]
				read_chr = read.reference_name
				read_start = max(1, read.reference_end - args.shift - args.extsize - 5 if read.is_reverse else read.reference_start + args.shift + 4)
				read_end = min(genome_size[read_chr], read.reference_end - args.shift - 5 if read.is_reverse else read.reference_start + args.shift + args.extsize + 4)
				read_qual = read.mapping_quality
				if read.is_reverse:
					read_orient = '-'
				else:
					read_orient = '+'
				#line_out = '\t'.join([read_chr, str(read_start), str(read_end), '{}_{}'.format(args.name,barcode), str(read_qual), read_orient])
				line_out = '\t'.join([read_chr, str(read_start), str(read_end), barcode, str(read_qual), read_orient])
				print(line_out, sep='\t', file=f)
			bamfile.close()
#	else:
#		raise FileNotFoundError('{} not found!'.format(rmdup_bam))
	
	# qc_metrics = {}
	# chr_names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
	# if os.path.isfile(md_bam):
	# 	bamfile = pysam.AlignmentFile(md_bam, 'rb')
	# 	for read in bamfile:
	# 		barcode = '{}_{}'.format(args.name, read.query_name.split('_')[0])
	# 		if barcode not in qc_metrics:
	# 			qc_metrics[barcode] = {}
	# 		if not read.is_duplicate:
	# 			if read.reference_name in chr_names: 
	# 				qc_metrics[barcode]['unique_usable_reads'] = qc_metrics[barcode].get('unique_usable_reads', 0) + 1
	# 			elif read.reference_name == 'chrM':
	# 				qc_metrics[barcode]['unique_mito_reads'] = qc_metrics[barcode].get('unique_mito_reads', 0) + 1
	# 			qc_metrics[barcode]['total_sequenced_reads'] = qc_metrics[barcode].get('total_sequenced_reads', 0) + 1
	# 		else:
	# 			qc_metrics[barcode]['duplicated_reads'] = qc_metrics[barcode].get('duplicated_reads', 0) + 1
	# 			qc_metrics[barcode]['total_sequenced_reads'] = qc_metrics[barcode].get('total_sequenced_reads', 0) + 1
	# else:
	# 	raise FileNotFoundError('{} not found!'.format(md_bam))
	# qc_metrics = pd.DataFrame.from_dict(qc_metrics, orient='index').fillna(0).astype(int)
	# macs2_cmd = ['macs2', 'callpeak', '-t', tagalign_file, '--outdir', args.output, '-n', args.name, '-p', '.05', '--nomodel', '--keep-dup', 'all', '--shift', '0', '--extsize', '200', '-g', 'hs']
	# with open(os.devnull, 'w') as f:
	# 	subprocess.call(macs2_cmd, stderr=f)
	# try:
	# 	os.remove(args.output_prefix + '_peaks.xls')
	# 	os.remove(args.output_prefix + '_summits.bed')
	# except:
	# 	pass
	# blacklist_cmd = subprocess.Popen(['bedtools', 'intersect', '-a' , args.output_prefix + '_peaks.narrowPeak', '-b', args.blacklist_file, '-v'], stdout=subprocess.PIPE)
	# intersect_cmd = subprocess.Popen(['bedtools', 'intersect', '-a', tagalign_file, '-b', '-' ], stdin=blacklist_cmd.stdout, stdout=subprocess.PIPE)
	# peak_counts = {bc:0 for bc in qc_metrics.index}
	# for line in intersect_cmd.stdout:
	# 	line = line.decode()
	# 	fields = line.rstrip().split('\t')
	# 	peak_counts[fields[3]] += 1
	# qc_metrics['reads_in_peaks'] = qc_metrics.index.map(peak_counts).fillna(0).astype(int)
	# tss_counts = {bc:0 for bc in qc_metrics.index}
	# tss_used = {bc:[] for bc in qc_metrics.index}
	# tss_cmd = subprocess.Popen(['bedtools', 'intersect', '-a', tagalign_file, '-b', args.promoter_file, '-wa', '-wb'], stdout=subprocess.PIPE)
	# for line in tss_cmd.stdout:
	# 	line = line.decode()
	# 	fields = line.rstrip().split('\t')
	# 	tss_counts[fields[3]] += 1
	# 	tss_used[fields[3]].append(fields[9])
	# qc_metrics['reads_in_promoters'] = qc_metrics.index.map(tss_counts).fillna(0).astype(int)
	# qc_metrics['tss_used'] = [len(set(tss_used[bc])) for bc in qc_metrics.index]
	# total_prom = len(open(args.promoter_file).read().splitlines())
	# qc_metrics['frac_reads_in_peaks'] = qc_metrics['reads_in_peaks'].div(qc_metrics['unique_usable_reads']).replace(np.inf, 0).fillna(0)
	# qc_metrics['frac_reads_in_promoters'] = qc_metrics['reads_in_promoters'].div(qc_metrics['unique_usable_reads']).replace(np.inf, 0).fillna(0)
	# qc_metrics['frac_promoters_used'] = qc_metrics['tss_used']/total_prom
	# qc_metrics['frac_mito_reads'] = qc_metrics['unique_mito_reads'].div(qc_metrics['unique_usable_reads'] + qc_metrics['unique_mito_reads']).replace(np.inf, 0).fillna(0)
	# qc_metrics['frac_duplicated_reads'] = qc_metrics['duplicated_reads'].div(qc_metrics['total_sequenced_reads']).fillna(0)
	# qc_metrics.to_csv(os.path.join(args.output_prefix + '.qc_metrics.txt'), sep='\t')
	return	

def generate_matrix(args):
	logging.info(f"generate_matrix started at {time.strftime('%H:%M:%S')}")
	tagalign_file = args.output_prefix + '.filt.rmdup.tagAlign.gz'
	#qc_metrics = pd.read_table(args.output_prefix + '.qc_metrics.txt', sep='\t', header=0, index_col=0)
	#pass_barcodes = qc_metrics.loc[qc_metrics['unique_usable_reads'] >= args.minimum_reads].index 
	pass_barcodes = open(args.keep_bc).read().splitlines()

	lf_mtx_file = args.output_prefix + '.long_fmt_mtx.txt.gz'
	barcodes_file = args.output_prefix + '.barcodes'
	regions_file = args.output_prefix + '.regions'
	mtx_file = args.output_prefix + '.mtx'
	
	windows_file = make_windows(args)
	window_intersect = intersect_regions(tagalign_file, windows_file)
	cut = subprocess.Popen(['cut', '-f', '4,10'], stdin=window_intersect.stdout, stdout=subprocess.PIPE)
	sort = subprocess.Popen(['sort', '-S', '{}G'.format(args.memory * args.threads)], stdin=cut.stdout, stdout=subprocess.PIPE)
	uniq = subprocess.Popen(['uniq', '-c'], stdin=sort.stdout, stdout=subprocess.PIPE)
	awk = subprocess.Popen(['awk', '''BEGIN{{OFS="\\t"}} {{print $2,$3,$1}}'''], stdin=uniq.stdout, stdout=subprocess.PIPE)
	with gzip.open(lf_mtx_file, 'wt') as f:
		subprocess.call(['gzip', '-c'], stdin=awk.stdout, stdout=f)

	lf_mtx = pd.read_table(lf_mtx_file, sep='\t', header=None, names=['barcode','region','count'])
	lf_mtx = lf_mtx.loc[lf_mtx['barcode'].isin(pass_barcodes)]
	lf_mtx.to_csv(lf_mtx_file, sep='\t', header=False, index=False, compression='gzip')

	tmp_R = args.output_prefix + '.tmp.R'
	with open(tmp_R, 'w') as tR:
		print('''library(Matrix)''', file=tR)
		print('''library(methods)''', file=tR)      
		print('''data <- read.table('{}', sep='\\t', header=FALSE)'''.format(lf_mtx_file), file=tR)
		print('''sparse.data <- with(data, sparseMatrix(i=as.numeric(as.factor(V1)), j=as.numeric(as.factor(V2)), x=V3, dimnames=list(levels(as.factor(V1)),levels(as.factor(V2)))))''', file=tR)
		print('''t <- writeMM(sparse.data, '{}')'''.format(mtx_file), file=tR)
		print('''writeLines(rownames(sparse.data), '{}')'''.format(barcodes_file), file=tR)
		print('''writeLines(colnames(sparse.data), '{}')'''.format(regions_file), file=tR)
	subprocess.call(['Rscript', tmp_R])
	subprocess.call(['gzip', mtx_file])
	os.remove(windows_file)
	os.remove(tmp_R)
	return

def intersect_regions(tagalign, regions):
	logging.info(f"intersect_regions started at {time.strftime('%H:%M:%S')}")
	awk_cmd = ['awk', '''BEGIN{{FS=OFS="\\t"}} {{peakid=$1":"$2"-"$3; gsub("chr","",peakid); print $1, $2, $3, peakid}}''', regions]
	intersect_cmd = ['bedtools', 'intersect', '-a', tagalign, '-b', '-', '-wa', '-wb']
	awk = subprocess.Popen(awk_cmd, stdout=subprocess.PIPE)
	intersect = subprocess.Popen(intersect_cmd, stdin=awk.stdout, stdout=subprocess.PIPE)
	return intersect

def make_windows(args):
	logging.info(f"make_windows started at {time.strftime('%H:%M:%S')}")
	makewindows_cmd = ['bedtools', 'makewindows', '-g', args.chrom_sizes, '-w', str(args.window_size * 1000)]
	filter_cmd = ['grep', '-v', '_']
	blacklist_cmd = ['bedtools', 'intersect', '-a', '-', '-b', args.blacklist_file, '-v']
	windows_file = '{}.{}kb_windows.bed'.format(args.output_prefix, args.window_size)
	with open(windows_file, 'w') as f:
		makewindows = subprocess.Popen(makewindows_cmd, stdout=subprocess.PIPE)
		filt = subprocess.Popen(filter_cmd, stdin=makewindows.stdout, stdout=subprocess.PIPE)
		subprocess.call(blacklist_cmd, stdin=filt.stdout, stdout=f)
	return windows_file

def main(args):
	logging.info('Start.')
	if not os.path.isdir(args.output):
		os.makedirs(args.output)
	args.output_prefix = os.path.join(args.output, args.name)
	if not args.skip_convert:
		convert_10X_bam(args)
	if not args.skip_rmdup:
		logging.info('Removing duplicate and mitochrondrial reads.'.format(args.minimum_reads))
		remove_duplicate_reads(args)
	if not args.skip_qc:                
		qc_metrics(args)
	if not args.skip_matrix:
		logging.info('Generating tagalign and chromatin accessibility matrix.')
		generate_matrix(args)
	logging.info('Finish.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Use 10X output to process snATAC-seq data.')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-b', '--input-bam', required=True, type=str, help='Position sorted bam from 10X output')
	io_group.add_argument('-o', '--output', required=False, type=str, default=os.getcwd(), help='Output directory to store processed files')
	io_group.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')

	align_group = parser.add_argument_group('Alignment arguments')
	align_group.add_argument('-t', '--threads', required=False, type=int, default=8, help='Number of threads to use for alignment [8]')
	align_group.add_argument('-m', '--memory', required=False, type=int, default=4, help='Maximum amount of memory (G) per thread for samtools sort [4]')
	align_group.add_argument('-q', '--map-quality', required=False, type=int, default=30, help='Mapping quality score filter for samtools [30]')
	align_group.add_argument('-ref', '--reference', required=False, type=str, default='/group/soranzo/paola.benaglio/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa', help='Path to the reference genome')
	
	#dup_group = parser.add_argument_group('Remove duplicates arguments')
	#dup_group.add_argument('--picard', required=False, type=str, default='/home/joshchiou/bin/picard.jar', help='Path to picard.jar')
	
	matrix_group = parser.add_argument_group('Matrix generation arguments')
	matrix_group.add_argument('--shift', required=False, type=int, default=-100, help='Read shift length')
	matrix_group.add_argument('--extsize', required=False, type=int, default=200, help='Read extension size')
	matrix_group.add_argument('--minimum-reads', required=False, type=int, default=500, help='Minimum number of reads for barcode inclusion')
	matrix_group.add_argument('--minimum-frip', required=False, type=float, default=0, help='Minimum frip for barcode inclusion')
	matrix_group.add_argument('--window-size', required=False, type=int, default=5, help='Size (kb) to use for defining windows of accessibility')
	matrix_group.add_argument('--chrom-sizes', required=False, type=str, default='/group/soranzo/paola.benaglio/references/hg38.chrom.sizes', help='Chromosome sizes file from UCSC')
	matrix_group.add_argument('--blacklist-file', required=False, type=str, default='/group/soranzo/paola.benaglio/references/hg38-blacklist.v2.bed', help='BED file of blacklisted regions')
	matrix_group.add_argument('--promoter-file', required=False, type=str, help='BED file of autosomal promoter regions')
	matrix_group.add_argument('--keep-bc', required=True, type=str, help='List of barcodes to keep from qc steps')

	skip_group = parser.add_argument_group('Skip steps')
	skip_group.add_argument('--skip-convert', required=False, action='store_true', default=False, help='Skip bam conversion step')
	skip_group.add_argument('--skip-rmdup', required=False, action='store_true', default=False, help='Skip duplicate removal step')
	skip_group.add_argument('--skip-qc', required=False, action='store_true', default=False, help='Skip quality metrics step')
	skip_group.add_argument('--skip-matrix', required=False, action='store_true', default=False, help='Skip matrix generation step')
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)
	args = process_args()
	main(args)

