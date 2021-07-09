import sys, os, subprocess
import argparse
from Bio import SeqIO

def parse_args():
	"""Parsing arguments via terminal"""
	parser = argparse.ArgumentParser(prog="remove_duplicates")
	parser = argparse.ArgumentParser(description='To remove duplicates from the dataset so that only one species is present')
	parser.add_argument('-file', '--File', help='Give this command a file with accession numbers marked with hashtag beforehand (anywwhere)')
	parser.add_argument('-csv', '--Dataset', help = 'Give this command the original dataset where the final results are')
	parser.add_argument('-genes', '--GeneDir', help = 'Give this command where the directory with the genes are (only will accept fasta files that end with ".fn"')
	parser.add_argument('-out', '--OutDir',
						help='Give this command the output directory (files will have the same name but with the prefix "uniques"')
	args = parser.parse_args()
	return args

def main(args):
	os.mkdir(args.OutDir)
	records = []
	original_dir = os.getcwd()
	with open(args.File) as f:
		lines = f.readlines()
		f.close()
	for line in lines:
		if '#' in line:
			records.append(line.strip('\n')[1:])
	os.chdir(args.GeneDir)
	for file in os.listdir():
		if file.endswith('.fn'):
			current_gene = SeqIO.parse(file, 'fasta')
			for record in current_gene:
				if record.id not in records:
					new_file = open(original_dir + '/' + args.OutDir + '/unique_{}.fn' .format(file[:-3]), 'a+')
					new_file.write('>{}\n{}\n' .format(record.id, str(record.seq)))
					new_file.close
	os.chdir(original_dir)
	with open(args.Dataset) as f:
		lines = f.readlines()
		for line in lines:
			if line.split(',')[1] not in records:
				new_file = open(original_dir + '/' + args.OutDir + '/unique_{}' .format(args.Dataset), 'a+')
				new_file.write(line)
				new_file.close()

if __name__ == "__main__":
	args = parse_args()
	main(args)