import sys,os, subprocess
import argparse

def parse_args():
	parser = argparse.ArgumentParser(prog="aggregate_blastn")
	parser = argparse.ArgumentParser(description='Given a txt file with the species names for every line, it reads the information and retrieves the sequences of interest')
	parser.add_argument('-blast_directory', '--Dir', help='Give this command a directory containing all blastn file outputs')
	parser.add_argument('-number_seqs', '--N_seqs', help='Give this command a file with all the species names')
	args = parser.parse_args()
	return args

def main(args):
	original_path = os.getcwd()
	os.chdir(args.Dir)
	final_file = open(os.path.join(original_path,'blastn.xlsx'),'a+')
	final_file.write('Sample\tContig\tSpecies Name\tKingdom\tTitle,Score\n')
	final_file.close()
	for file in sorted(os.listdir()):
		current_blastn = open(file,'r')
		i = 0
		j = 1
		final_file = open(os.path.join(original_path,'blastn.xlsx'),'a+')
		final_file.write(file[:-7])
		for lines in current_blastn:
			lines = lines.replace(',',';')
			i += 1
			if i == 1:
				final_file.write('\tContig ' + str(j) + '\t' + lines)
			if i == int(args.N_seqs):
				j += 1
				i = 0
				final_file.write('' + '\t ' + '\t' + lines)
			elif i != args.N_seqs and i != 1: 
				final_file.write('' + '\t ' + '\t' + lines)
		current_blastn.close()
		final_file.close()

if __name__ == "__main__":
	args = parse_args()
	main(args)
