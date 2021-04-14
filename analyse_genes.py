import sys, os, subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio import Entrez
Entrez.email='joao.edu.t@hotmail.com'

def main():
	parser = argparse.ArgumentParser(prog="gene_alignment")
	parser = argparse.ArgumentParser(description='Given a txt file with the accession numbers, it opens the corresponding genbank files and extract the gene features and proceeds to align each of the genes')
	parser.add_argument('-file', '--TxtFile', help='Give this command a txt file with all the accession numbers to extract straight seperated by paragraphs')
	parser.add_argument('-genes','--NewSequences', help='Give the directory where all the new genes from the assembled sequences reside')
	parser.add_argument('-align','--AlignmentAlgorithm', help='Give the multiple sequence alignment algorithm of your choosing: muscle, mafft or clustalw')
	parser.add_argument('-sequences','--FullSequences', help='Give the directory where all the complete set of sequences are.')
	parser.add_argument('-annotations','--GenomeAnnotations', help='Give the directory where all the annotations of the respective genes are')
	args = parser.parse_args()
	return args

def singleEntry(singleID): 
    '''Extracts the genbank file, with an accession number as input'''
    handle = Entrez.efetch(db='nucleotide',id=singleID, rettype = 'gb', retmode= 'text')
    file_name = '%s.gb' % singleID
    f = open('%s.gb' % singleID, 'w')
    f.write(handle.read())
    handle.close()
    f.close()
    return file_name

def create_gene_fasta(gene,sequence,name=None,direct='no'):
	'''It creates and continuosly appends individual gene sequences from genbank files or fasta files'''
	ATP6_list = ['atp6','ATP6','ATPase 6','ATPase6']
	ATP8_list = ['atp8','ATP8','ATPase 8','ATPase8']
	COI_list = ['cox1','COI','COX1','COXI']
	COII_list = ['cox2','COII','COX2']
	COIII_list = ['cox3','COIII','COX3']
	CYTB_list = ['cob','CYTB','cyt b','Cyt b']
	rest_list = ['ND1','ND2','ND3','ND4','ND4L','ND5','ND6']
	individual_list = ['nad1','nad2','nad3','nad4','nad4l','nad5','nad6']
	orfs = ['f-orf','F-ORF','F-orf','F-Orf','FORF','forf','Forf','HORF','h-orf','H-Orf','Horf','horf','H-orf']
	rnas = ['rrnL','rrnS']
	if gene in ATP6_list: #no case-switch in python...
		gene = ATP6_list[0]
	elif gene in ATP8_list:
		gene = ATP8_list[0]
	elif gene in COI_list:
		gene = COI_list[0]
	elif gene in COII_list:
		gene = COII_list[0]
	elif gene in COIII_list:
		gene = COIII_list[0]
	elif gene in CYTB_list:
		gene = CYTB_list[0]
	elif gene in rest_list:
		gene = gene
	elif gene in individual_list:
		#gene = gene[0] + gene[2:]
		#gene = gene.upper()
		gene = gene
	elif gene in rnas:
		gene = gene
	elif gene in orfs:
		gene = orfs[0]
	else:
		sys.exit("There seems to be a gene read (%s) that does not match the following: ATP8, ATP6, COI, COII, COIII, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6" % gene)
	if direct == 'no':
		file = open(gene + '.fn','a+')
		file.write(">%s\n%s\n" % (name,sequence))
	elif direct == 'yes':
		file = open(gene + '.fn','a+')
		file.write(">%s\n%s\n" % (name.id,sequence))

def check_ambiguity(seq,seq_id,fasta='no'):
	'''Replacing all the IUPAC ambiguity codes into N, so GUIDANCE can work properly'''
	if fasta == 'no':
		for i in range(len(seq)):
			if seq[i] != 'A' and seq[i] != 'T' and seq[i] != 'G' and seq[i] != 'C' and seq[i] != 'N':
				print('Processing ambiguities on %s, found the ambiguity code: %s; replacing it with N ...' % (seq_id,seq[i]))
				seq[i] = 'N'
	elif fasta == 'yes':
		uppercase_seq = seq.seq.upper()
		seq = uppercase_seq.tomutable()
		for i in range(len(seq)):
			if seq[i] != 'A' and seq[i] != 'T' and seq[i] != 'G' and seq[i] != 'C' and seq[i] != 'N' and seq[i] != '\n':
				print('Processing ambiguities on %s, found the ambiguity code: %s; replacing it with N ...' % (seq_id,seq[i]))
				seq[i] = 'N'
	return seq

def get_ncbi_cds(initial_file):
	'''accesses the txt file input and accesses the necessary functions to retrieve the gb file from acession number'''
	initial_file = open(initial_file).read()
	acc_numbs = list(initial_file.split('\n'))
	for seq_records in acc_numbs:
		gb_file = singleEntry(seq_records)
		for rec in SeqIO.parse(gb_file, "genbank"):
			print('Working with Genbank record %s' %  rec.id)
			if rec.features:
				for feature in rec.features:
					if feature.type == "CDS":
						gene_sequence = check_ambiguity(feature.location.extract(rec).seq.tomutable(),rec.id)
						create_gene_fasta(gene=feature.qualifiers['gene'][0].strip(''),
							name=rec.name,
							sequence=gene_sequence)
			else:
				print('There seems to be no features for the accession number %s' % rec.id)

def individual_genes(folder):
	'''From the output of the script CreateGeneSequences.pl, it reads the directory and finds the seperated
	gene files and appends then on the respective individual gene fasta created for all sequences.'''
	gene_list = ['atp6.fn','atp8.fn','cox1.fn','cox2.fn','cox3.fn','cob.fn','nad1.fn','nad2.fn','nad3.fn','nad4.fn','nad4l.fn','nad5.fn','nad6.fn','rrnL.fn','rrnS.fn']
	aa_list = ['atp6.fa','atp8.fa','cox1.fa','cox2.fa','cox3.fa','cob.fa','nad1.fa','nad2.fa','nad3.fa','nad4.fa','nad4l.fa','nad5.fa','nad6.fa']	
	print('Working on individual gene fastas')
	for i in range(len(gene_list)):
		for filename in os.listdir(folder):
			if filename == gene_list[i]:
				for record in SeqIO.parse(os.path.join(str(folder),filename),'fasta'):
					final_seq = check_ambiguity(record,record.id,fasta='yes')
					if str(gene_list[i][:-3]) + '.fa' in aa_list:
						analyse_seq(final_seq, gene_list[i][:-3], record.id)
					create_gene_fasta(gene=filename[:-3],
						sequence=final_seq,
						name=record,
						direct='yes')

def mafft(aa='no',gene=None,sequence=None):
	'''Performs multiple alignment algorithm mafft into the individual gene files'''
	if aa == 'no':
		print('\n\nInitiating mafft multiple sequence alignment!')
		file_list = ['cox1','cox2','cox3','cob','atp6','atp8','nad1','nad2','nad3','nad4','nad4l','nad5','nad6','rrnL','rrnS']
		for name in file_list:
			cline = MafftCommandline(input = name + '.fn')
			stdout, stderr = cline()
			with open(name + '.aln', "w") as handle:
				handle.write(stdout)
			align = AlignIO.read(name + ".aln", "fasta")

def clustalw(aa='no'):
	'''Performs multiple sequence alignment algorithm clustalw into the individual gene files'''
	if aa == 'no':
		print('Initiating clustalw multiple sequence alignment!')
		file_list = ['cox1','cox2','cox3','cob','atp6','atp8','nad1','nad2','nad3','nad4','nad4l','nad5','nad6','rrnL','rrnS']
		for name in file_list:
			cline = ClustalwCommandline("clustalw", infile=name + '.fn')
			child = subprocess.call(str(cline),
	  		  stdout=subprocess.PIPE,
	  		  shell=(sys.platform!="win32"))
			align = AlignIO.read(name + '.aln','clustalw')

def analyse_seq(gene_sequence,gene_name,name):
	boolean = 0
	if len(gene_sequence) % 3 != 0:
		print('!WARNING!: gene %s from %s appears to have length not multiple of 3' % (gene_name, name))
		boolean = 1
	rna = transcribe(gene_sequence)
	aminoacids = translate_rna(rna)
	count = 0
	for i in range(len(aminoacids)-1):
		if aminoacids[i] == '*':
			count += 1
	if count != 0:
		print('!WARNING!: gene %s from %s appears to have %d stop codon(s) midway' % (gene_name, name, count))
		boolean = 1
	#if boolean == 1:
		#print('Performing pairwise alignment of the gene with its genome')
		#mafft()

def transcribe(sequence):
	sequence = str(sequence)
	rna_seq = sequence.replace('T', 'U')
	return(rna_seq)

def translate_rna(sequence):
    codon2aa = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N", 
                "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T", 
                "AGA":"S", "AGC":"S", "AGG":"S", "AGU":"S", 
                "AUA":"M", "AUC":"I", "AUG":"M", "AUU":"I", 

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H", 
                "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P", 
                "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R", 
                "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L", 

                "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D", 
                "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A", 
                "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G", 
                "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V", 

                "UAA":"*", "UAC":"Y", "UAG":"*", "UAU":"T", 
                "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S", 
                "UGA":"W", "UGC":"C", "UGG":"W", "UGU":"C", 
                "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}
    protein_seq = ''
    for n in range(0, len(sequence), 3):
        if sequence[n:n+3] in codon2aa:
            protein_seq += codon2aa[sequence[n:n+3]]
    return protein_seq

if __name__ == "__main__":
	args = main()
	if args.TxtFile != None:
		get_ncbi_cds(args.TxtFile)
	if args.NewSequences != None:
		individual_genes(args.NewSequences)
	if args.AlignmentAlgorithm == 'mafft' or args.AlignmentAlgorithm == 'Mafft' or args.AlignmentAlgorithm == 'MAFFT':
		mafft()
	elif args.AlignmentAlgorithm == 'clustalw' or args.AlignmentAlgorithm == 'ClustalW':
		clustalw()
	elif args.AlignmentAlgorithm == 'muscle' or args.AlignmentAlgorithm == 'Muscle' or args.AlignmentAlgorithm == 'MUSCLE':
		muscle()
	elif args.AlignmentAlgorithm == None:
		sys.exit('.. All done! Consult the files in the directory')
	else:
		print('The alignment algorithm chosen is not valid, please select one of the following: mafft, clustalw or muscle.')
	print('.. All done! Consult the files in the directory')
