import sys,os, subprocess
import argparse 
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(prog="Reorganize_Sequence")
parser = argparse.ArgumentParser(description='Performs a ClustalW alignment between the whole mitogenomic sequence and the COI sequence and reorganizes the first to start at the first position of the alignment. Example use: python3 /home/edu/Desktop/Mitogenomas/Reorganize_Sequence.py -mito AnoAnaIran.fas -coi /home/edu/Desktop/Mitogenomas/COI.fas')
parser.add_argument('-mito', '--WholeMitogenome', help='use this argument to insert the whole mitochondrial genome sequence')
parser.add_argument('-coi', '--COIsequence', help='insert in this argument the COI/COX1 sequence')
parser.add_argument('-rc', '--ReverseComplement', help='Does the reverse complement if the option is added')
args = parser.parse_args()

mitogenome = open(args.WholeMitogenome).read()
coi = open(args.COIsequence,'r')

def reverse_complement(file_):
	file__ = open(file_,'r')
	seq_ = SeqIO.read(file_,'fasta')
	seq_ = str(seq_.seq)
	seq = Seq(seq_)
	file = open('RC_'+str(args.WholeMitogenome),'a+')
	file.write(file__.readlines()[0])
	file.write(str(seq.reverse_complement()))
	file__.close()
	return file.close()

def create_infile(seq1,seq2):
	infile = open('infile.fasta','w')
	infile.write(str(seq1))
	infile.close()
	infile = open('infile.fasta','a+')
	infile.write('\n>COI\n')
	seq2 = [x[:-1] for x in seq2.readlines() if '>' not in x]
	seq2 = ''.join(seq2)
	infile.write(seq2)
	return infile

if args.ReverseComplement == None:
	infile = create_infile(mitogenome,coi)
else:
	rc_mitogenome = reverse_complement(args.WholeMitogenome)
	rc_mitogenome = open('RC_'+str(args.WholeMitogenome)).read()
	infile = create_infile(rc_mitogenome,coi)

infile.close()
coi.close()

cline = ClustalwCommandline("clustalw", infile="infile.fasta")
child = subprocess.call(str(cline),
  stdout=subprocess.PIPE,
  shell=(sys.platform!="win32"))

align = AlignIO.read("infile.aln", "clustal")
count = 0
for record in align:
	if record.id == 'COI':
		for letter in record.seq:
			if letter == '-':
				count += 1
			else:
				break


mitogenome = open(args.WholeMitogenome,'r')
final_file = open('R_' + str(args.WholeMitogenome), 'w')
final_file.write(mitogenome.readlines()[0])
final_file.close()
mitogenome.close()

if args.ReverseComplement == None:
	record = SeqIO.read(args.WholeMitogenome,'fasta')
else:
	record = SeqIO.read('RC_'+str(args.WholeMitogenome),'fasta')

final_file = open('R_' + str(args.WholeMitogenome), 'a+')
final_file.write(str(record.seq[count::]))
final_file.write(str(record.seq[0:count]))
final_file.close()



