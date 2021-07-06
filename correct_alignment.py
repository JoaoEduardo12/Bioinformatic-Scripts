import sys, os, subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

def parse_args():
    """Parsing arguments"""
    parser = argparse.ArgumentParser(prog="correct_alignment")
    parser = argparse.ArgumentParser(
        description='Corrects gene sequences by deleting the stop codons, and replacing ambiguity codes with'
                    'N')
    parser.add_argument('-gene', '--GeneFile',
                        help='Give this command a file with the gene sequence')
    parser.add_argument('-file', '--File',
                        help='Give this command a txt file with each line with a file name')
    parser.add_argument('-dir', '--Directory',
                        help='Give this command the directory where the intended files are')
    parser.add_argument('-genetic_code', '--GeneCode',
                        help='Give this command the genetic code of the input sequences. Defaults to 11')
    parser.add_argument('-out', '--OutDir',
                        help='Give this command the output directory where the new file will be created')
    parser.add_argument('-mito', '--Mitochondrial',
                        help='Argument switch to tell the program its for mitochondrial genes',
                        default = 'standard')
    args = parser.parse_args()
    return args

def main(args):
    os.mkdir(args.OutDir)
    warning_gene, warning_stop, ambiguity_codes = [], [], []
    if args.Mitochondrial == 'standard':
        pcg = ['atp6','atp8','cox1','cox2','cox3','cob','cytb', 'nad1','nad2','nad3','nad4','nad4l','nad5','nad6']
    else: pcg = args.Mitochondrial
    if args.File:
        with open(args.File) as f:
            lines = f.readlines()
            f.close()
        for line in lines:
            for record in SeqIO.parse(os.path.join(args.Directory, line.strip('\n')), 'fasta'):
                final_seq = check_ambiguity(record, record.id, ambiguity_codes)
                if line.split('.')[0] in pcg:
                    final_seq = analyse_pcg(final_seq, line.split('.')[0], record.id, warning_gene, warning_stop)
                create_corrected_gene(gene = line.strip('\n'),
                                      sequence = final_seq,
                                      name = record.id,
                                      dir = args.OutDir)
    write_log(warn_gene = warning_gene, warn_stop = warning_stop, amb_code = ambiguity_codes, dir = args.OutDir)

def check_ambiguity(seq,seq_id,ambiguity_code):
    '''Replacing all the IUPAC ambiguity codes into N, so GUIDANCE can work properly'''
    corrected_seq = ''
    for i in range(len(seq)):
        if seq[i] != 'A' and seq[i] != 'T' and seq[i] != 'G' and seq[i] != 'C' and seq[i] != 'N':
            ambiguity_code.append('Processing ambiguities on %s, found the ambiguity code: %s; replacing it with N ...' % (seq_id,seq[i]))
            print('Processing ambiguities on %s, found the ambiguity code: %s; replacing it with N ...' % (seq_id,seq[i]))
            corrected_seq += 'N'
        else: corrected_seq += seq[i]
    return corrected_seq

def analyse_pcg(gene_sequence,gene_name,name, warning_gene, warning_stop):
    if len(gene_sequence) % 3 != 0:
        print('!WARNING!: gene %s from %s appears to have length not multiple of 3' % (gene_name, name))
        warning_gene.append('!WARNING!: gene %s from %s appears to have length not multiple of 3' % (gene_name, name))
        gene_sequence = correct_pcg(gene_sequence)
    aminoacids = str(Seq(gene_sequence).translate(table=args.GeneCode))
    count = 0
    for i in range(len(aminoacids)-1):
        if aminoacids[i] == '*':
            count += 1
    if count != 0:
        warning_stop.append('!WARNING!: gene %s from %s appears to have %d stop codon(s) midway' % (gene_name, name, count))
        print('!WARNING!: gene %s from %s appears to have %d stop codon(s) midway' % (gene_name, name, count))
    return gene_sequence


def correct_pcg(gene_sequence):
    if (len(gene_sequence) + 1) % 3 == 0:
        gene_sequence += 'A'
    elif (len(gene_sequence) + 2) % 3 == 0:
        gene_sequence += 'AA'
    return gene_sequence

def create_corrected_gene(gene = '', sequence = '', name = '', dir = ''):
    with open(os.path.join(dir, gene), 'a+') as f:
        f.write('>{}\n{}\n' .format(name, sequence[:-3]))
        f.close()

def write_log(warn_gene, warn_stop, amb_code, dir):
    with open(os.path.join(dir, 'log'), 'a+') as f:
        f.write('### WARNINGS ###\n\n\n\n#PROTEIN-CODING GENE LENGTH:\n\n{}\n\n'
                '#STOP CODONS:\n\n{}\n\n'
                '#AMBIGUITY CODES:\n\n{}\n\nScript runned successfully, check the corrected genes in {}'
                .format('\n'.join(map(str,warn_gene)), '\n'.join(map(str,warn_stop)),
                        '\n'.join(map(str, amb_code)), dir))
        f.close()

if __name__ == "__main__":
    args = parse_args()
    main(args)
