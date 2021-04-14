import sys, os, subprocess
import argparse

def main():
	parser = argparse.ArgumentParser(prog="tab2gb_tbl")
	parser = argparse.ArgumentParser(description='Given a file *geneAnnotations.tab it transforms it into another tabular file for submission in genbank')
	parser.add_argument('-file', '--TabFile', help='Give this command a tab file to be converted into genbanks submission tbl')
	args = parser.parse_args()
	return args

def preprocess(file):
	with open(args.TabFile,'r') as f:
		lines = f.readlines()
	lines = [line.replace(' ', '') for line in lines]
	return lines

def get_arguments(lines):
	tbl_file = open(str(args.TabFile)[:-20] + '_gb.tbl','a+')
	tbl_file.write('>Feature ' + str(args.TabFile)[:-20] + '\n')
	tbl_file.close()
	for line in lines:
		if line.count('\t') == 5 and '---' not in line:
			molecule, type_, gene, start, stop, strand = line.split('\t')
		elif line.count('\t') != 5 and '---' not in line:
			sys.exit('There seems to be excessive hidden characters in this file. Please remove the unnecessary tabulations and run again')
		elif '---' in line:
			sys.exit('\nAll done!')
		gene, product = get_names(gene, type_)
		make_tbl(molecule, type_, gene, start, stop, strand, product)

def make_tbl(molecule,type_,gene,start,stop,strand,product):
	if strand == '1\n':
		tbl_file = open(str(args.TabFile)[:-20] + '_gb.tbl','a')
		tbl_file.write(start + '\t' + stop + '\t' 'gene' + '\n')
		tbl_file.write('\t'*3 + 'gene' + '\t' + gene + '\n')
		tbl_file.write(start + '\t' + stop + '\t' + type_ + '\n')
		tbl_file.write('\t'*3 + 'product' + '\t' + product + '\n')
		tbl_file.close()
	else:
		tbl_file = open(str(args.TabFile)[:-20] + '_gb.tbl','a')
		tbl_file.write(stop + '\t' + start + '\t' 'gene' + '\n')
		tbl_file.write('\t'*3 + 'gene' + '\t' + gene + '\n')
		tbl_file.write(stop + '\t' + start + '\t' + type_ + '\n')
		tbl_file.write('\t'*3 + 'product' + '\t' + product + '\n')
		tbl_file.close()

def get_names(gene,type_):
	if type_ == 'tRNA' or type_ == 'rRNA':
		rna_dict = {'trnD':'trnD(gac)', 'trnG':'trnG(gga)', 'trnL2':'trnL2(tta)', 'trnV':'trnV(gta)', 'trnI':'trnI(atc)',
		'trnC':'trnC(tgc)', 'trnQ':'trnQ(caa)', 'trnF':'trnF(ttc)', 'trnP':'trnP(cca)', 'trnN':'trnN(aac)', 'trnL1':'trnL1(cta)',
		'trnY':'trnY(tac)', 'trnT':'trnT(aca)', 'trnK':'trnK(aaa)', 'trnR':'trnR(cga)', 'trnW':'trnW(tga)', 'trnM':'trnM(atg)',
		'trnE':'trnE(gaa)', 'trnS1':'trnS1(aga)', 'trnS2':'trnS2(tca)', 'trnA':'trnA(gca)' ,'trnH':'trnH(cac)',
		'rrnS':'12S ribosomal RNA', 'rrnL':'16S ribosomal RNA'}
		gene = rna_dict[gene]
	product_dict = {'cox1':'cytochrome c oxidase subunit 1', 'cox3':'cytochrome c oxidase subunit 3',
	'atp6':'ATP synthase F0 subunit 6', 'trnD(gac)':'tRNA-ASP', 'atp8':'ATP synthase F0 subunit 8',
	'nad4l':'NADH dehydrogenase subunit 4L', 'nad4':'NADH dehydrogenase subunit 4', 'nad6':'NADH dehydrogenase subunit 6',
	'trnG(gga)':'tRNA-GLY', 'nad1':'NADH dehydrogenase subunit 1', 'trnL2(tta)':'tRNA-LEU2', 'trnV(gta)':'tRNA-Val',
	'trnI(atc)':'tRNA-Ile', 'trnC(tgc)':'tRNA-Cys', 'trnQ(caa)':'tRNA-Gln', 'nad5':'NADH dehydrogenase subunit 5',
	'trnF(ttc)':'tRNA-Phe', 'cob':'cytochrome b', 'trnP(cca)':'tRNA-Pro', 'trnN(aac)':'tRNA-Asn', 'trnL1(cta)':'tRNA-Leu',
	'16S ribosomal RNA':'l-rRNA', 'trnY(tac)':'tRNA-Tyr', 'trnT(aca)':'tRNA-Thr','trnK(aaa)':'tRNA-Lys', '12S ribosomal RNA':'s-rRNA',
	'trnR(cga)':'tRNA-Arg', 'trnW(tga)':'tRNA-Trp', 'trnM(atg)':'tRNA-Met', 'nad2':'NADH dehydrogenase subunit 2', 'trnE(gaa)':'tRNA-Glu',
	'trnS1(aga)':'tRNA-Ser', 'trnS2(tca)':'tRNA-Ser', 'trnA(gca)':'tRNA-Ala', 'trnH(cac)':'tRNA-His','nad3':'NADH dehydrogenase subunit 3',
	'cox2':'cytochrome c oxidase subunit 2'}
	product = product_dict[gene]
	return gene, product

if __name__ == "__main__":
	args = main()
	cleared_lines = preprocess(args)
	get_arguments(cleared_lines)
