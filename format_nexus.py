import sys,os, subprocess
import argparse 

parser = argparse.ArgumentParser(prog="FormatNexus")
parser = argparse.ArgumentParser(description='With the output nexus file from SequenceMatrix, it creates an output for usage in PartitionFinder and IQTree; example usage: python format_nexus.py -file \'Sequence Matrix Block.nex\'')
parser.add_argument('-file','--SequenceMatrixFile',help = "insert the output file from SequenceMatrix")
parser.add_argument('-outdir','--Directory',help = "(optional) Name the directory you wish to output the pretended files")
args = parser.parse_args()

if args.Directory == None:
	args.Directory = os.getcwd() 

file=open(args.SequenceMatrixFile,'r')
partition_finder = open(os.path.join(str(args.Directory),'PartitionFinder_block.txt'),'a+')
output_iqtree = open(os.path.join(str(args.Directory),'IQTree_block.txt'),'a+')
output_iqtree.write('#Nexus\nbegin sets;')

for lines in file:
    if 'CHARSET a' in lines or 'CHARSET c' in lines or 'CHARSET n' in lines:
    	charset, gene, equal, position = lines.split(' ')
    	first, second = position.split('-')
    	output_iqtree.write('\ncharset '+gene+'_1' + ' = ' + str(int(first) + 0) + '-' + second[:-2] + '\\3;')
    	output_iqtree.write('\ncharset '+gene+'_2' + ' = ' + str(int(first) + 1) + '-' + second[:-2] + '\\3;')
    	output_iqtree.write('\ncharset '+gene+'_3' + ' = ' + str(int(first) + 2) + '-' + second[:-2] + '\\3;')
    	partition_finder.write(gene+'_1'+' = ' + str(int(first) + 0) + '-' + second[:-2] + '\\3;\n')
    	partition_finder.write(gene+'_2'+' = ' + str(int(first) + 1) + '-' + second[:-2] + '\\3;\n')
    	partition_finder.write(gene+'_3'+' = ' + str(int(first) + 2) + '-' + second[:-2] + '\\3;\n')
    elif 'CHARSET r' in lines:
    	charset, gene, equal, position = lines.split(' ')
    	output_iqtree.write('\ncharset ' + gene + ' = ' + position[:-2] + ';')
    	partition_finder.write(gene + ' = ' + position)

output_iqtree.write('\nend;')
print('\n\n..done!')
print('Output files written to: ',str(args.Directory))
file.close()
output_iqtree.close()
partition_finder.close()
