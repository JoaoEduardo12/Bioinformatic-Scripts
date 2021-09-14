import sys, os, subprocess
import argparse
from ete3 import EvolTree

tree = EvolTree("tree.nw", binpath = "/home/edu/miniconda3/envs/ete3/bin/ete3_apps/bin")
tree.link_to_alignment("infile.phy", alg_format = "phylip")
tree.workdir = os.getcwd()

print(tree)

print('running model M0, for comparison with branch-site models...')

tree.run_model('M0', keep = True)
#tree.link_to_evol_model("/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Phylogenetic_Tree","M0")
chimaeriformes =tree.get_common_ancestor("HM147138.1","HM147135.1")
#chimaeriformes =tree.get_common_ancestor("Human_ECP","Goril_ECP")

for leaf in chimaeriformes:
    tree.mark_tree([leaf.node_id], marks = ["#1"])
#tree.run_model("bsA." + chimaeriformes)
#tree.mark_tree([leaf.node_id], marks = ["#1"])
print("Running")
print(tree.write())
tree.run_model('bsA.Chimaeriformes')
tree.run_model("bsA1.Chimaeriformes")

print ('p-value of positive selection for sites on this branch is: ')
ps = tree.get_most_likely('bsA.Chimaeriformes' , 'bsA1.Chimaeriformes')
print(str(ps))
rx = tree.get_most_likely('bsA1.Chimaeriformes', 'M0')
print(str(rx))
model = tree.get_evol_model("bsA.Chimaeriformes")
if ps < 0.05 and float(model.classes['foreground w'][2]) > 1:
    print ('we have positive selection on sites on this branch')
    tree.show(histfaces=['bsA1.Chimaeriformes'])
elif rx<0.05 and ps>=0.05:
    print ('we have relaxation on sites on this branch')
else:
    print ('no signal detected on this branch, best fit for M0')
#tree.show(histfaces=['bsA1.'])

for models in tree._models:
    print(tree.get_evol_model(models))

from _pickle import dump

#out = open('my_tree.pik', 'w')
#dump(tree, out)
#out.close()

#tree.mark_tree(map(lambda x: x.node_id, tree.get_descendants()),
#                    marks=[''] * len(tree.get_descendants()), verbose=True)

print("The End.")
