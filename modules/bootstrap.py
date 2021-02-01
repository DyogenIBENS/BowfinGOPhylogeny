"""

Module with function to generate and manipulate bootstrap replicates of the adjacency data.

"""

import os

import re
import itertools
import random

import numpy as np

from Bio import Phylo
from Bio.Phylo.Consensus import get_support
from ete3 import Tree



def bootstrap_matrix(all_adj, all_species, sp_dict):
    """
    Bootstrap the binary matrix 100 times and transform to 100 distance matrix.
    """
    bin_mat = np.array(all_adj, dtype=int)
    for k in range(100):

        dmat = np.zeros((len(all_species), len(all_species)))
        all_vect = range(len(all_adj))
        if (k+1)%10 == 0:
            print(f'{k+1} bootstrap replicates done')
        boot_indices = [random.choice(all_vect) for _ in range(len(all_adj))]
        boot_replicate = bin_mat[boot_indices, :]

        for i, j in list(itertools.combinations(all_species, 2)):
            i = all_species.index(i)
            j = all_species.index(j)


            nb_sim = len(np.logical_and(boot_replicate[:, i], boot_replicate[:, j]).nonzero()[0])

            max1 = len([a for a in boot_replicate[:, i] if a == 1])
            max2 = len([a for a in boot_replicate[:, j] if a == 1])
            dmat[i, j] = 1 - nb_sim/float(min(max1, max2))
            dmat[j, i] = dmat[i, j]

        os.makedirs("output/bootstrap/", exist_ok=True)

        #save matrix
        np.savetxt('output/bootstrap/dist_mat.txt', dmat,
                   header=' '.join([sp_dict[i.split()[0]].replace(' ', '_') for i in all_species]))

        #Format matrix for ape R package
        with open('output/bootstrap/dist_mat.txt', 'r') as infile,\
             open('output/bootstrap/dist_mat_'+str(k)+'.txt', 'w') as out:
            res = ''
            for i, line in enumerate(infile):
                if i == 0:
                    res += line[2:]
                else:
                    res += line
            out.write(res)


def add_bootstrap_support(tree, bootstrap_trees, root=True, branch_length_diff=None):
    """
    Add bootstrap support to the NJ tree.
    """
    diff_bl = {}


    #root all trees
    for treename in [tree] + bootstrap_trees:

        #reformat newick (sometimes R functions use comma for branch length values)
        with open(treename, 'r') as infile:
            mytree = infile.readlines()[0].strip()
            mytree = re.sub('0,([0-9])', '0.\\1', mytree)

        mytree = Tree(mytree)

        lca = mytree.get_common_ancestor(['Chicken', 'Xenopus'])
        if root:
            mytree.set_outgroup(lca)

        if not mytree.check_monophyly(values=['Gar', 'Bowfin'], target_attr="name")[0]:
            print(mytree)

        if branch_length_diff:
            for species_pair in branch_length_diff:
                sp1, sp2 = species_pair
                ancestor = mytree.get_common_ancestor([sp1, sp2])
                branch_sp1 = mytree.get_distance(ancestor, sp1)
                branch_sp2 = mytree.get_distance(ancestor, sp2)
                diff = branch_sp2 - branch_sp1
                diff_bl[species_pair] = diff_bl.get(species_pair, [])
                diff_bl[species_pair].append(diff)

        mytree.write(outfile=treename+'_reformatted')

    mytrees = [Phylo.read(mytree+'_reformatted', "newick") for mytree in bootstrap_trees]

    target = Phylo.read(tree+'_reformatted', 'newick')
    tree_with_support = get_support(target, mytrees)

    # Phylo.draw_ascii(majority_tree)
    Phylo.write(tree_with_support, 'output/bootstrap_tree.nwk', 'newick')

    #remove Bio.Phylo formatting that does not work with ete3 conventions...
    with open('output/bootstrap_tree.nwk', 'r') as infile:
        species_tree = infile.readlines()[0]
    species_tree = species_tree.replace('):', ')')

    with open('output/ete3_formatted_bootstrap_tree.nwk', 'w') as out:
        out.write(species_tree)

    return diff_bl
