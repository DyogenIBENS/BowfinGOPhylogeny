"""

Module with function to generate and manipulate bootstrap replicates of the adjacency data.

"""

import os

import re
import itertools

import numpy as np

from Bio import Phylo
from Bio.Phylo.Consensus import get_support
from ete3 import Tree



def bootstrap_matrix(all_adj, all_species, sp_dict, outfolder='output'):

    """
    Bootstraps the binary matrix 100 times and builds 100 corresponding distance matrices.
    Distance matrixes are saved to file in `output/bootstrap/`.

    Args:
        all_adj (list of list) : binary adjacency matrix, each sublist is a presence/absence vector,
                                 species index in columns are given by theh `all_species` argument.
        all_species (list) : list of species in the matrix (i.e column names)
        sp_dict (dict) : dict giving latin (key) to common (value) species name (for plot)

    Note (#FIXME):
        It is not very elegant that different functions compute the distance matrix for the original
        data (see `modules/matrix.make_distance_matrix`) and the boostrapped matrices (here). I've
        checked that implementations are identical (i.e same results for a given binary matrix),
        but having a single function would be better.
    """

    bin_mat = np.array(all_adj, dtype=int)
    for k in range(100):

        dmat = np.zeros((len(all_species), len(all_species)))
        # all_vect = range(len(all_adj))
        if (k+1)%10 == 0:
            print(f'{k+1} bootstrap replicates done')

        boot_indices = np.random.choice(len(all_adj), len(all_adj))
        boot_replicate = bin_mat[boot_indices, :]

        for i, j in list(itertools.combinations(all_species, 2)):
            i = all_species.index(i)
            j = all_species.index(j)


            nb_sim = len(np.logical_and(boot_replicate[:, i], boot_replicate[:, j]).nonzero()[0])

            max1 = sum(boot_replicate[:, i])
            max2 = sum(boot_replicate[:, j])

            dmat[i, j] = 1 - nb_sim/float(min(max1, max2))
            dmat[j, i] = dmat[i, j]

        os.makedirs(f"{outfolder}/bootstrap/", exist_ok=True)

        #save matrix
        np.savetxt(f'{outfolder}/bootstrap/dist_mat_tmp.txt', dmat,
                   header=' '.join([sp_dict[i.split()[0]].replace(' ', '_') for i in all_species]))

        #Format matrix for ape R package
        with open(f'{outfolder}/bootstrap/dist_mat_tmp.txt', 'r') as infile,\
             open(f'{outfolder}/bootstrap/dist_mat_'+str(k)+'.txt', 'w') as out:
            res = ''
            for i, line in enumerate(infile):
                if i == 0:
                    res += line[2:]
                else:
                    res += line
            out.write(res)
    print("DONE")


def add_bootstrap_support(tree, bootstrap_trees, root=True, branch_length_diff=None, outfolder='output', outgroups=['Chicken', 'Xenopus']):

    """
    Adds bootstrap support to a target phylogenetic tree, using bootstraped trees.
    The resulting tree with support values is written in `output/ete3_formatted_bootstrap_tree.nwk`

    Args:
        tree (str) : target tree file name
        boostrap_trees (list of str) : path  to the bootstrapped trees
        root (boolean, optional) : Whether to root the tree using Xenopus and Chicken.
        branch_length_diff (list of tuples, optional): species pair for which to extract branch
                                                       length difference from bootstrap trees

    Returns:
        (dict) : for each species pair (key) a list of branch length differences (value)
                 in bootstraped trees
    """

    diff_bl = {}


    #root all trees
    for treename in [tree] + bootstrap_trees:

        #reformat newick (sometimes R functions use comma for branch length values)
        with open(treename, 'r') as infile:
            mytree = infile.readlines()[0].strip()
            mytree = re.sub('0,([0-9])', '0.\\1', mytree)

        mytree = Tree(mytree)

        lca = mytree.get_common_ancestor(outgroups)
        if root:
            mytree.set_outgroup(lca)

        # if not mytree.check_monophyly(values=['Gar', 'Bowfin'], target_attr="name")[0]:
        #     print(mytree)

        if branch_length_diff:
            for species_pair in branch_length_diff:
                sp1, sp2 = species_pair
                ancestor = mytree.get_common_ancestor([sp1, sp2])
                branch_sp1 = mytree.get_distance(ancestor, sp1)
                branch_sp2 = mytree.get_distance(ancestor, sp2)
                diff = branch_sp1 - branch_sp2
                diff_bl[species_pair] = diff_bl.get(species_pair, [])
                diff_bl[species_pair].append(diff)

        mytree.write(outfile=treename+'_reformatted')

    mytrees = [Phylo.read(mytree+'_reformatted', "newick") for mytree in bootstrap_trees]

    target = Phylo.read(tree+'_reformatted', 'newick')
    tree_with_support = get_support(target, mytrees)

    # Phylo.draw_ascii(majority_tree)
    Phylo.write(tree_with_support, f'{outfolder}/bootstrap_tree.nwk', 'newick')

    #remove Bio.Phylo formatting that does not work with ete3 conventions...
    with open(f'{outfolder}/bootstrap_tree.nwk', 'r') as infile:
        species_tree = infile.readlines()[0]
    species_tree = species_tree.replace('):', ')')

    with open(f'{outfolder}/ete3_formatted_bootstrap_tree.nwk', 'w') as out:
        out.write(species_tree)

    return diff_bl
