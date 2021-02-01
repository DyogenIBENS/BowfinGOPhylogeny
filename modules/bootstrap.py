import os

import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style("white")

from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus,get_support
from ete3 import Tree, NodeStyle, TreeStyle,TextFace, AttrFace,faces

import random
random.seed(1234)


def bootstrap_matrix(all_adj, all_species, sp_dict):
    """
    Bootstrap the binary matrix 100 times and transform to 100 distance matrix.
    """
    A = np.array(all_adj, dtype=int)
    for k in range(100):

        dmat = np.zeros((len(all_species),len(all_species)))
        all_vect = range(int(len(all_adj)))
        if (k+1)%10 == 0:
            print(f'{k+1} bootstrap replicates done')
        boot_indices = [random.choice(all_vect) for _ in range(int(len(all_adj)))]
        boot_replicate = A[boot_indices,:]

        for i,j in list(itertools.combinations(all_species, 2)):
            i = all_species.index(i)
            j = all_species.index(j)


            nb_sim = len(np.logical_and(boot_replicate[:,i],boot_replicate[:,j]).nonzero()[0])

            max1 = len([a for a in boot_replicate[:,i] if a==1])
            max2 = len([a for a in boot_replicate[:,j] if a==1])
            dmat[i,j] = 1 - nb_sim/float(min(max1,max2))
            dmat[j,i] = dmat[i,j]
        try:
            os.mkdir('output/bootstrap/')
        except OSError:
            pass
        
        #save matrix
        np.savetxt('output/bootstrap/dist_mat.txt', dmat, header=' '.join([sp_dict[i.split()[0]].replace(' ', '_') for i in all_species]))
        
        #Format matrix for ape R package
        with open('output/bootstrap/dist_mat.txt','r') as f, open('output/bootstrap/dist_mat_'+str(k)+'.txt', 'w') as fw:
            res = ''
            for i, line in enumerate(f):
                if i == 0:
                    res += line[2:]
                else:
                    res+=line
            fw.write(res)
    print(f'{k+1} bootstrap replicates done')


def add_bootstrap_support(tree, bootstrap_trees):
    """
    Add bootstrap support to the NJ tree.
    """
    diff_bl = []


    #root all trees
    for treename in [tree] + bootstrap_trees:
        t = Tree(treename)
        lca = t.get_common_ancestor(['Chicken', 'Xenopus'])
        t.set_outgroup(lca)
        if not t.check_monophyly(values=['Gar', 'Bowfin'], target_attr="name")[0]:
            print(t)
            
        else:
            ac = t.get_common_ancestor(['Gar', 'Bowfin'])
            b_gar = t.get_distance(ac, 'Gar')
            b_bow = t.get_distance(ac, 'Bowfin')
            diff_bl.append(b_bow-b_gar)

        t.write(outfile=treename)

    mytrees = [Phylo.read(t, "newick") for t in bootstrap_trees]


    target = Phylo.read(tree, 'newick')
    majority_tree = get_support(target, mytrees)

    # Phylo.draw_ascii(majority_tree)
    Phylo.write(majority_tree, 'output/bootstrap_consensus.nwk', 'newick')

    #remove Bio.Phylo formatting that does not work with ete3 conventions...
    with open('output/bootstrap_consensus.nwk','r') as f:
        t = f.readlines()[0]
    t = t.replace('):',')')

    with open('output/bootstrap_consensus.nwk','w') as f:
        f.write(t)
        
    # return diff_bl
    return