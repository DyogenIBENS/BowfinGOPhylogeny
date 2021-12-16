"""

Module with function to build a distance matrix from adjacencies data.

"""

import itertools

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_style("white")


def indx_of_gene(ch, gene):

    """
    Returns the relative position of a gene on a chromosome.
    None if absent.

    Args:
        ch (list of int): chromosome as a list of genes
        gene (int): gene ID

    Returns:
        (int): index of gene on chromosomes.
    """

    indx = None

    if gene in ch or -1*gene in ch:
        if gene in ch:
            indx = ch.index(gene)

        else:
            indx = ch.index(gene*-1)

    return indx


def adj_of_neighbors(ch, gene, gene_r, gene_l):

    """
    Helper function to search for relaxed adjacencies between a non-duplicated species and a
    duplicated species. Checks whether an oriented triplet gene_r, gene, gene_l in a duplicated
    species is in the same order in a non-duplicated specied species, and is less than 10 genes
    apart.

    Args:
        ch (list of int): considered chromosome of the non-duplicated species, as a list of genes
        gene, gene_r, gene_l (str): genes of the oriented triplet to search

    Returns:
        (int): whether the triplet exists in the same order in the non-duplicated species
    """

    ok = False

    indx, indx_r, indx_l = None, None, None

    indx = indx_of_gene(ch, gene)

    indx_r = indx_of_gene(ch, gene_r)

    indx_l = indx_of_gene(ch, gene_l)

    if indx_l is not None and indx_r is not None and indx is not None:

        adj = ()
        if indx_r < indx < indx_l and abs(indx_l-indx_r) < 10:
            adj = (ch[indx_r], ch[indx], ch[indx_l])

        elif indx_l < indx < indx_r and abs(indx_l-indx_r) < 10:
            adj = (ch[indx_l], ch[indx], ch[indx_r])


        if adj in ((gene_l, gene, gene_r), (gene_r * -1, gene * -1, gene_l * -1)):

            ok = True

    return ok


def correct_frac_bias(adj_list, adj_list_rev, all_species, d_seq, non_dup):

    """
    Searches for adjacencies absent in a duplicated species because of post-WGD fractionation.
    For each adj present in a non-duplicated species, if gene losses are responsible for breaking
    the adjacency, rather than rearrangements, then the adj is counted as present in the duplicated
    species.

    More precisely, for an adjacency B-C in a non-duplicated species which corresponds
    to I-B-K and X-C-Y triplets in a duplicated species, if the triplets are found in the same order
    and at a distance < 10 genes in the duplicated species, then the adjaency is considered as
    broken by gene losses.

    The adjacency dicts `adj_list` and `adj_list_rev` are updated in-place.

    Args:
        adj_list, adj_list_rev (dicts): for each species (key), the list of adjacencies in each
                                        direction (value)
        all_species (list) : list of species to consider
        d_seq (dict): for each species (key) a list of lists (value), giving,
                      for each chromomsome (or scaffold), an ordered list of genes.
        non_dup (list): list of non-duplicated species
    """

    for non_dup_sp in non_dup:
        for adj in adj_list[non_dup_sp]:

            for dup_sp in [i for i in all_species if i not in non_dup]:
                all_ch_0 = []
                all_ch_1 = []
                if adj not in adj_list[dup_sp] and adj not in adj_list_rev[dup_sp]:

                    ok = False

                    #find neighbors of gene_0 and neighbors of gene_1
                    gene_0, gene_1 = adj
                    ch_ind = 0
                    for ch in d_seq[dup_sp]:

                        #neighbors of gene_0
                        indx = indx_of_gene(ch, gene_0)
                        if indx != 0 and indx is not None and indx < len(ch)-1:
                            all_ch_0.append(ch)
                            ch1 = ch_ind
                            gene_0_l, gene_0_r = ch[indx-1], ch[indx+1]
                            found_0 = True

                        #neighbors of gene 1
                        indx = indx_of_gene(ch, gene_1)
                        if indx != 0 and indx is not None and indx < len(ch)-1:
                            all_ch_1.append(ch)
                            ch2 = ch_ind
                            gene_1_l, gene_1_r = ch[indx-1], ch[indx+1]
                            found_1 = True

                        ch_ind += 1

                    if len(all_ch_0) == 1 and len(all_ch_1) == 1 and found_1 and found_0\
                       and ch1 != ch2:

                        #search if 'adj' gene_0_l, gene_0, gene_0_r exists in non dup
                        for ch in d_seq[non_dup_sp]:

                            if adj_of_neighbors(ch, gene_0, gene_0_r, gene_0_l):

                                if adj_of_neighbors(ch, gene_1, gene_1_r, gene_1_l):

                                    ok = True

                    if ok:

                        adj_list[dup_sp].add(adj)

                        adj_list_rev[dup_sp].add((adj[1]*-1, adj[0]*-1))


def make_distance_matrix(adj_list, adj_list_rev, sp_dict, all_species, show=True, outfolder='output'):

    """
    Makes a distance matrix from adjacencies data using a normalized breakpoint distance
    (1 - proportion of shared adjacencies) and draws it.

    Args:
        adj_list, adj_list_rev (dicts): for each species (key), the list of adjacencies in each
                                        direction (value)
        sp_dict (dict) : dict giving latin (key) to common (value) species name (for plot)
        all_species (list) : list of species to consider
        show (boolean): whether or not to show the plot
    """

    all_species_pairs = list(itertools.combinations(all_species, 2))
    mat = np.zeros((len(all_species), len(all_species)))
    for pair in all_species_pairs:
        i = all_species.index(pair[0])
        j = all_species.index(pair[1])
        for adj in adj_list[pair[0]]:
            if adj[0] != adj[1]:

                #check if exist in any of the two equivalent orientations
                if adj in adj_list[pair[1]] or adj in adj_list_rev[pair[1]]:
                    mat[i, j] += 1

        mat[i, j] = 1 - mat[i, j]/float(min(len(adj_list[pair[0]]), len(adj_list[pair[1]])))
        mat[j, i] = mat[i, j]


    df_data = pd.DataFrame(mat,
                           columns=[sp_dict[i.split()[0]] for i in all_species],
                           index=[sp_dict[i.split()[0]] for i in all_species])


    #set up the matplotlib figure
    _, _ = plt.subplots(figsize=(8, 8))

    #draw the heatmap
    sns.heatmap(df_data, cmap=sns.cm.rocket, square=True, vmin=0, vmax=1)

    plt.tight_layout()
    plt.savefig(f'{outfolder}/distance_matrix.svg', dpi=300)

    if show:

        plt.show()

    plt.close('all')

    #save the matrix as text
    np.savetxt(f'{outfolder}/dist_mat.txt', mat,
               header=' '.join([sp_dict[i.split()[0]].replace(' ', '_') for i in all_species]))

    #format it for R ape package
    with open(f'{outfolder}/dist_mat.txt', 'r') as infile, open(f'{outfolder}/dist_mat', 'w') as out:
        res = ''
        for i, line in enumerate(infile):
            if i == 0:
                res += line[2:]
            else:
                res += line
        out.write(res)
