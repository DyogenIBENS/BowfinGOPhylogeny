import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style("white")



def indx_of_gene(ch, gene):
    """
    Return the relative position of a gene on a chromosome.
    None if absent.
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
    Search if an ortiented triplet gene_r, gene, gene_l in a duplicated species is in the same order
    in a non-duplicated species.
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


def correct_frac_bias(adj_list, adj_list_rev, all_speciesecies, non_dup):
    
    """
    Try to correct adjacencies absent of a duplicated species because of post-WGD fractionation.
    """
    
    for SP in NON_D:
        for ADJ in adj_list[SP]:

            for sp in [i for i in all_speciesecies if i not in non_dup]:
                all_ch_0 = []
                all_ch_1 = []
                if ADJ not in adj_list[sp] and ADJ not in adj_list_rev[sp]:

                    ok = False

                    #find neighbors of gene_0 and neighbors of gene_1
                    gene_0, gene_1 = ADJ
                    ch_ind =0
                    for ch in d_seq[sp]:
                        
                        #neighbors of gene_0
                        indx = indx_of_gene(ch, gene_0)
                        if indx != 0 and indx != None and indx < len(ch)-1:
                            all_ch_0.append(ch)
                            ch1 = ch_ind
                            gene_0_l, gene_0_r = ch[indx-1], ch[indx+1]
                            found_0 = True

                        #neighbors of gene 1
                        indx = indx_of_gene(ch, gene_1)
                        if indx != 0 and indx != None and indx < len(ch)-1:
                            all_ch_1.append(ch)
                            ch2 = ch_ind
                            gene_1_l, gene_1_r = ch[indx-1], ch[indx+1]
                            found_1 = True

                        ch_ind +=1

                    if len(all_ch_0) ==1 and len(all_ch_1)==1 and found_1 and found_0 and ch1!=ch2:

                        #search if 'adj' gene_0_l, gene_0, gene_0_r exists in non dup
                        for ch in d_seq[SP]:

                            if adj_of_neighbors(ch, gene_0, gene_0_r, gene_0_l):

                                    if adj_of_neighbors(ch, gene_1, gene_1_r, gene_1_l):

                                        ok = True
                                        
                    if ok:

                        adj_list[sp].append(ADJ)

                        adj_list_rev[sp].append((ADJ[1]*-1, ADJ[0]*-1))
    return



def make_distance_matrix(adj_list, adj_list_rev, sp_dict, all_species):
    
    """
    Make a distance matrix and draw it.
    """
    all_species_pairs = list(itertools.combinations(all_species, 2))
    mat = np.zeros((len(all_species), len(all_species)))
    d_adj_all = {}
    for pair in all_species_pairs:
        i = all_species.index(pair[0])
        j = all_species.index(pair[1])
        for adj in adj_list[pair[0]]:
            if adj[0] != adj[1]:

                #check orientation
                if adj in adj_list[pair[1]] or adj in adj_list_rev[pair[1]]:
                    mat[i, j] += 1
                    mat[j, i] += 1

        mat[i, j] = 1 - mat[i, j]/float(min(len(adj_list[pair[0]]), len(adj_list[pair[1]])))
        mat[j, i] = mat[i, j]
    

    DF_data = pd.DataFrame(mat, 
                            columns = [sp_dict[i.split()[0]] for i in all_species],
                            index = [sp_dict[i.split()[0]] for i in all_species])


    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(8, 8))

    # Draw the heatmap 
    sns.heatmap(DF_data, cmap=sns.cm.rocket, square=True, vmin=0, vmax=1)

    plt.tight_layout()
    plt.savefig('output/distance_matrix.svg', dpi=300)

    plt.show()
    
    #save the matrix
    np.savetxt('output/dist_mat.txt', mat, header=' '.join([sp_dict[i.split()[0]].replace(' ', '_') for i in all_species]))
    
    #format it for R ape package
    with open('output/dist_mat.txt','r') as f, open('output/dist_mat','w') as fw:
        res = ''
        for i, line in enumerate(f):
            if i == 0:
                res += line[2:]
            else:
                res += line
        fw.write(res)
    return
