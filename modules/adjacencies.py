"""
Module with fonctions to work with gene adjacencies data. 
"""

import itertools

def load_genomes(infile):
    
    """
    Load a multi-fasta of gene order and store all gene adjacencies.
    """
    adj = []
    d_seq = {}
    with open(infile,'r') as f:
        for line in f:
            if line[0] == '>' and adj:
                d_seq[sp] = adj
                sp = line[1:].split('_')[0]
                adj = []

            elif line[0] == '>':
                sp = line[1:].strip().split('_')[0]
            else:
                tmp_adj = line.strip().split()[:-1]

                if len(tmp_adj) > 1:
                    adj.append([int(i) for i in tmp_adj])

        #don't forget last species
        d_seq[sp] = adj
        
    return d_seq


def save_all_adj(d_seq,to_ignore=[]):
    """
    Save two dict of adjacencies per species one with adjecnecies in each sens.
    For instance, (-1,2) is equivalent to (-2,1).
    """
    adj_list = {}
    adj_list_rev = {}
    for sp in d_seq:
        adj_list[sp]=[]
        adj_list_rev[sp]=[]
        for chrom in d_seq[sp]:
            for i in range(len(chrom)-1):
                if chrom[i] != chrom[i+1]:
                    if sp == 'Tetraodon':
                        if ((int(chrom[i])),(int(chrom[i+1]))) not in to_ignore and ((int(chrom[i+1])),(int(chrom[i]))) not in to_ignore:
                            adj_list[sp].append((int(chrom[i]),int(chrom[i+1])))
                            adj_list_rev[sp].append((int(chrom[i+1])*-1,int(chrom[i])*-1))
                    else:
                        adj_list[sp].append((int(chrom[i]),int(chrom[i+1])))
                        adj_list_rev[sp].append((int(chrom[i+1])*-1,int(chrom[i])*-1))
    return adj_list,adj_list_rev

def indx_of_gene(ch,gene):
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


def adj_of_neighbors(ch, gene,gene_r,gene_l):
    
    """
    Search if an ortiented triplet gene_r,gene,gene_l in a duplicated species is in the same order
    in a non-duplicated species.
    """
    
    ok = False

    indx,indx_r,indx_l = None,None,None

    indx = indx_of_gene(ch,gene)

    indx_r = indx_of_gene(ch,gene_r)

    indx_l = indx_of_gene(ch,gene_l)

    if indx_l!=None and indx_r!= None and indx != None:

        adj = ()
        if indx_r < indx < indx_l and abs(indx_l-indx_r) < 10:
            adj = (ch[indx_r],ch[indx],ch[indx_l])

        elif indx_l < indx < indx_r and abs(indx_l-indx_r) < 10:
            adj = (ch[indx_l],ch[indx],ch[indx_r])


        if adj == (gene_l,gene,gene_r) or adj == (gene_r*-1,gene*-1,gene_l*-1):

            ok = True

    return ok


def make_matrix(adj_list, adj_list_rev, all_species, unRAND):
    """
    Make a binary matrix af gene adjacencies absence/presence. Will only be used for bootstrap.
    """
    ALL_SPECIES_PAIRS = list(itertools.combinations(all_species, 2))

    d_adj_all ={}
    for pair in ALL_SPECIES_PAIRS:
        for adj in adj_list[pair[0]]:
            if adj[0] != adj[1]: #no tandem duplication

                #adjacency in species 1
                
                #save unique name for adjacencies in the same orientation
                a1 = abs(adj[0])
                b1 = abs(adj[1])
                if a1<b1:
                    adj_to_save = (adj[0],adj[1])
                else:
                    adj_to_save = (adj[1]*-1,adj[0]*-1)

                #check orientation in species 2
                if adj in adj_list[pair[1]] or adj in adj_list_rev[pair[1]]:
                    if adj_to_save not in d_adj_all:
                        d_adj_all[adj_to_save] = []
                        
                    #save adjacency for both species if same orientation
                    d_adj_all[adj_to_save].append(pair)

                # if adjacency only in species 1 save it for species 1
                else:
                    if adj_to_save not in d_adj_all:
                        d_adj_all[adj_to_save] = []
                    if [pair[0]] not in d_adj_all[adj_to_save]:
                        d_adj_all[adj_to_save].append([pair[0]])

        #save adjacencies unique to species 2
        for adj in adj_list[pair[1]]:
            if adj not in adj_list[pair[0]] and adj not in adj_list_rev[pair[0]]:
                a1 = abs(adj[0])
                b1 = abs(adj[1])

                if a1<b1:
                    adj_to_save = (adj[0],adj[1])

                else:
                    adj_to_save = (adj[1]*-1,adj[0]*-1)

                if adj_to_save not in d_adj_all:
                        d_adj_all[adj_to_save] = []
                        
                if [pair[1]] not in d_adj_all[adj_to_save]:
                        d_adj_all[adj_to_save].append([pair[1]])

    #transform to a matrix and print some numbers
    all_adj = []
    gt1 = 0
    bt1 = 0
    t = 0
    gar = 3
    bowfin=2
    to_remove = []
    for fam in d_adj_all:
        vector = [0]*len(all_species)
        for hit in d_adj_all[fam]:

            for sp in hit:

                i = all_species.index(sp)
                vector[i]=1

        all_adj.append(vector)    
        if vector[gar]==1 and vector[bowfin]==0 and vector[0]==0 and vector[1]==0:
            for i in range(len(vector)):

                if i!=gar and i!=bowfin and vector[i]!=0 and i>gar:
                    gt1+=1
                    break


        #print families in support of Holostei monophyly
        if vector[gar] == 1 and vector[bowfin] == 1:
            ok=True
            for i in range(len(vector)):
                if i!=gar and i!=bowfin and vector[i]!=0 :
                    ok = False
                    break
            if ok:
                # print(fam)
                t+=1
                
                
        if vector[9] == 1:
            ok=True
            for i in range(len(vector)):
                if i!=9 and vector[i]!=0 :
                    ok = False
                    break
            if ok:
                if str(abs(fam[0])) in unRAND:
                    to_remove.append(fam)
                    all_adj = all_adj[:-1]

        if vector[gar] ==0 and vector[bowfin] == 1 and vector[0] == 0 and vector[1] == 0:
            ok = False
            for i in range(len(vector)):
                if i!=gar and i!=bowfin and vector[i]!=0 and i> gar:
                    bt1+=1
                    break

    print ('Number of total adjacencies: '+str(len(d_adj_all)))
    print( 'Number of derived adjacencies shared only by Bowfin and Gar: '+str(t) )
    print ('Number of derived adjacencies shared by Bowfin and Teleosts only: '+str(bt1))
    print ('Number of derived adjacencies shared by Gar and Teleosts only: '+str(gt1))
        
    return all_adj,to_remove