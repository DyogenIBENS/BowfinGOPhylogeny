"""
Module with fonctions to work with gene adjacencies data.
"""

import itertools


def load_genomes(input_file):

    """
    Loads a multi-fasta of gene order and stores all gene adjacencies.

    Args:
        input_file (str): name of the input file

    Returns:
        (dict): for each species (key) a list of lists (value), giving, for each chromomsome (or
                scaffold), an ordered list of genes.
    """

    adj = []
    d_seq = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            if line[0] == '>' and adj:
                d_seq[sp] = adj
                sp = line[1:].split('_')[0]
                adj = []

            elif line[0] == '>':
                sp = line[1:].strip().split('_')[0]
            else:
                tmp_adj = line.strip().split()[:-1] #discard telomere end (would bias poor assembly)

                if len(tmp_adj) > 1:
                    adj.append([int(i) for i in tmp_adj])

        #don't forget last species
        d_seq[sp] = adj

    return d_seq


def save_all_adj(d_seq, to_ignore=None, sp_with_ign=None):

    """
    Saves two dict of adjacencies per species each with adjecencies in each direction.
    For instance, (-1,2) is equivalent to (-2,1) stored in the second dict.

    Args:
        d_seq (dict): for each species (key) a list of lists (value), giving,
                      for each chromomsome (or scaffold), an ordered list of genes.
        to_ignore (list of tuples): adjacency to filter in species sp_with_ign
        sp_with_ign (str): name of the species with adjacencies to filter out

    Returns:
        (dicts) : for each species (key), the list of adjacencies.
    """

    assert (to_ignore and sp_with_ign) or (not to_ignore and not sp_with_ign), "ambiguous arguments"

    adj_list, adj_list_rev = {}, {}

    for sp in d_seq:
        adj_list[sp], adj_list_rev[sp] = set(), set()

        for chrom in d_seq[sp]:

            for i in range(len(chrom)-1):

                if chrom[i] != chrom[i+1]: #ignore tandem duplications (can also be gene splits)

                    #since all tetraodon scaffolds are randomly assebled into a single contig
                    #discard those adjacencies on random contig that exist only in Tetraodon
                    if to_ignore and sp_with_ign and sp == sp_with_ign:
                        if ((int(chrom[i])), (int(chrom[i+1]))) not in to_ignore\
                        and ((int(chrom[i+1]))*-1, (int(chrom[i]))*-1) not in to_ignore:

                            adj_list[sp].add((int(chrom[i]), int(chrom[i+1])))
                            adj_list_rev[sp].add((int(chrom[i+1])*-1, int(chrom[i])*-1))

                    else:
                        adj_list[sp].add((int(chrom[i]), int(chrom[i+1])))
                        adj_list_rev[sp].add((int(chrom[i+1])*-1, int(chrom[i])*-1))

    return adj_list, adj_list_rev



def get_adj(adj_list, adj_list_rev, d_adj_all, pair):

    """
    Browses a pair of species and stores shared and unique gene adjacencies. Fills `d_adj_all`
    in-place.

    Args:
        adj_list, adj_list_rev (dicts): for each species (key), the list of adjacencies in each
                                        direction (value)
        d_adj_all (dict): resultung dict, filed-in place with the current species pair : for each
                          adjacency (key) the set of species that have it (value)
        pair (tuple of strings): species pair considered

    """

    for adj in adj_list[pair[0]]: #adjacency in species 1

        if adj[0] != adj[1]: #no tandem duplication

            #save unique name for adjacencies in the same orientation
            #FIXME do this at extraction to avoid having two dicts...
            gene1_id = abs(adj[0])
            gene2_id = abs(adj[1])
            if gene1_id < gene2_id:
                adj_to_save = (adj[0], adj[1])
            else:
                adj_to_save = (adj[1]*-1, adj[0]*-1)

            #check orientation in species 2
            if adj in adj_list[pair[1]] or adj in adj_list_rev[pair[1]]:
                if adj_to_save not in d_adj_all:
                    d_adj_all[adj_to_save] = set()

                #save adjacency for both species if same orientation
                d_adj_all[adj_to_save].add(pair[0])
                d_adj_all[adj_to_save].add(pair[1])

            # if adjacency only in species 1 save it for species 1
            else:
                if adj_to_save not in d_adj_all:
                    d_adj_all[adj_to_save] = set()
                d_adj_all[adj_to_save].add(pair[0])

    #also save adjacencies unique to species 2
    for adj in adj_list[pair[1]]:
        if adj not in adj_list[pair[0]] and adj not in adj_list_rev[pair[0]]:
            gene1_id = abs(adj[0])
            gene2_id = abs(adj[1])
            if gene1_id < gene2_id:
                adj_to_save = (adj[0], adj[1])
            else:
                adj_to_save = (adj[1]*-1, adj[0]*-1)

            if adj_to_save not in d_adj_all:
                d_adj_all[adj_to_save] = set()

            d_adj_all[adj_to_save].add(pair[1])



def dict_to_mat(d_adj_all, all_species, unrandom, sp_with_rm, gar=None, bowfin=None, chicken=None,
                xenopus=None):

    """
    Transforms the dict of adjacencies to a matrix and optionally printout some numbers on shared
    adjacencies.

    Args:
        d_adj_all (dict) : for each adjacency (key) the set of species that have it (value)
        all_species (list) : list of species to consider
        unrandom (list of tuple): list of adjacencies on random contig to discard
        sp_with_rm (string): name of species for which to discard adj on random contig
        gar, bowfin, chicken, xenopus : index of relevant species in `all_species` (to print counts)

    Returns:
        (list of list) : binary adjacency matrix, each sublist is a presence/absence vector, species
                         index in columns are given by theh `all_species` input argument.
        (list of tuples) : list of random contig adjacencies filtered out
    """

    all_adj = []
    to_remove = []
    count = False

    if bowfin is not None and gar is not None and chicken is not None and xenopus is not None:
        shared_gar_tel = 0
        shared_bow_tel = 0
        shared_gar_bow = 0
        count = True
        teleosts = [i for i in range(len(all_species)) if i not in [gar, chicken, bowfin, xenopus]]

    for fam in d_adj_all:
        vector = [0]*len(all_species)
        for sp in d_adj_all[fam]:

            i = all_species.index(sp)
            vector[i] = 1

        all_adj.append(vector)


        #filter out tetraodon adj
        idx = all_species.index(sp_with_rm)
        if vector[idx] == 1:
            ok = True
            for i, val in enumerate(vector):
                if i != idx and val != 0:
                    ok = False
                    break
            if ok:
                if str(abs(fam[0])) in unrandom:
                    to_remove.append(fam)
                    all_adj = all_adj[:-1]

        #count shared derived adjacencies for the parsimony analysis
        if count:

            if vector[gar] == 1 and vector[bowfin] == 0 and vector[chicken] == 0\
               and vector[xenopus] == 0:

                if count:
                    for i, val in enumerate(vector):

                        if i != gar and i != bowfin and val != 0 and i > gar:
                            shared_gar_tel += 1
                            break

            #families in support of Holostei monophyly
            if vector[gar] == 1 and vector[bowfin] == 1:
                ok = True
                for i, val in enumerate(vector):
                    if i != gar and i != bowfin and val != 0:
                        ok = False
                        break
                if ok:
                    shared_gar_bow += 1

            if vector[gar] == 0 and vector[bowfin] == 1 and vector[chicken] == 0\
               and vector[xenopus] == 0:

                ok = False
                for i, val in enumerate(vector):
                    if i != gar and i != bowfin and val != 0 and i in teleosts:
                        shared_bow_tel += 1
                        break

    if count:

        print(f'Number of total adjacencies: '+str(len(d_adj_all)))
        print(f'Number of derived adjacencies shared only by Bowfin and Gar: {shared_gar_bow}')
        print(f'Number of derived adjacencies shared by Bowfin and Teleosts only: {shared_bow_tel}')
        print(f'Number of derived adjacencies shared by Gar and Teleosts only: {shared_gar_tel}')

    return all_adj, to_remove


def make_binary_adj_matrix(adj_list, adj_list_rev, all_species, unrandom=None, sp_with_rm=None,
                           print_count=True):
    """
    Makes a binary matrix af gene adjacencies absence/presence. Will be used for bootstrap.

    Args:
        adj_list, adj_list_rev (dicts): for each species (key), the list of adjacencies in each
                                        direction (value)
        all_species (list) : list of species to consider
        unrandom (list) : list of families on random contig
        sp_with_rm (str) : species with random contig (currently support for only one species)
        print_count (boolean) : whether to print stats on shared adjacencies

    Returns:
        (list of list) : binary adjacency matrix, each sublist is a presence/absence vector, species
                         index in columns are given by theh `all_species` input argument.
        (list of tuples) : list of random contig adjacencies filtered out
    """
    all_species_pairs = list(itertools.combinations(all_species, 2))

    d_adj_all = {}
    for pair in all_species_pairs:
        get_adj(adj_list, adj_list_rev, d_adj_all, pair)

    #transform to a matrix and print some numbers
    if print_count:
        gar = all_species.index("Lepisosteus oculatus")
        bowfin = all_species.index("Amia calva")
        chicken = all_species.index("Gallus gallus")
        xenopus = all_species.index("Xenopus tropicalis")
        all_adj, to_remove = dict_to_mat(d_adj_all, all_species, unrandom, sp_with_rm, gar, bowfin,
                                         chicken, xenopus)
    #else only make matrix
    else:
        dict_to_mat(d_adj_all, all_species, unrandom, sp_with_rm)

    return all_adj, to_remove
