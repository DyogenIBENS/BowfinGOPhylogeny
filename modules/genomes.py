"""
Module with fonctions to work with gene families and gene order data.
"""

import bz2

def extract_all_genes(genefile):

    """
    Get all genes of a species.

    Arg:
        genefile (str): name of the input gene coordinate file in DYOGEN format.

    Returns:
        (set) : full genes set
        (set) : subset of genes on random contig

    """

    sp_genes = set()
    genes_on_random = set()
    with bz2.open(genefile, 'rt') as infile:
        for line in infile:
            line = line.strip().split('\t')
            sp_genes.add(line[-1])
            if 'random' in line[0]: #working here but not pretty nor general
                genes_on_random.add(line[-1])
    return sp_genes, genes_on_random



def filter_families(genes_allsp, non_dup_sp, input_file, output_file):

    """
    For gene families (here defined by the Euteleostomi ancestral gene), extract gene families with
    exactly one member in non-duplicated species and one or two in duplicated species.
    Intermediary output files with these families are written for further inspection.

    Args:
        genes_allsp (dict): for each species (key) its genes as a set (value)
        non_dup_sp (list of str): name of non-duplicated species
        input_file (str): name of the ancgene family file (DYOGEN format).
        output_file (str): path to the output file where to write the subset of gene families
    """

    count_genes = {}
    tmp_count_genes = {}
    markers = 0
    with bz2.open(input_file, "rt") as infile, open(output_file, 'w') as outfile,\
         open('output/Families_1-to-1.tsv', 'w') as fw1,\
         open('output/Families_1-to-2.tsv', 'w') as fw2,\
         open('Families_all.tsv', 'w') as fw3:

        #go through all families
        for line in infile:

            ok = True
            fam1to1 = True
            fam1to2 = True

            #check gene content for each species
            for sp in genes_allsp:

                familie_genes = set(line[1:].split())

                tmp_count_genes[sp] = genes_allsp[sp].intersection(familie_genes)

                if sp.replace(".", " ") in non_dup_sp:
                    if len(genes_allsp[sp].intersection(familie_genes)) != 1:
                        ok = False

                else:
                    if len(genes_allsp[sp].intersection(familie_genes)) != 1 and\
                       len(genes_allsp[sp].intersection(familie_genes)) != 2:

                        ok = False

                    if len(genes_allsp[sp].intersection(familie_genes)) != 1:
                        fam1to1 = False

                    if len(genes_allsp[sp].intersection(familie_genes)) != 2:
                        fam1to2 = False

            #write filtered families and a summary of strict 1:1 and strict 1:2 families
            if ok:
                markers += 1
                outfile.write(line)
                for sp in tmp_count_genes:

                    #summary: family name (ancestor gene), modern species and gene
                    for gene in list(tmp_count_genes[sp]):
                        fw3.write(line.split()[0]+'\t'+gene+'_'+sp+'\tAll_families\n')
                        if fam1to1:
                            fw1.write(line.split()[0]+'\t'+gene+'_'+sp+'\tFamily_1-to-1\n')
                        elif fam1to2:
                            fw2.write(line.split()[0]+'\t'+gene+'_'+sp+'\tFamily_1-to-2\n')

                    if sp not in count_genes:
                        count_genes[sp] = 0
                    count_genes[sp] += len(tmp_count_genes[sp])

    #print final number of genes per species
    print(f"Number of marker gene families selected: {markers}.")
    print("Number of corresponding genes in each species:", count_genes)


def read(input_file, genes_on_random):

    """
    Loads a families file.

    Args:
        input_file (str): family file to load
        genes_on_random (set): genes on a random contig

    Returns:
        (dict): families defined by ancgene name (key) and a corresponding unique ID (value)
        (list): unique ID of families with a member on random contig
    """

    families = {}
    un_random = []
    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):
            line = line.strip().split()

            for gene in line[1:]:
                families[gene] = str(i+1)

                if gene in genes_on_random:
                    un_random.append(str(i+1))

    return families, un_random


def write_genomes(name_families, genefile, all_species, outfile):

    """
    Write genomes in a multiple fasta format with genes ordered along chromsomes, with chromosome
    ends indicated by `$`. Also returns a dictionary storing all gene adjacencies for each species.

    Args:
        name_families (dict) : selected ancgenes (key) and their unique ID (value)
        genefile (str) : path to gene coordinates files in DYOGEN format
        all_species (list of str): name of species to consider, used to find corresponding
                                   genefiles.
        outfile (str): path to the output multi-fasta genome files.
    """

    #start by emptying the outfile in case it exists, since we will append to it in a loop.
    open(outfile, 'w').close()

    for sp in all_species:
        i = 0
        #Load genes file
        fisha = genefile % (sp.replace(' ', '.'))
        with bz2.open(fisha, 'rt') as infile, open(outfile, 'a') as out:
            prev = ''
            i += 1
            #write species header
            out.write('>'+sp.split('.')[0]+'_'+str(i)+'\n')

            #sort genes by chromomomes and by position
            data = [line.strip().split('\t') for line in infile]
            data.sort(key=lambda x: (x[0], int(x[1])))

            for line in data:

                chrom = line[0]
                sens = line[3]

                #save direction as '-' if on reverse strand
                if sens == '-1':
                    sens = '-'
                else:
                    sens = ''

                gene = line[4]

                #Write gene if in families to keep
                if gene in name_families:

                    nb = sens + name_families[gene] # convert gene to a gene family ID

                    if prev in (chrom, ''):
                        out.write(nb+' ')

                    else:
                        out.write('$\n'+nb+' ') #new chromosome

                    prev = chrom

            out.write('$\n')
