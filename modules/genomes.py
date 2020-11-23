"""
Module with fonctions to work with gene families and gene order data. 
"""

import bz2

def extract_all_genes(genefile):
    """
    Get all genes of a species.


    """
    sp_genes = set()
    genes_on_random = set()
    with bz2.open(genefile, 'rt') as f:
        for line in f:
            line = line.strip().split('\t')
            sp_genes.add(line[-1])
            if 'random' in line[0]: #not pretty
                genes_on_random.add(line[-1])
    return sp_genes, genes_on_random



def filter_families(genes_allsp, non_dup_sp, input_file, output_file):
    """
    For gene families defined by the Euteleostomi ancestor gene, extract gene families with exactly
    one member in non-duplicated species and one or two in duplicated species.
    Intermediary oupputs contain these families for further inspection.
    """
    
    count_genes = {}
    tmp_count_genes = {}
    markers = 0
    with bz2.open(input_file, "rt") as infile, open(output_file,'w') as outfile,\
         open('output/Families_1-to-1','w') as fw1, open('output/Families_1-to-2','w') as fw2:
        
        #go through all families
        for line in infile:
            
            ok = True
            fam1to1 = True
            fam1to2 = True
            
            at_least_one = False
            
            #check gene content for each species
            for sp in genes_allsp:

                familie_genes = set(line[1:].split())
                
                tmp_count_genes[sp] = genes_allsp[sp].intersection(familie_genes)

                if sp.replace(".", " ") in non_dup_sp:
                    if len(genes_allsp[sp].intersection(familie_genes)) != 1 :
                        ok = False
                    
                else:
                    if (len(genes_allsp[sp].intersection(familie_genes)) != 1 and
                        len(genes_allsp[sp].intersection(familie_genes)) != 2 ):

                        ok = False
                    
                    if len(genes_allsp[sp].intersection(familie_genes)) != 1 :
                        fam1to1 = False
                    
                    if len(genes_allsp[sp].intersection(familie_genes)) != 2 :
                        fam1to2 = False
                        
                if len(genes_allsp[sp].intersection(familie_genes)) > 0:
                    at_least_one = True
                    
            #write filtered families and a summary of strict 1:1 and strict 1:2 families
            if ok:
                markers += 1
                outfile.write(line)
                for sp in tmp_count_genes:

                    #summary: family name (ancestor gene), modern species and gene
                    for gene in list(tmp_count_genes[sp]):
                        if fam1to1:
                            fw1.write(line.split()[0]+'\t'+gene+'_'+sp+'\tFamily_1-to-1\n')
                        elif fam1to2:
                            fw2.write(line.split()[0]+'\t'+gene+'_'+sp+'\tFamily_1-to-2\n')
                            
                    if sp not in count_genes:
                        count_genes[sp]=0
                    count_genes[sp] += len(tmp_count_genes[sp])

    #print final number of genes per species
    print(f"Number of marker genes selected : {markers}.")
    print("Number of genes in each species:", count_genes)
    return


def read(infile, genes_on_random):
    """
    Loads a families file.
    """
    d = {}
    unRAND = []
    with open(infile,'r') as f:
        for i, line in enumerate(f):
            line = line.strip().split()

            for gene in line[1:]:
                d[gene] = str(i+1)
                
                if gene in genes_on_random:
                    unRAND.append(str(i+1))
                
    return d, unRAND


def write_genomes(name_families, genefile, all_species, outfile):
    """
    Write genomes in a multiple fasta format with genes ordered along chromsomes, with chromosome
    ends indicated by `$`. Returns a dictionary storing all gene adjacencies for each species.


    """

    #start by emptying the outfile in case it exists, since we will append to it in a loop.
    with open(outfile, 'w'): pass

    d = {}
    for sp in all_species:
        i = 0
        #Load genes file
        fisha = genefile % (sp.replace(' ', '.'))
        with bz2.open(fisha,'rt') as f, open(outfile,'a') as fw:
            prev=''
            i+=1
            #write species header
            fw.write('>'+sp.split('.')[0]+'_'+str(i)+'\n')

            #sort genes by chromomomes and by position
            a = [line.strip().split('\t') for line in f]
            a.sort(key=lambda x:( x[0],int(x[1])))

            first_gene = ''
            
            for line in a:

                chrom = line[0]
                sens = line[3]
                
                #save direction as '-' if on reverse strand
                if sens =='-1':
                    sens = '-'
                else:
                    sens = ''

                gene = line[4]

                #Write gene if in families to keep
                if gene in name_families:

                    nb = sens + name_families[gene]# convert gene to a gene family ID

                    if chrom == prev or prev == '':
                        fw.write(nb+' ')
                        
                    else:
                        fw.write('$\n'+nb+' ') #new chromosome


                    prev = chrom
                    prev_gene = gene

            fw.write('$\n')

