## SOURCE DATA for Fig. 1c and 1d (Bowfin vs Spotted gar and Bowfin vs Medaka)

- `data/supplements1/ancGenes.Neopterygii.list.bz2` (orthology data): All neopterygii gene families extracted from TreeBeST gene trees, as a 2-columns file (column 1: family name, column 2: all genes in the family). TreeBeST gene trees are in `data/supplements1/gene_trees.nhx.bz2`.

- `data/gene_order_phylogeny/genes.*.list.bz2` (gene coordinates files): 5-columns coordinates files (column 1: chromosome, column 2: start, column 3: stop, column 4 : orientation, column 5: gene name)

## SOURCE DATA for Fig. 1e and 1f (gene-order NJ phylogeny and gene-order parsimony phylogeny)

- `data/gene_order_phylogeny/ancGenes.Euteleostomi.list.bz2` (orthology data): All Euteleostomi gene families extracted from TreeBeST gene trees, as a 2-columns file (column 1: family name, column 2: genes in family). TreeBeST gene trees are in `data/supplements1/gene_trees.nhx.bz2`.

- `data/gene_order_phylogeny/genes.*.list.bz2` (gene coordinates files, as above).

See the notebook `gene_order_phylogeny.ipynb` for implementation details.

## SOURCE DATA for Extended Data Fig. 2

- As for Fig. 1c and Fig. 1d, for comparisons between bowfin, gar and medaka.

- `data/gene_order_phylogeny/ancGenes.Euteleostomi.list.bz2` (orthology data, same as above): for comparisons with chicken.

## SOURCE DATA for Supp Fig. 5

- `data/supplements2/orthologs_bichir_*.txt.bz2` (orthology data): pairwise reciprocal best blast hits for bichir vs bowfin and bichir vs spotted gar, 3-column files (column 1: orthology id, column 2: bichir gene, column 3: orthologous gene).

- `data/supplements2/genes.*.list.bz2` (gene coordinates files, as above): gene coordinates files for bichir, bowfin and spotted gar.