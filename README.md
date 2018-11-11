# CoSBI

I. Compile and install CoSBI

There is a make file in the package. To compile CoSBI, type
make 

II. Quickstart Guide


CoSBI –tf Tricluster_File –o Output_Directory –g min_g –s min_s –a alpha_value –b beta_value [ -af Annotation_File ]

-tf   Tricluster file path.
-o    Output directory.
-g    Min_g value.
-s    Min_s value.
-a    Alpha value.
-b    Beta value.
-af   Annotation file path. (optional)

Two example files are included in the “Data” folder. 
“T_cell_enhancers.tricluster” contains ChIP-Seq data for 213 known enhancers and equal number of random genomic windows in human CD4+ T cells. Each locus is represented by a 5K bp window centered on the center of the enhancer/random sequence. Each window is represented by 25 data points. Each data point is the ChIP-Seq tag counts of one of 39 histone modifications in human CD4+ T cells. 

“T_cell_enhancers.annotation” contains genomic coordinates for each of the 426 genomic loci in the “T_cell_enhancer.tricluster” file.
