rm(list=ls())
load("Monocytes/data/expression_data_monocytes_and_bcells_LYZ_region.RData")


# Extract data using only genes in LYZ region identified with FDR<0.05 ----------------------

ls()

dim(expr_ifn) # Currently 1073 genes

# Get the 0.05 FDR threshold values 
fdr.thresh = mat_fdr_thres[,1]
names(fdr.thresh)[2] = 'Unstim' # Make capital letter to match other data frames

# Get the relevant data
df_mono = df_hits[df_hits$Condition!='B-cells',] # Do not use B-cell data
df_mono = df_mono[df_mono$PPI>=fdr.thresh[df_mono$Condition],] # Only above 0.05 FDR threshold

# Get the relevant gene names
genes_id = unique(df_mono$Gene_symbol)
length(genes_id) # now 381 genes

# Reduce dimension of data sets to only these genes
expr_ifn_mono = expr_ifn[,genes_id]
expr_lps2_mono = expr_lps2[,genes_id]
expr_lps24_mono = expr_lps24[,genes_id]
expr_unstim_mono = expr_unstim[,genes_id]

# Save objects
expr_all = list(expr_ifn_mono,expr_lps2_mono,expr_lps24_mono,expr_unstim_mono)
save(expr_all,genes_id, file="Monocytes/data/expression_data_monocytes_LYZ_region_FDR_5.RData")







