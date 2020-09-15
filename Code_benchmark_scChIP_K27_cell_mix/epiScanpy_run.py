import anndata as ad
import episcanpy.api as epi
import scanpy as sc
import numpy as np
import pandas as pd
import re
import pyranges as pr

sc.set_figure_params(scanpy=True, dpi=80, dpi_save=250,
                     frameon=True, vector_friendly=True,
                     color_map="YlGnBu", format='pdf', transparent=False,
                     ipython_format='png2x')

DATADIR = './EpiScanpy/'

HBCx22 = epi.read_text("Matrices/HBCx_22.tsv", delimiter="\t",first_column_names="regions")
Jurkat = epi.read_text("Matrices/Jurkat.tsv", delimiter="\t",first_column_names="regions")
Ramos = epi.read_text("Matrices/Ramos.tsv", delimiter="\t",first_column_names="regions")
MM468 = epi.read_text("Matrices/MM468.tsv", delimiter="\t",first_column_names="regions")

MM468.var_names = "MM468_" + MM468.var_names
HBCx22.var_names = "HBCx22_" + HBCx22.var_names
Jurkat.var_names = "Jurkat_" + Jurkat.var_names
Ramos.var_names = "Ramos_" + Ramos.var_names

Ramos = Ramos.T
HBCx22 = HBCx22.T
Jurkat = Jurkat.T
MM468 = MM468.T

adata = MM468.concatenate(HBCx22, Jurkat, Ramos, join="inner", index_unique=None)

cells_passing_QC = pd.read_csv("cells_passing_QC.txt",names=["cells"], header=None) 

l=list()
for index, row in cells_passing_QC.iterrows():
    l.append(row['cells'])
    print(row['cells'])

exclude_regions = pr.read_bed("Additional_Files/MM468_identified_CNA.bed")
regions =  pd.DataFrame(adata.var)
regions= pd.DataFrame(regions.index)
regions[regions.columns[0]]

chrom = [re.sub(":.*", "",x) for x in regions[regions.columns[0]] ]
start  = [re.sub("-.*", "",x) for x in regions[regions.columns[0]] ]
start  = [re.sub(".*:", "",x) for x in start ]
end = [re.sub(".*-", "",x) for x in regions[regions.columns[0]] ]

chrom = np.array(chrom)
start = np.array(start)
end = np.array(end)

region_range = pr.PyRanges(chromosomes=chrom,starts=start,ends=end)

to_excl = region_range.overlap(other=exclude_regions)
to_excl = to_excl.df["Chromosome"].astype(str) +":" + to_excl.df["Start"].astype(str) + "-" + to_excl.df["End"].astype(str)

adata.obs_names = [re.sub("HBCx22", "HBCx_22",x) for x in adata.obs_names]

adata = adata[adata.obs_names.isin(l),]
adata = adata[:,~adata.var_names.astype(str).isin(pd.Series.tolist(to_excl))]

epi.pp.binarize(adata, copy=False)
epi.pp.coverage_cells(adata, binary=True, bins=100)
epi.pp.commonness_features(adata, binary=True)
epi.pp.filter_features(adata, min_cells=1)
epi.pp.commonness_features(adata, binary=True, log=True)
adata50 = epi.pp.select_var_feature(adata, nb_features=50000, copy=True)
adata50
epi.pp.regress_out(adata50, "nb_features")

epi.pp.lazy(adata50)

epi.pp.pca(adata50, n_comps=50, svd_solver='arpack')

reduced_dim = adata50.obsm.get("X_pca").T
reduced_dim = pd.DataFrame(reduced_dim)
reduced_dim.columns = adata50.obs_names
reduced_dim.to_csv(DATADIR + "reduced_dim_EpiScanpy.csv")




