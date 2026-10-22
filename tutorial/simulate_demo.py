"""
Tutorial step 1: simulate the demo dataset with known gene expression programs.
Same generative model and parameters as the paper's simulations.py (Splatter-style), scaled down:
default 20,000 cells x 10,000 genes, 20 cell-identity groups + 1 activity program (groups 1-4,
35% of their cells, usage 0.1-0.7, last 1,000 genes) = 21 true programs.
Usage: python simulate_demo.py <out_dir> [n_cells] [n_genes] [n_groups] [seed]
Outputs: <out_dir>/demo_counts.h5ad (raw counts, sparse CSR), <out_dir>/ground_truth_programs.tsv (gene, program),
         <out_dir>/cellparams.csv (group, has_program, program_usage per cell)
"""
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

out_dir = sys.argv[1]
ncells = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
ngenes = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
ngroups = int(sys.argv[4]) if len(sys.argv) > 4 else 20
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 9458
os.makedirs(out_dir, exist_ok=True)

# parameters identical to simulations.py
libloc, libscale = 7.64, 0.78
mean_rate, mean_shape = 7.68, 0.34
expoutprob, expoutloc, expoutscale = 0.00286, 6.15, 0.49
diffexpprob, diffexpdownprob, diffexploc, diffexpscale = 0.025, 0., 1.0, 1.0
bcv_dispersion, bcv_dof = 0.448, 22.087
nproggenes, progdownprob, progdeloc, progdescale = 1000, 0., 1.0, 1.0
progcellfrac, proggroups = .35, [1, 2, 3, 4]
minprogusage, maxprogusage = .1, .7
np.random.seed(seed)

# cells: group and library size
groupid = np.random.choice(np.arange(1, ngroups + 1), size=ncells, p=[1 / ngroups] * ngroups)
libsize = np.random.lognormal(mean=libloc, sigma=libscale, size=ncells)
cellnames = ['Cell%d' % i for i in range(1, ncells + 1)]
cellparams = pd.DataFrame({'group': groupid, 'libsize': libsize}, index=cellnames)

# genes: base means with outliers
basegenemean = np.random.gamma(shape=mean_shape, scale=1. / mean_rate, size=ngenes)
is_outlier = np.random.choice([True, False], size=ngenes, p=[expoutprob, 1 - expoutprob])
gene_mean = basegenemean.copy()
gene_mean[is_outlier] = np.random.lognormal(mean=expoutloc, sigma=expoutscale, size=is_outlier.sum()) * np.median(basegenemean)
genenames = ['Gene%d' % i for i in range(1, ngenes + 1)]
geneparams = pd.DataFrame({'gene_mean': gene_mean}, index=genenames)

# activity program: last nproggenes genes up-regulated
prog_gene = np.zeros(ngenes, bool); prog_gene[-nproggenes:] = True
DEratio = np.random.lognormal(mean=progdeloc, sigma=progdescale, size=nproggenes)
DEratio[DEratio < 1] = 1 / DEratio[DEratio < 1]
ratio = np.ones(ngenes); ratio[-nproggenes:] = DEratio
prog_mean = gene_mean * ratio
cellparams['has_program'] = False
cellparams['program_usage'] = 0.0
for g in proggroups:
    cells = cellparams.index[cellparams['group'] == g]
    has = np.random.choice([True, False], size=len(cells), p=[progcellfrac, 1 - progcellfrac])
    cellparams.loc[cells[has], 'has_program'] = True
    cellparams.loc[cells[has], 'program_usage'] = np.random.uniform(minprogusage, maxprogusage, size=has.sum())

# identity groups: DE genes per group (program genes excluded)
group_mean, truth = {}, []
for g in range(1, ngroups + 1):
    isDE = np.random.choice([True, False], size=ngenes, p=[diffexpprob, 1 - diffexpprob])
    isDE[prog_gene] = False
    r = np.random.lognormal(mean=diffexploc, sigma=diffexpscale, size=isDE.sum())
    r[r < 1] = 1 / r[r < 1]
    gr = np.ones(ngenes); gr[isDE] = r
    group_mean[g] = gene_mean * gr
    truth += [(x, f"group{g}") for x in np.array(genenames)[gr > 1]]
truth += [(x, "activity") for x in np.array(genenames)[prog_gene]]

# cell x gene means: group profile, mixed with the program profile by usage, scaled to library size
G = np.vstack([group_mean[g] / group_mean[g].sum() for g in range(1, ngroups + 1)])
P = prog_mean / prog_mean.sum()
u = cellparams['program_usage'].values[:, None]
mu = (1 - u) * G[groupid - 1] + u * P[None, :]
mu *= (libsize / mu.sum(1))[:, None]

# biological coefficient of variation, then Poisson counts
bcv = bcv_dispersion + 1 / np.sqrt(mu)
bcv = bcv * np.sqrt(bcv_dof / np.random.chisquare(bcv_dof, size=ngenes))
mu = np.random.gamma(shape=1 / bcv ** 2, scale=mu * bcv ** 2)
counts = sparse.csr_matrix(np.random.poisson(lam=mu).astype(np.float32))
del mu, bcv

adata = ad.AnnData(counts, obs=pd.DataFrame(index=cellnames), var=pd.DataFrame(index=genenames))
adata.write_h5ad(os.path.join(out_dir, "demo_counts.h5ad"), compression="gzip")
pd.DataFrame(truth, columns=["gene", "program"]).to_csv(os.path.join(out_dir, "ground_truth_programs.tsv"),
                                                        sep="\t", index=False)
cellparams.to_csv(os.path.join(out_dir, "cellparams.csv"))
print(f"Simulated {ncells:,} cells x {ngenes:,} genes | {ngroups} groups + 1 activity program | "
      f"activity cells: {int(cellparams['has_program'].sum()):,} | counts sparse: {sparse.issparse(adata.X)}")
