# PhyloWGS Witness

## Description
Witness is a web-based viewer for PhyloWGS results. It uses the JSON files
written by `write_results.py`. Refer to steps 4 and 5 in the [PhyloWGS
README](../README.md) for instructions on running `write_resuls.py` and
starting Witness.


## Using Witness
Upon opening Witness and selecting a run on the left side, you will see two
tabs in the upper-right corner, named "Tree summaries" and "Tree viewer".

If you see no runs available to load, double-check that the
`witness/index.json` file is not empty. If it is, you likely erred in following
the naming conventions listed in step 4 of the [PhyloWGS
README](../README.md)'s "Running PhyloWGS" section.

### Tree summaries
The tree summaries section provides summary statistics characterizing all the
trees sampled by PhyloWGS for your data. Note that several of these
visualizations are currently incorrect in multisample cases. Visualizations
include the following:

  1. VAF histogram: in single-sample runs, this will often show clear clusters.
     For multisample runs, the VAF given is simply the mean across samples, and
     so is meaningless and can be misleading.
  2. Tree cluster summary table
  3. Tree structure scatter plot. Ellipses are contour plots of Gaussian
     distributions used to cluster trees. See the "Description of Reported
     Values" section below for how to interpret the axes.
  4. Cellular prevalence histogram. As with the VAF histogram, for multisample
     runs, the values shown are simply the means across samples, and so can be
     misleading.
  5. Histograms of number of SSMs per subpopulation. Subpopulations in each
     tree are ranked by cellular prevalence.
  6. Histograms of number of subpopulations per tree.

### Tree viewer
The tree viewer lets you examine individual trees sampled by PhyloWGS. This
includes a rendering of the tree showing structure, in which the area of each
tree node corresponds to the number of SSMs in the associated subpopulation.
Moreover, you can see a line chart showing how the cancer cell fraction changes
from sample to sample, and a table providing statistics about each
subpopulation.

On the right side, you can select particular tree samples to view. For each
tree, you will see the sample index, nlgLH, number of nodes, cluster
membership, LI (linearity index), BI (branching index), and CI (clustering
index). See the next section to understand these values.

For each tree, we also list the nlgLH (normalized log likelihood). Given how we
normalize the log likelihood, it can be understood as "average reconstruction
error in bits per SSM, per sample". Consequently, lower nLgLH values indicate
that a tree better fits the data. If you must pick only one tree, the tree with
the lowest nLgLH is a good candidate; it is better, however, to examine the
representative trees reported for each tree cluster to gain a better
understanding of the full range of trees your data supports.

Tree indices
------------
To summarize trees, we generate three indices for each, then cluster trees based on these values:

1. LI (linearity index): proportion of mutation pairs are in linear relations --
  i.e., given mutations `A` and `B`, the pair `(A, B)` is in a linear
  relationship if `A` is in a population ancestral to `B` or vice versa.
2. BI (branching index): proportion of mutation pairs that occur in different branches of the tree.
3. CI (clustering index): proportion of mutation pairs that are placed in the same cluster (i.e., population).

As every mutation pair must be in one of the above three relationships, we have
`LI + BI + CI = 1`.
