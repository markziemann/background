# Two subtle problems with over-representation analysis

Over-representation analysis is frequently used to determine the functionally
important themes from gene lists.
These gene lists are classified in different ways, but are commonly defined
by differential regulation in omics experiments like RNA-seq, proteomics or
micro-array.
Many different tools are available to conduct this type of analysis including
[DAVID](https://david.ncifcrf.gov/) and [Enricher](https://maayanlab.cloud/Enrichr/).

For the results to be reliable, the gene list needs to be provided together
with a background list.
The background list is the full list of genes that were reliably detected by
the assay.
This is important, because these assays typically cannot detect all genes
reliably due to technical limitations.
As different tissues express genes very differently, the right background list
is critical to obtain the functional categories that are specific to the
experiment conducted.
Using a generic background list of all annotated genes could lead to a major
bias.

Additionally, as such analyses involve the analysis of hundreds of
functional categories in parallel, the resulting p-values cannot be taken at
face value due to the high rate of false positives.
To mitigate this, false discovery rate correction is applied to adjust the
p-values.

In this work we highlight two subtle issues with some popular implementations 
of over-representation analysis.

## Genes without annotated categories are discarded

The first issue related to how the background list is defined.
If a user does the right thing and provides a background list along with the
foreground list, then some tools will exclude any genes from the background
that do not have any annotated functional categories.
This results in a large fraction of genes in the background being excluded,
and the enrichment ratio being unreliable.
This issue is more severe for analyses that involve smaller libraries of
functional annotations.
For example KEGG has annotations for a few thousand genes, meaning most of
the background genes get discarded.

## p-value adjustment is applied incorrectly

The second issue is that the false discovery rate is implemented incorrectly.
ORA enrichment tools intersect the gene list with each functional category
and return results for each functional category with at least one common
gene.
But there are some functional categories with zero genes in common with the
gene list.
These are excluded from the end results.
The process of attempting to intersect these categories with the submitted
gene list is itself a test, and so statistically speaking, the result of
that failed test should be reported in the final results with no overlaps and
a p-value of 1.
This may not be a problem when running enrichment analysis on a gene list 
with 3000 members, but it might be serious when it is small, like <200
members.
This results in adjusted p-values being slightly smaller than they should be.
