**Web**
* DAVID 2021 (Dec. 2021), DAVID Knowlegebase (v2023q4, updated quarterly): Has background problem according to "Pop Total" values changing between KEGG and DAVID results. Fold enrichment is reported. FDR values are reported. Although pathways without overlaps are not reported, an analysis of the results is consistent with no FDR problem.
* PANTHER v19: reported background size is consistent with no background problem. The full results included sets with no overlap. Fold enrichment is reported. FDR values are reported. No FDR problem.
* ENRICHR ( June 8, 2023 update): p-values and reported odds ratios are consistent with no background problem. I submitted a shortlist of 200 genes. The FDR values were consistent with discarding pathways with no overlap.
*  KOBAS-i (KOBAS 3.0): No enrichment score provided. But p-values are much mmore significant than when run on DAVID, so it looks like un-annotated genes are not removed from the background. After loading p and fdr values into R it is clear that pathways without any overlaps are not reported and are not included in FDR correction.
* WebGestalt 2024: results similar to DAVID. Enrichment score reported (ratio). Comparison with oratool shows genes with no annotation are removed from the background. Only gene sets with overlaps reported, but analysis of adjusted p-values suggests no FDR problem.
* g:Profiler: No enrichment score. Background problem as "effective domain size" changes between databases like KEGG and Reactome. FDR values provided but no raw p-values. Analysis of p-values suggests no FDR problem.
* eggNOG: Not an enrichment tool. Exclude.
* STRING: Could not ascertain from the results of the tool, so I have emailed the admins.
* Shinygo: Enrichment scores are consistent with not discarding genes without any categories. According to the code, the number of genes in pathways in background genes is found before discarding. Therefore no FDR problem.

**Propriety web**
* IPA: Ingenuity Pathway Analysis

**Cytoscape plugins**
* BINGO: Stated universe size is different. Therefore unannotated genes are discarded. Analysis of the pvalues with a fg shortlist indicates that sets with no overlaps are discarded and not included in the FDR correction.
* ClueGO: When running analysis with different gene set libraries, the universe size was changing, indicating that unannotated genes were ddiscarded. The FDR values look to have a bug. Here are the top p-values I got: 0.01 0.01 0.01 0.01 0.02 0.02 0.02 0.02 0.02 0.02, and here are the top FDR values I got 0.84 1.00 0.84 0.91 0.68 0.82 0.61 0.68 0.62 0.82. The FDR values I derived were all 0.8227273. Results are consistent with excluding sets with no overlaps. (724 Reactome terms from 7787 genes only)

**R/Bioconductor**
* clusterprofiler/DOSE: According to documentation, results and forum posts background genes without annotations are discarded. Results are consistent with FDR problem. No enrichment score reported. Nominal p-values and FDR values are reported.
* topGO: No enrichment score, clumsy workflow, Background problem (vignette), does not report FDR adjusted significance values, although all gene sets are reported even ones with no overlap.
* GOseq: No enrichment score, Background problem but has option to force behaviour `use_genes_without_cat=TRUE`.  All gene sets are reported, even those with no overlap. 
* goana/kegga: No enrichment score reported.  All gene sets are reported, even those with no overlap. No FDR values provided. According to the documentation and source code, the background problem is only an issue when using default of whole genome background; when using custom background, unannotated genes are not discarded.
* FGSEA::FORA no enrichment score. p-values and FDR values. No exclusion of unannotated genes. No exclusion of pathways with no overlap.
