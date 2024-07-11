**Web**
* DAVID 2021 (Dec. 2021), DAVID Knowlegebase (v2023q4, updated quarterly): Has background problem according to "Pop Total" values changing between KEGG and DAVID results. Fold enrichment is reported. FDR values are reported. Although pathways without overlaps are not reported, an analysis of the results is consistent with no FDR problem.
* PANTHER v19: reported background size is consistent with no background problem. The full results included sets with no overlap. Fold enrichment is reported. FDR values are reported. No FDR problem.
* ENRICHR (8 June 2023 update): p-values and reported odds ratios are consistent with no background problem. I submitted a shortlist of 200 genes. The FDR values were consistent with discarding pathways with no overlap.
* KOBAS-i (KOBAS 3.0): No enrichment score provided. But p-values are much mmore significant than when run on DAVID, so it looks like un-annotated genes are not removed from the background. After loading p and fdr values into R it is clear that pathways without any overlaps are not reported and are not included in FDR correction.
* WebGestalt 2024: results similar to DAVID. Enrichment score reported (ratio). Comparison with oratool shows genes with no annotation are removed from the background. Only gene sets with overlaps reported, but analysis of adjusted p-values suggests no FDR problem.
* g:Profiler (13 Feb 2024 update): No enrichment score. Background problem as "effective domain size" changes between databases like KEGG and Reactome. FDR values provided but no raw p-values. Analysis of p-values suggests no FDR problem.
* STRING: Could not ascertain from the results of the tool, so I have emailed the admins. They report "Proteins without annotations are included in the default backgorund." and "We do not FDR correct or test for terms with 1 or less protein annotated, but that is for backgorund not for the (foreground) tested set."
* Shinygo v0.80: Enrichment scores are consistent with not discarding genes without any categories. According to the code, the number of genes in pathways in background genes is found before discarding. Therefore no FDR problem.

**Propriety web**
* IPA: Ingenuity Pathway Analysis

**Cytoscape plugins**
* BINGO (v3.0.5): Stated universe size is different. Therefore unannotated genes are discarded. Analysis of the pvalues with a fg shortlist indicates that sets with no overlaps are discarded and not included in the FDR correction.
* ClueGO (v2.5.10): When running analysis with different gene set libraries, the universe size was changing, indicating that unannotated genes were ddiscarded. The FDR values look to have a bug. Here are the top p-values I got: 0.01 0.01 0.01 0.01 0.02 0.02 0.02 0.02 0.02 0.02, and here are the top FDR values I got 0.84 1.00 0.84 0.91 0.68 0.82 0.61 0.68 0.62 0.82. The FDR values I derived were all 0.8227273. Results are consistent with excluding sets with no overlaps. (724 Reactome terms from 7787 genes only). No Enrichment score.

**R/Bioconductor**
* clusterprofiler (v4.12.0): According to documentation, results and forum posts background genes without annotations are discarded. Results are consistent with FDR problem. No enrichment score reported. Nominal p-values and FDR values are reported.
* topGO (v2.56.0): No enrichment score, clumsy workflow, Background problem (vignette), does not report FDR adjusted significance values, although all gene sets are reported even ones with no overlap.
* GOseq (v1.56.0): No enrichment score, Background problem but has option to force behaviour `use_genes_without_cat=TRUE`.  All gene sets are reported, even those with no overlap. 
* goana/kegga (limma v3.60.3): No enrichment score reported.  All gene sets are reported, even those with no overlap. No FDR values provided. According to the documentation and source code, the background problem is only an issue when using default of whole genome background; when using custom background, unannotated genes are not discarded.
* FORA (fgsea v1.30.0): no enrichment score. p-values and FDR values. No exclusion of unannotated genes. No exclusion of pathways with no overlap.

| Tool | Version | Type | Provides FDR values | Provides Enrichment Score | Proper background handling | Proper FDR | Reference |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DAVID | v2023q4 | Web | Yes | Yes | No | Yes | [@Sherman2022-jn] |
| Panther | v19 | Web | Yes | Yes | Yes | Yes | [@Mi2019-ak] |
| Enrichr | June 2023 | Web | Yes | Yes | Yes | No | [@Kuleshov2016-og] |
| KOBAS-i | KOBAS 3.0 | Web | Yes | No | Yes | No | [@Bu2021-ak] |
| WebGestalt | 2024 | Web | Yes | Yes | No | Yes | [@Elizarraras2024-nr] |
| g:Profiler | Feb 2024 | Web | Yes | No | No | Yes | [@Kolberg2023-nr] |
| STRING-DB | v12.0 | Web | Yes | Yes | Yes | Yes | [@Szklarczyk2023-an] |
| ShinyGO | v0.80 | Web | Yes | Yes | Yes | Yes | [@Ge2020-lu] |
| Ingenuity Pathway Analysis | ? | ? | ? | ? | ? | ? | [@Kramer2014-kd] |
| BinGO | v3.0.5 | Cytoscape | Yes | No | No | No | [@Maere2005-gq] |
| ClueGO | v2.5.10 | Cytoscape | Yes | No | No | No | [@Bindea2009-jn] |
| clusterProfiler | v4.12.0 | R package | Yes | No | No | No | [@Wu2021-wy] |
| topGO | v2.56.0 | R package | No | No | No | Yes+ | [@Adrian_Alexa2017-nf] |
| GOseq | v1.56.0 | R package | No | No | No* | Yes+ | [@Young2010-iw] |
| goana/kegga | limma v3.60.3 | R package | No | No | No** | Yes+ | [@Ritchie2015-oz] |
| fora | fgsea v1.30.0 | R package | Yes | No | Yes | Yes | [@Korotkevich2016-gd] |




