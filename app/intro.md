
We have noticed two subtle problems with some ORA tools.
The first problem we call the "Background problem", and involves genes without
annotated categories being discarded erroneously.
The second problem we call the "FDR problem" because some tools incorrectly
apply false discovery rate correction of p-values.
For more information, refer to our [README](https://github.com/markziemann/background/tree/main).
This tool is released under an MIT Licence and comes without any warranty.

### How to use the tool

To get started, you will need to upload two `.txt` files containing gene symbols;
one for the foreground genes and one for the background.
The foreground list contains genes that were identified by some 'omics
analysis as being interesting.
The background list contains all the genes that were detected robustly in the
omics assay.
Here "robustly" means that the detection was good enough that it *could* have
been part of the foreground list.

After uploading the lists, select the gene set library you want to test.
We have sourced some popular options from [MSigDB v2023.2](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp).
Then select the problem/error you would like to characterise; the options are
"Background error", "FDR error" and "Both errors".

Then click on the "Data" tab which will summarise the gene lists that have been
uploaded.
It is a good idea to check these values to ensure that the list was uploaded in
its entirety.
Then click on the "Comparative Analysis" tab which will kick off the enrichment
analysis.
This part could take up to 30 seconds to complete, so be patient and avoid
clicking away or reloading the page.
The table will eventually appear.
By default it is sorted with sets that have divergent FDR values at the top.
Then click on the "Charts" tab to view some plots.

The "Download Report" button can then be used to generate a full report for 
provenance.
Again it can take 30-60 seconds to generate so be patient and don't click away
until the file has been downloaded.