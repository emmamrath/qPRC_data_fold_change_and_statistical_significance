qPCR data fold change and statistical significance
==================================================

Calculate the fold change and whether the fold change is statistically significant, for qPCR data.  

This code has been validated against output of Windows program REST for the same input data.  

This code was used in the following paper which can be considered the citation for this code:  
BAMLET kills chemotherapy-resistant mesothelioma cells, holding oleic acid in an activated cytotoxic state.
Rath EM, Cheng YY, Pinese M, Sarun KH, Hudson AL, Weir C, Wang YD, Håkansson AP, Howell VM, Liu GJ, Reid G, Knott RB, Duff AP, Church WB.
PLoS One. 2018 Aug 29;13(8):e0203003. doi: 10.1371/journal.pone.0203003. eCollection 2018.
PMID: [30157247](https://www.ncbi.nlm.nih.gov/pubmed/?term=30157247) PMCID: [PMC6114908](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6114908/) DOI: [10.1371/journal.pone.0203003](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0203003)

Citation for [REST](https://www.gene-quantification.de/rest-2009.html) Windows software:  
Relative expression software tool (REST) for group-wise comparison and statistical analysis of relative expression results in real-time PCR.
Pfaffl MW, Horgan GW, Dempfle L.
Nucleic Acids Res. 2002 May 1;30(9):e36.
PMID: [11972351](https://www.ncbi.nlm.nih.gov/pubmed/?term=11972351) PMCID: PMC113859 DOI: [10.1093/nar/30.9.e36](https://academic.oup.com/nar/article/30/9/e36/1089004)

### Command line to call this program

```
Rscript qPRC_data_fold_change_and_statistical_significance.R <input_file> > <output_file>
Rscript qPRC_data_fold_change_and_statistical_significance.R example_data.txt > example_output.txt
```

Please note that the following are currently hard-coded and need to be converted into input parameters so as to generalise this script:  
* treatment vs control cell-lines
* target genes vs control gene

### Format of input qPCR data file

```
cell_line	gene	threshold_cycle	Ct_mean	Ct_SD
MeT5A	ATP5G1	27.673	27.139	0.559
MeT5A	ATP5G1	27.186	27.139	0.559
MeT5A	ATP5G1	26.558	27.139	0.559
REN	ATP5G1	27.516	27.202	0.505
REN	ATP5G1	26.619	27.202	0.505
REN	ATP5G1	27.471	27.202	0.505
MeT5A	RNA18S1	7.651	7.805	0.223
MeT5A	RNA18S1	8.060	7.805	0.223
MeT5A	RNA18S1	7.704	7.805	0.223
REN	RNA18S1	5.199	5.167	0.200
REN	RNA18S1	4.953	5.167	0.200
REN	RNA18S1	5.350	5.167	0.200
```

