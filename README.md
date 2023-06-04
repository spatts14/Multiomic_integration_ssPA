# Multiomic_integration_ssPA_network_approach
Imperial College London, MRes Data stream (2022-2023)
Semester one project: Multi-omic integration using pathway level analysis

**OBJECTIVES AND AIMS:**
*The aim of this project is to develop a network-based approach to integrate multi-omics
datasets using pathway scores, with a particular focus on incorporating metabolomic platforms.*
1. Generate correlation networks using Spearman correlation between pathway scores for both metabolomics and proteomics datasets.
2. Generate overlap network with calculated overlap coefficient, which represents the overlap between molecules in pathways for metabolomics and proteomics to remove pathway correlations that result solely from high overlap.
3. Filter correlation networks to remove edges from highly overlapping pathways.
4. Integrate metabolomic correlation and proteomic correlation networks into a single
correlation network.
5. Simulate pathway pair enrichment using semi-synthetic data to validate workflow’s ability
to detect statistically significant correlations.

## General Workflow 
 A. Analyze metabolomic and proteomic dataset individually
 
 B. Calculate overlap coefficent 
 
 C. Integrate datasets and vizualize datasets using a network-based approach


<img width="911" alt="Screen Shot 2023-04-12 at 10 09 16 AM" src="https://user-images.githubusercontent.com/116185211/231410385-82553f3b-a295-4457-8342-f96a39f1c765.png">

## Packages and software used:
*Python packages:*
ssPA package ([GitHub py-ssPA](https://github.com/cwieder/py-ssPA))
NetworkX ([GitHub NetworkX](https://github.com/networkx/networkx))

*Software:*
Cytoscape 


## Code workflow
1. Import data
2. Pre-process data (remove outliers, remove samples with high % missing data, impute missiing data, log2 transform data, standardize data)
3. Confirm protein/metabolite/transcript name is in the correct format (i.e protein is represented by UniPort ID, metabolites are ChEBI ID)
4. Perform ssPA analysis - choose correct ssPA function (ssGESA, kPCA, ssClusterPA, etc.) 
```ruby
kpca_scores = sspa.sspa_kpca(uniport_proteomic, reactome_uniprot)
```
To confirm the correct ssPA function for your analysis, please see Wieder C in the reference.
For further information on ssPA package, see ([ssPA package](https://github.com/cwieder/py-ssPA))

5. Calculate correlation between pathways 
```ruby
correlation_t_test, pval = stats.spearmanr(kpca_scores, axis=0, nan_policy='propagate', alternative='two-sided')

#Put into a dataframe and pval
correlation_t_test = pd.DataFrame(correlation_t_test, columns = kpca_scores_corr.columns, index=kpca_scores_corr.index)
pval_corr = pd.DataFrame(pval, columns = kpca_scores_corr.columns, index=kpca_scores_corr.index)
```
6. Perform multiple testing correction
```ruby
pval_corr_list = pval_corr.stack().reset_index()
pvals_corrected_stats = statsmodels.stats.multitest.multipletests(pval_corr_list["p_val"], alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
reject = pvals_corrected_stats[0]
pvals_corrected = pvals_corrected_stats[1]
alphacBonf = pvals_corrected_stats[2]

# Put p value into a dataframe
pvals_corrected_df = pd.DataFrame(pvals_corrected, columns = ["p_val_adjust"])
```
8. Calculate pathway overlap using with the overlap coefficient (OC) (see notebooks for source code)
9. Remove edges with correlations with p > 0.005 and OC > 0.5
10. Generate networks using NetoworkX
11. Format Networks using Cytoscape

For network integration
1. Keep pathways found in both molecular networks (i.e remove pathways only present in one molecular network)
   * Of note, use networks without correlations or OC values removed (i.e skip step 9)
3. Combine pathway scores from each molecular network (e.g average pathway scores if you would like to weight each network evenly)
4. Remove edges with correlations with p > 0.05 and OC > 0.5
5. Generate networks using NetoworkX
6. Format Networks using Cytoscape

## Abstract
Rapid advances in high-throughput technologies have resulted in the generation of vast amounts of omics data, demanding innovations in analysis to provide the richest interpretation and understanding of the biological system. Pathway analysis (PA) has emerged as an exciting tool to gain novel insights and give meaning to underlying biology from these large omics datasets. PA overcomes the challenge of interpreting large datasets by combining biological knowledge from datasets with various mathematical analysis, computational algorithms, and statistical testing to reduce complexity. Single sample pathway analysis (ssPA), a non-conventional PA technique, derives pathway scores from high-throughput omics data, such as proteomics and metabolomics, and estimates the level of pathway intensity for each sample in a dataset. Correlations of pathway scores between pathways are computed and visualized in a network to improve interpretability. In such networks, pathways are represented as nodes and correlations scores are represented as edges. While there are many strategies already available to analyze high-throughput omics data, the methodology described in this thesis is the first of its kind to integrate multiple omics at the pathway level using a network-based approach. Using this proposed method, we have shown that we can integrate proteomic and metabolomic datasets at the pathway level and identify clusters of highly correlated pathways from the resulting networks. We hypothesize this approach will aid researchers in interpreting and visualizing omics data to identify unique and intriguing patterns. Additionally, this approach could be used as a visualization tool to better understand which pathways are altered in a given disease. 

## Methods:
The multi-omics data used in this work was sourced from publicly available metabolomic and proteomic datasets from Su et al., 2020 23. The data can be accessed as detailed in Su et al. in “Data and Code Availability”. The patient cohort contains samples collected from 258 healthy controls and 139 COVID19 patients of varying severity shortly after diagnosis. Disease severity was represented using WHO status as follows: 0 (healthy), 1-2 (mild), 3-4 (moderate), 5-7 (severe). Metabolomic and proteomic data in the form of n × m  matrices, where rows represent n samples and columns represent m features, was processed prior to pathway analysis. Processing included: (1) removal of outliers (2) removal of features with more than 60% (metabolomic) or 40% (proteomic) of data missing and donor samples with more than 50% of features (3) impute missing data (4) 〖log〗_2-transformed to normalize their distribution (5) standardized the feature columns to a ( μ=0 and σ=1). 

## Results:
Correlation between pathway scores were computed using the Spearman correlation algorithm. We also calculated the p-value for each correlation and corrected them using the Bonferroni multiple correction method. Szymkiewicz-Simpson Overlap Coefficient (OC) was used to calculate the overlap between pathway metabolites and proteins. An OC of 0 indicates there is no overlap between pathway “A” and pathway “B”, whereas an overlap of 1 indicates that the smaller set is a subset of the larger set or the pathways contain the same metabolites or proteins. The NetworkX Python package was used to create overlap and correlation networks. Cytoscape was used to visualize the network using an edge-weighted spring embedded layout. 
The overlap network was created using the calculated OC for pathway pairs. A pairwise comparison was computed for each pathway resulting in an OC value, which was used as an edge in the network. Nodes represent pathways and weighted undirected edges represent the OC between a pair of pathways. Edge color denotes OC, where light colors have a low OC and dark colors have a high OC. The edge threshold for the network in both the proteomic and metabolomic overlap network was set on OC ≥ 0.5, therefore any edges with OC ≥ 0.5 were shown. This was done to reduce the density of the network and improve visualization. 
The final correlation network is constructed where nodes represent pathways and weighted undirected edges represent the absolute value of the Spearman correlation between pathway pairs. For all networks, edges with an OC ≥ 0.5 were removed from the network in order to remove edges with high overlap. Additionally, Spearman correlations with Bonferroni corrected p ≥ 0.005 (proteomic and metabolomic) or p ≥ 0.05 (integrated) were removed from the network. This thresholds was set to remove non-statistically significant correlations and to improve visualization of networks. 
In the integrated correlation network, metabolomic and proteomic networks were combined in such a way that pathways present in both datasets were kept, while pathways that appeared in only one dataset were discarded. Correlation scores, OC, and PS for each WHO status for the integrated network were calculated using respective values from the proteomic and metabolomic networks.
The ultimate goal of this methodology is to aid scientists in interpreting and visualizing omics data to identify unique and intriguing patterns in their datasets that may elucidate biological mechanisms. Thus, we validated that this approach could identify relevant pathways in a given disease. We confirmed that pathways identified in this case study were relevant by consulting available COVID19 literature. The innate immune system was identified by our method as a relevant pathway in the dataset. While the innate immune system is quite broad, it has been widely implicated in COVID19 infection. This is rather unsurprising that the innate immune system is involved, as SARS-CoV-2 is a virus, and the innate immune system is the first line of defense against a viral infection.  Nevertheless, we have shown that we can integrate proteomic and metabolomic datasets at the pathway level and identify clusters of highly correlated pathways from the resulting networks.


## Conclusions: 
The methodology described in this work, to the best of our knowledge, is the first of its kind to integrate multiple omics at the pathway level using a network-based approach. In particular, our approach focuses on metabolomics as one of the omics platforms, where many other methods are developed using genomic models. Additionally, we propose that this approach is more user-friendly to the average biologist, making pathway level analysis more accessible to the community. Lastly, integrating multi-omics using this method can aid in elucidating biologically relevant pathways and may provide insights into the underlying mechanisms of disease through highlighting which pathways are highly correlated. The individual aspects of this method are not novel, as multi-omic data integration, PA, and network-based approaches are widely used. However, we believed that this approach is the first time that multi-omic data integration, PA, and network-based approaches have been combined in this specific way with a particular focus on incorporating metabolomic datasets.

## References:
Wieder C, Lai RPJ, Ebbels TMD. Single sample pathway analysis in metabolomics: performance evaluation and application. BMC Bioinformatics. 2022;23(1).



