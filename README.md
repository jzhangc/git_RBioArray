# git_RBioArray
Simple to use package for both microarray and RNAseq data analysis (limma based)

Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        source("https://bioconductor.org/biocLite.R")
      
        biocLite()
    
  - Install stable release
  
        devtools::install_github("jzhangc/git_RBioArray/RBioArray", repos = BiocInstaller::biocinstallRepos())   

  - Install development build
  
        devtools::install_github("jzhangc/git_RBioArray/RBioArray", repos = BiocInstaller::biocinstallRepos(), ref = "beta")  

Update log

    0.4.7 (feature preview)
    (ICEBOX)
      - New clustering funcions:
        - Bayesian biclustering function rbioarray_bbc()
        - K-mean clustering function rbioarray_kmean()
        
      - Updates to RNAseq functions:
        - rbioseq_ImportCount() function added to import HTseq counted files to R environment
        - rbioseq_transform() function added for transforming data for clustering analysis and any feature selection/classcification processes
        - rbioseq_hcluster() function re-written with rbioseq_transform() function incorporated
        - compositional analysis methods added to the limma DE method due to the availability of clr transfromation in rbioseq_transform()
        - DESeq2-based method added to rbioseq_DE() function
          - the current understanding is that DESeq2 method contains a compositional analysis mode, see rlog() function from the DESeq2 package
          - RBioArray implementation of DESeq2 method uses compositional analysis by default
        
      - Updates to microarray functions:
        - MA plot option added for functions rbioarray_flt() and rbioarray_PreProc()
        
      - Updates to correlation functions:
        - rbioarray_corcluster_super() now supports VLR (Log-Ratio Variance) correlation and differential proportionality for NGS compositional data analysis 
      
    (ADDED)
      - updates to microarray functions:
        - rbioarray_DE() now has the options to produce csv files for either (i) all probes, (ii) significant probes, (iii) all probes with gene name, (iv) significant probes with gene name. However, DE reuslts for all probes will be exported to the R environment regardless of these settings. Similarly, F stats is also always exported to the working directory and the R environment regardless of these settings. 
        - rbioarray_DE() with the "FORK" cluster module re-written for foreach style parallel computing
        - To keep things consistent with limma's Elist, the "target" component list output changed to "targets"
      
      - Bug fixes
      
      
    0.4.6 (4.26.2018)
      - Centered log transformation and isometric log transformation function rbioseq_clr_ilr_transfo() added
      - Normalizatin method arugment norm.method added to rbioseq_DE() function
        - options are: "TMM","RLE","upperquartile","none". NOTE: clr and ilr options to be added in the next release
      - Correlation function rbioarray_corcluster_super() now outputs a summary csv file for the significance test
      - FDR corretion option added for correlation p values
      - The "q.value" argument for all the applicable functions changed to "sig.p"
      - DE method argument adjusted for rbioarray_hcluster_super()
      - pcutoff argument for rbioarray_hcluster_super set to "DE.p.sig"
      - DE argument for rbioarray_DE() adnd rbioarray_venn_DE changed to "sig.method"
      - Clustering functions separated into cluster.R file
      - Stats variable argument statsVar added to rbioGS_kegg(), so that P value and FC can be used. 
      - Additional argument check added for rbioGS_kegg()
      - Rightside y-axis now uses a function from RBioplot pakcage, which now is a dependency
      - Bug fixes
      

    0.4.5
      - Correlation p value function cor_pvalue() added
      - Updates made to rbioarray_corcluster_super()
        - Pearson correlation significance plot added to rbioarray_corcluster_super() via argument sigPlot = TRUE
        - Spearman correlation method added, controlled via argument method
        - Various arguments added for the significance plot
        - Alpha can be set for the p values via argument cor.sig
        - The function now outputs a p-value matrix
        - Output file name to the significant plot is ".sigplot.pdf"
        - Output file name to the correlation hclust heatmap change to ".corheatmap.pdf"
        - The "none" option added for DE argument for rbioarray_venn_DE()
        - Values for DE options are no longer case sensitive
        - A bug fixed where the control variable check screens DE dataframe
      - Additional argument checking mechanisms added for functions:
        - rbioarray_hcluster_super()
      - Other bug fixes
      
    
    0.4.4
      - Gene repeats processing methods added via argument combineGeneDup for rbioarray_flt()
      - rbioarray_corcluster_super() now also outputs a correlation matrix file to the working directory
      - Additional arguments added to rbioarray_hclust_super() to enhance data compatibility
      - Additional arguments added the non-Agilent datasets and datasets without a control type variable for functions:
          - rbioarray_flt()
          - rbioarray_DE()
          - rbioarray_hclust()
      - Additional argument checking mechanisms added for functions:
          - rbioarray_flt()
          - rbioarray_DE()
          - rbioarray_hcluster_super()
          - rbioarray_corcluster_super()
      - Unified argument fltlist set for all the functions that take the filtered data list
      - Filtering without control probes method added to rbioarray_flt()
      - The argument for annotation files changed from "anno" to "annot"
      - annot = NULL is now acceptable to rbioarray_DE()
      - DE method "none" added for rbioarray_DE()
      - objTitle argument value now added as prefix to the output file names for rbioarray_DE()
      - A bug fixed for rbioarray_DE() where 1 was set for pcutoff when no significant target found under FDR
      - Codes cleaned up
      - Other bug fixes
      

    0.4.3
      - A bug fixed for rbioarray_venn_DE() when p value threshold method set to "fdr"
      - Other bug fixes
      

    0.4.2
      - Annotation matrix now can be accpeted by rbioarray_corcluster_super()
      - Removing Agilent microarray control funtionality added to rbioarray_corcluster_super()
      - rbioarray_hclust_super() updated for a better data format compatability
      - Bug fixes
      

    0.4.1
      - rbioarray_hclust() updated for a better data format compatability
      - Bug fixes
      

    0.4.0
      - rbioarray_corcluster_super() added: Supervised Pearson correlation clustering analysis and heatmap
      - rbioarray_hcluster_super() updated with more concise codes
      - rbioseq_de() now outputs filtered and normalized read counts matrix to the environment
      - Bug fixes
      

    0.3.4
      - GS functions updated
      - Bug fixes
      

    0.3.2 - 0.3.3
      - Medaka added for rbioGS_sp2hsaEntrez function
      - EntrezGeneID option added for rbioGS_sp2hsaEntrez function
      - Bug fixes
      
    
    0.3.0
      - RNAseq data DE analysis function rbioseq_DE() added
      - Unsupervised hcluster and heatmap functioin rbioseq_hcluster() added
      - RNAseq data preprocessing function added to rbioseq_DE()
      - Message display added to the DE functions
      - heatmap and Venn functions updated with correct variable arguments
      - Bug fixes
      

    0.2.7 - 0.2.10
      - Boxplot and Scatter plot functions for GS optimized
      - Bug fixes
      
    
    0.2.6
      - Both PSOCK and FORK cluster types are now avaiable for all the functions featuring parallel computing
      - Small "quality of life" changes
      - Bug fixes
      

    0.2.5
      - Direct loading GS database files functionality added for GS functions
      - Bug fixes
      
    
    0.2.2 - 0.2.4
      - mmu/rno to hsa entrez ID converesion function overhaul
      - Arguments for input DE data list for GS functions unified as "DElst"
      - Arguments for parallel computing unified for all the relevant functions as "parallelComputing" and "clusterType"
      - Bug fixes
      

    0.2.1
      - Venn diagram function now outputs a csv containing overlapping genes/probes
      - Multicore support added for venn diagram function
      - Code optimization for venn diagram function
      - Removing probes without gene symbol functionality added for unsupervised hcluser function
      - DE function now outputs F stats into both the environment and work directory
      - Bug fixes
      
    
    0.2.0 
      - Venn diagram function added
      - Array weight added to data preprocessing function for Elist objects
      - Missing depndencies for parallel computing added for various functions
      - Bug fixes
      
    
    0.1.34 
      - Elist objects now can be properly recognized by all the functions
      

    0.1.33
      - Help page for unsupervised heat map updated
      

    0.1.32
      - Users can now choose to keep or remove the control probes (Agilent microarray platform) for rbioarray_hcluster() function
      
    
    0.1.31
      - FC filter added for supervised hclust heatmap function
      

    0.1.30
      - Volcano dots colours are now user customizable
      - Volcano dots annotation functionality added
      - KEGG visualization function updated
      - Bug fixes and other improvements
      
      
    0.1.29
      - Code redundancy further reduced for array DE function
      

    0.1.28
      - Array DE function re-written with better parallel computing (doParallel and foreach), 
      greatly reduced code redundancy, and backend prepration for plot annotation functionality
      - Bug fixes
      
      
    0.1.27
      - Supervised hclust heatmap function updated for label display functionality
      - Bug fixes
      
      
    0.1.26
      - Supervised hclust heatmap function added
      - Bug fixes
      
      
    0.1.25
      - Stand alone GSA and GSA plotting functions added
      - Bug fixes
      
      
    0.1.23 - 0.1.24
      - All-in-one GSA function re-written
      - Bug fixes
      
      
    0.1.22
      - GS plotting function updated
      - Bug fixes
      
      
    0.1.21
      - rbioGS_boxplot, rbioGS_scatter and rbioGS_all functions updated
      - bug fixes
      
      
    0.1.20
      - GSA methods for GS functions are now user customizable
      - Bug fixes
      
      
    0.1.19
      - Preliminary all-in-one GS function added
      - Bug fixes
      
      
    0.1.17 - 0.1.18
      - GS plotting functions added
      - Bug fixes
      
      
    0.1.15 - 0.1.16
      - GS function re-written
      - Bug fixes
      
      
    0.1.14
      - Part of the array DE function re-written
      - Bug fixees
      
      
    0.1.13
      - Plot section of the array DE function improved
      - Bug fixes
      
      
    0.1.10 - 0.1.12
      - Log transformation for array functionality added
      - hclust heatmap function added
      - Bug fixes
      
      
    0.1.9
      - Volcano plot functionality added for array DE function
      - Bug fixes
      
      
    0.1.8
      - Fit and contrasts re-written for array DE function
      
      
    0.1.7
      - Preparation work for the next update
      - Bug fixes
      
      
    0.1.6
      - Array DE function updated
      - Bug fixes
      
      
    0.1.5
      - Array DE function added
      
    
    0.1.2 - 0.1.4
      - Preliminary microarray functions added
      - Bug fixes
      
      
    0.1.0 - 0.1.1
      - Initial seq functions added
      - Initial GS functions added
      - Initial commit
