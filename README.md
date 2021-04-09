# git_RBioArray
Simple to use package for both microarray and RNAseq data analysis (limma based)

To cite in publication:
  
    Zhang J, Hadj-Moussa H, Storey KB. 2020. Marine periwinkle stress-responsive microRNAs: a potential factor to reflect anoxia and freezing survival adaptations. GENOMICS. 2020 Jul 27: S0888-7543(20)30169-5. doi: 10.1016/j.ygeno.2020.07.036.
    Zhang J, Wallace SJ, Shiu MY, Smith I, Rhind SG, Langlois VS. 2017. Human hair follicle transcriptome profiling: a minimally invasive tool to assess molecular adaptations upon low-volume, high-intensity interval training. Physiological Reports. 5(23) pii: e13534. doi: 10.14814/phy2.13534.

Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        if (!requireNamespace("BiocManager"))
            install.packages("BiocManager")
            
        BiocManager::install()
    
  - Install stable release
  
        devtools::install_github("jzhangc/git_RBioArray/RBioArray", repos = BiocManager::repositories())   

  - Install development build
  
        devtools::install_github("jzhangc/git_RBioArray/RBioArray", repos = BiocManager::repositories(), ref = "beta")  

Update log

    0.5.5 (feature roadmap)
    (ICEBOX)
      - General updates
        - splines pacakge added as a depedencing for continuous outcome support
        
      - New microarray functions:
        - New DE analysis function added: rbioarray_de_analysis()
          - The fnction has export options:
            (i) all features, 
            (ii) features with an annotation name, 
            (iii) significant features with or without annotation name depending on the setting for argument "gene_symbol"
            
      - New clustering funcions:
        - Bayesian biclustering function rbio_bb()
        - Propotionality correlation function added for CLR transformed RNAseq data
        - rbio_kmeans() supports methods for rbioarray_de and rbioseq_de objects for automatic K means cluster
      
      - New network analysis functions:
        - SNF network and network fusion functions
        - KNN graph construction
      
      - Updates to microarray functions
        - MA plot option added for functions rbioarray_transfo_normalize() and rbioarray_filter_combine()
        - Relevant functions updated with continuous outcome support
        
      - Updates to RNAseq functions
        - Relevant functions updated with continuous outcome support
        - rnaseq_de() output "rbioseq_de" class now includes "voom_output" and "calcNormFactors_outpout" (filtered raw count with lib size)
        
      - Revamp GS functions
      
      - Other updates
        - Functions updated for R Notebook/Markdown compatibility
        - Combine documentation for methods of different classes

    (ADDED)
      - New ML functions:
        - rbio_randomforest_fs() added for recursive randome forest feature selection analysis
          - RBioFS package is now part of the dependencies
        
      - Updates to clustering funcions:
        - rbio_supervised_hcluster() code base updated for better data compatibility
      
      - Updates to DE significant test functions
        - sig() updated so that the volcano plots would have the export.name as part of the file name
        

    0.5.4 (March.24.2021)
      - General updates
        - Typo fixed
        - Help page updated for data input functions regarding data column (sample) order/name and the sample/name order from the annotation data
      
      - New RNAseq functions
        - rbioseq_gtf(): a faster function to replace rbioseq_import_gtf()

      - New clustering funcions:
        - K means clustering function rbio_kmeans()
          - K means clustering plotting function with PCA functionality rbio_kmeans_plot() added

      - New network analysis functions:
        - rbio_tom() added for TOM (topological overlap measure) analysis
          - The output "rbio_tom_graph" object contains an igraph object for network visualization
          - The function also supports multiple automated tree cuting techniques, such as WGCNA's dynamic tree cutting
        - rbio_network() added for network constuction and visualization
      
      - Updates to GSA functions
        - A bug fixed for rbioGS_sp2hsaEntrez() due to updated BiomaRt package
        - Manual page updated for rbioGS_sp2hsaEntrez() to reflect the new S3 DE classes
      
      - Updates to microarray functions
        - rbioarray_transfo_normalize() now only supports quantile method for normalization
        - rbioarray_filter_combine() now can deal with "none or all" filtered situations

      - Updates to RNAseq functions
        - rbioseq_import_count() re-written:
          - rbioseq_import_count() function now requre sample id variable name
          - rbioseq_import_count() function now accepts data.frames
          - The "rbioseq_count" object now have the same sample order for raw_read_count and targets
        - Updates to rbioseq_import_gtf()
          - it is now a legacy function for comptability
          - rbioseq_import_gtf() not outputs a data frame
          - rbioseq_import_gtf() now outputs feature length
        - A bug fixed for rnaseq_de() where filtering not working properly
        - A bug fixed for sig.rbioseq_de() where the voom output was not used as normalized E matrix
      
      - Updates to DE significant test functions
        - The "sig" objects from sig() function now includes F_stats in the input_data item
        - Manual updated with added explanation on the "thresholding_summary" 
        - When gene_symbol = TRUE, sig() will only remove genes without a symbol from the volcano plot. 
          All other results items still contain the results for all genes. 

      - Updates to cluster funcions:
        - rbio_unsupervised_hcluster() now supports "rbioarray_de" objects
        - rbio_unsupervised_hcluster() now outputs a list object including distance and cluster results (for network analysis)
        - rbio_unsupervised_hcluster() now accepts user defined export name prefix via the export.name argument
        - rbio_supervised_hcluster() now outputs a list object including distance and cluster ressults (for network analysis)
        - rbio_unsupervised_corcluster() and rbio_supervised_corcluster() now export correlation results into the environment
          - With suffixes of "_cor_unsuper" and "_cor_sig", respectively
        - rbio_unsupervised_corcluster() and rbio_supervised_corcluster()'s heatmap (and sigplot) now can be turned off
        - rbio_supervised_hcluster() how automatically cluster f_stats significant results 
        - hcluster methods added for rbio_unsupervised_corcluster() and rbio_supervised_corcluster()
        - Manual page updated for rbio_unsupervised_corcluster(), along with its S3 class types
        - Manual page updated for rbio_unsupervised_hcluster(), along with its S3 class types
        - Fixed a bug for rbio_supervised_hcluster() where RNAseq data cannot be processed
      
      - Other fixes
        - Typos fixed for manual pages
        

    0.5.3 (Aug.6.2020)
      - General updates
        - match.arg() method added to relevant functions for better user experience
        - Manual pages combined for S3 methods
        - All functions updated to be compatible with R version 4.0
        - Citation added

      - Updates to microarray functions
        - A bug fixed for rbioarray_rlist() where the function crashes when no gene annotation data frame is provided
        - Information regarding the design colnames specification added for rbioarray_transfo_normalize()'s manual page
        - A bug fixed for sig() where the cutoff was set to "raw p value less than the FDR resulted p cutoff"

      - Updates to RNAseq functions
        - A bug fixed for rnaseq_de default method where it doesn't automatically load cpm() from edgeR
        - A bug fixed for rnaseq_de default method where the function has an unused arugment "normalization"
        - A bug fixed for rnaseq_de default method where the function fails when filtering resulted in singular (i.e. all filtered in or out) results
        - manual page updated for rnaseq_de()
        
      - Update(s) to cluster functions
        - rbio_supervised_hcluster() export file name suffix changed to "_sig_heatmap.pdf"

      - Updates to legacy functions
        - Suppressmessage method added for legacy functions figure export dev.off()
        - rbioarrary_DE() now also exports the DE summary data frame into the environment
        - A bug fixed for rbioarrary_DE() where the cutoff was set to "raw p value less than the FDR resulted p cutoff"
        - verbose argument added to rbioarrary_DE(), rbioarray_hclust_super
        - rbioarrary_DE() updated with suppport for continuous outcome
        - rbioarray_hcluster() and rbioarray_hcluster_super() can now hide top heatmap strip by setting ColSideCol = FALSE
        - rbioarray_hcluster_super() export file name suffix changed to "sig.pdf"
        
      - Other updates
        - Small fixes
        - Dependency ggplot2 now requires version 3.0.0
      
    
    0.5.2 (Jan.17.2019)
      - New clustering funcions:
        - Fixes made to rbio_supervised_hcluster() so that it skips the comparisons without significant result
        
      - Updates to shared functions:
        - Full DE gene-level stats are now stored in "sig" class
        - Docmentation updated for sig() function
      
      - Updates to RNAseq functions:
        - The "between.samples.norm.method" argument changed to "library.size.scale.metho" for rnaseq_de()
        - The "between-sample" item changed to "Library size-scaling" for the pirnt method for "rbioseq_de"
        - Print method for "rbioseq_de" now shows "Library size-scaling (between-sample)" normalization method first
        - Updated documentation for voom process for rnaseq_de()

      - Updates to microarray functions:
        - The gene.annot.dataframe argument changed to extra.gene.annot.dataframe for rbioarray_rlist()
        - Better doumentation for rbioarray_transfo_normalize(), rbioarray_filter_combine() and microarray_de()

      - Other updates
        - New Venn diagram function rbio_venn_de() for the S3 class "sig"
          - Due to the S3 revamp, the old venn diagram function rbioarray_venn_DE is now a legacy function
        - Package description updated
        - New bioconductor installation instructions added
      
      - Other bug fixes
       
    
    0.5.1
      - Updates to RNAseq functions
        - TMM scaling now correctly identified as "between-sample" normalization, whereas Voom as "between-gene" normalization method
      
      - Bug fixes
      
        
    0.5.0
      - New RNAseq functions
        - rbioseq_import_gtf() function to import and parse gtf/gff annotation files
        - rbioseq_import_count() function added to import read count files (e.g. HTseq-count files) to R environment. The function outputs an "rbioseq_count" object
        - rbioseq_clr_ilr_transfo() function added for transforming data for clustering analysis and any feature selection/classcification processes
        - New RNAseq differential expressin function rbioseq_de_analysis() with OOP elements - using S3 classes
          - The fnction has export options:
            (i) all features, 
            (ii) features with a annotation name, 
            (iii) significant features with or without annotation name depending on the setting for argument "gene_symbol".
        - rnaseq_de function added and produces "rbioseq_de" object. The function inlcudes methods for the following classes: "rbioseq_count", "mir_count"
        - Significance test function added and produces "sig" object
        - S3 print methods added for the classes: "rbioseq_de", "rbioseq_count", "sig"
        
      - New microarray functions
        - Function rbioarray_rlist() to import and produce a raw data class "rbioarray_rlist"
          - The function also supports "EListRaw" class objects from limma package
        - Function rbioarray_transfo_normalize() for log transforming and normalizing raw data from "rbioarray_rlist". The function produces an "rbioarray_plist" class object
        - Function rbioarray_filter_combine() for filtering, averaging and (if set) combining transcripts from the same gene/genomic feature. The function produces an "rbioarray_flist" object
        - S3 print method for "rbioarray_rlist", "rbioarray_plist", "rbioarray_flist" classes

      - New cluster functions
        - rbio_unsupervised_hcluster, with "rbioarray_flist" and "rbioseq_de" classes as input
        - rbio_supervised_hcluster, with "sig" class as input
        - rbio_unsupervised_corcluster, with "rbioarray_de" and "rbioseq_de" classes as input 
        - rbio_supervised_corcluster, with "sig" class as input

      - Updates to RNAseq functions:
        - Due to the overhual of RNAseq and microarray DE functions, rbioarray_DE() is now considered as a "legacy function". However, it is still functional for compatibility.
      
      - Updates to microarray functions:
        - Due to incorporating OOP functions, rbioarray_PreProc() is now considered as a "legacy function". However, it is still functional for compatibility
        - Due to incorporating OOP functions, rbioarray_flt() is now considered as a "legacy function". However, it is still functional for compatibility
        - Due to the overhual of microarray DE functions, rbioarray_DE() is now considered as a "legacy function". However, it is still functional for compatibility

      - Updates to cluster functions:
        - Due to S3 method implementation, all cluster functions are now legacy functions: rbioarray_hcluster(), rbioseq_hcluster(), rbioarray_hcluster_super(), rbioarray_corcluster_super(). They are still functional
        
      - Updates to legacy functions:
        - rbioarray_DE() now has the options to produce csv files for either 
          (i) all probes, 
          (ii) significant probes, 
          (iii) all probes with gene name, 
          (iv) significant probes with gene name. 
        However, DE reuslts for all probes will be exported to the R environment regardless of these settings. Similarly, F stats is also always exported to the working directory and the R environment regardless of these settings
        - rbioarray_DE() with the "FORK" cluster module re-written for foreach style parallel computing
        - To keep things consistent with limma's Elist, the "target" component list output changed to "targets" for all microarray functions
        
      - Verbose option added for all the functions
        
      - Wording adjustment for clearer documentation
      
      - Version bumped to 0.5.0
      
      - Bug fixes
      
      
    0.4.6
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
      

    0.1.1 - 0.1.34
      - Elist objects now can be properly recognized by all the functions
      - Help page for unsupervised heat map updated
      - Users can now choose to keep or remove the control probes (Agilent microarray platform) for rbioarray_hcluster() function
      - FC filter added for supervised hclust heatmap function
      - Volcano dots colours are now user customizable
      - Volcano dots annotation functionality added
      - KEGG visualization function updated
      - Code redundancy further reduced for array DE function
      - Array DE function re-written with better parallel computing (doParallel and foreach), 
      - Supervised hclust heatmap function updated for label display functionality
      - Supervised hclust heatmap function added
      - Stand alone GSA and GSA plotting functions added
      - All-in-one GSA function re-written
      - rbioGS_boxplot, rbioGS_scatter and rbioGS_all functions updated
      - GSA methods for GS functions are now user customizable
      - Preliminary all-in-one GS function added
      - GS plotting functions added
      - GS function re-written
      - Part of the array DE function re-written
      - Plot section of the array DE function improved
      - Log transformation for array functionality added
      - hclust heatmap function added
      - Volcano plot functionality added for array DE function
      - Fit and contrasts re-written for array DE function
      - Array DE function updated
      - Array DE function added
      - Preliminary microarray functions added
      - Bug fixes
      
      
    0.1.0
      - Initial seq functions added
      - Initial GS functions added
      - Initial commit
