# git_RBioArray
Simple to use package for both microarray and RNAseq data analysis (limma based)

Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        source("https://bioconductor.org/biocLite.R")
      
        biocLite()
    
  - Install the package
  
        devtools::install_github("jzhangc/git_RBioArray/RBioArray", repos = BiocInstaller::biocinstallRepos())   



Update log

    0.2.6 (March.26.2017)
      - Both PSOCK and FORK cluster types are now avaiable for all the functions featuring parallel computing
      - Small "quality of life" changes
      - Bug fixes
    

    0.2.5
      - Direct loading GS database files functionality added for GS functions
      - Bug fixes
      
    
    0.2.4
      - mmu/rno to hsa entrez ID converesion function overhaul
      - Arguments for input DE data list for GS functions unified as "DElst"
      - Arguments for parallel computing unified for all the relevant functions as "parallelComputing" and "clusterType"
      - Bug fixes
      

    0.2.2 - 0.2.3 
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
      
    0.1.24
      - All-in-one GSA function re-written
      - Bug fixes
      
    0.1.23
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
      
    0.1.18
      - GS plotting functions added
      - Bug fixes
      
    0.1.17
      - Bug fixes
      
    0.1.16
      - GS function re-written
      - Bug fixes
      
    0.1.15
      - Bug fixes
      
    0.1.14
      - Part of the array DE function re-written
      - Bug fixees
      
    0.1.13
      - Plot section of the array DE function improved
      - Bug fixes
      
    0.1.12
      - Log transformation for array functionality added
      - hclust heatmap function added
      - Bug fixes
      
    0.1.10 - 0.1.11
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
      - Initial release
    
