# HiCNAtra: An analytical and visualization tool for CNV discovery and contact map correction of Hi-C/3C-seq data of cancer cell lines. 

`HiCNAtra` is a MATLAB-based tool that accepts HDF5 files, the output of [hiclib](https://mirnylab.bitbucket.io/hiclib/index.html?) after applying the iterative-mapping technique, as input. `HiCNAtra` pipeline is divided into three modules: 1) computation of the read depth (RD) signal from Hi-C or 3C-seq reads (RD calculator), 2) RD-based detection of copy number events based on [CNAtra](https://github.com/AISKhalil/CNAtra) approach (CNV caller) and 3) bias correction of chromatin interaction matrix introduced by CNVs and other systematic biases (Contact map normalization).  
`HiCNAtra` generates many output files providing the detailed characterization of the copy number profile, the raw contact map, and the `HiCNAtra`-corrected contact map of the input Hi-C/3C-seq data. It saves BED format files of both large-scale copy number variations (LCVs) and focal alterations (FAs) that can be uploaded into UCSC Genome Browser. Also, it saves text files of the cis and trans interaction frequencies before and after HiCNAtra correction.  
`HiCNAtra` also provides a visual platform that allows manual review of the CNV profile and the contact maps pre- and post-normalization. The detailed description of `HiCNAtra` inputs, parameters, methods, and outputs will be provided in **"HiCNAtra_User_Guide.pdf"** file.  
**For a full description of the method and applications, please visit [HiCNAtra Manuscript](https://www.biorxiv.org/content/10.1101/798710v1).**
  
## Contents
- [Download](#Download)
- [Directory Setup](#directory_setup)
- [Annotations](#annotations)
- [Parameters](#parameters)
- [Usage](#usage)

## <a name="Download"></a>Download:









## HiCNAtra software... Coming soon!  
