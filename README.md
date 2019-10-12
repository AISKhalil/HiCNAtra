# HiCNAtra: An analytical and visualization tool for CNV discovery and contact map correction of Hi-C/3C-seq data of cancer cell lines. 

**HiCNAtra** is a MATLAB-based tool that accepts HDF5 files, the output of [hiclib](https://mirnylab.bitbucket.io/hiclib/index.html?) after applying the iterative-mapping technique, as input. **HiCNAtra** pipeline is divided into three modules: 1) computation of the read depth (RD) signal from Hi-C or 3C-seq reads (RD calculator), 2) RD-based detection of copy number events based on [CNAtra](https://github.com/AISKhalil/CNAtra) approach (CNV caller) and 3) bias correction of chromatin interaction matrix introduced by CNVs and other systematic biases (Contact map normalization).  
**HiCNAtra** generates many output files providing the detailed characterization of the copy number profile, the raw contact map, and the HiCNAtra-corrected contact map of the input Hi-C/3C-seq data. It saves BED format files of both large-scale copy number variations (LCVs) and focal alterations (FAs) that can be uploaded into UCSC Genome Browser. Also, it saves text files of the cis and trans interaction frequencies before and after HiCNAtra correction.  
**HiCNAtra** also provides a visual platform that allows manual review of the CNV profile and the contact maps pre- and post-normalization. The detailed description of HiCNAtra inputs, parameters, methods, and outputs will be provided in **"HiCNAtra_User_Guide.pdf"** file.  
**For a full description of the method and applications, please visit [HiCNAtra Manuscript](https://www.biorxiv.org/content/10.1101/798710v1).**
  
## Contents
- [Download](#Download)
- [Directory Setup](#directory_setup)
- [Annotations](#annotations)
- [Parameters](#parameters)
- [Input Preparation](#input_preparation)
- [Usage](#usage)

### <a name="Download"></a>Download
```bash
cd ~
git clone https://github.com/AISKhalil/HiCNAtra.git
```

### <a name="directory_setup"></a>Directory Setup
After downloading the **HiCNAtra** directory, `HiCNAtra` annotation files (reference genome sequence, mappability tracks, and GC tracks) should be downloaded and allocated to their corresponding sub-directories inside the **HiCNAtraTool** directory:
- The annotations directory structure will look like this:

```
    HiCNAtraTool/
    +- @HiCNAtra/
    +- Annotations/
    |  +- hg19/
    |  |  +- UCSC_chromFa/
    |  |  |  +- chr1.fa
    |  |  |  +- chr2.fa
    |  |  |  +- . . .
    |  |  |
    |  |  +- Anshul_UniqueMappability/
    |  |  |  +- globalmap_k20tok81/
    |  |  |  |  +- chr1.uint8.unique
    |  |  |  |  +- chr2.uint8.unique    
    |  |  |  |  +- . . .
    |  |  |  |
    |  |  |  +- globalmap_k101tok101/    
    |  |  |  +- . . .    
    |  |  |
    |  |  +- ChrisaMiller_GCContents/
    |  |  |  +- gcWinds.readLength.100/
    |  |  |  |  +- chr1.gc
    |  |  |  |  +- chr2.gc
    |  |  |  |  +- . . .    
    |  |  |  |
    |  |  |  +- gcWinds.readLength.50/
    |  |  |  +- . . .    
    |  |  +- UCSC_Centromeres.txt
    |  |  +- UCSC_Telomeres.txt
    |  |  +- UCSC_gapRegions.txt
    |  |  +- Anshul_wgEncodeHg19ConsensusSignalArtifactRegions.bed
    |  |  
    |  +- hg18/
    |  +- hg38/
    |  +- . . .
```

### <a name="annotations"></a>Annotations
**1)** **HiCNAtra** uses the reference genome sequence (e.g. hg19) for computing the effective lengths that are used for correcting the Hi-C/3C-seq contact map based on the experiment (restriction enzyme).
Please download and extract the reference genome sequence [hg19 reference genome sequence](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz) in `...HiCNAtra/HiCNAtraTool/Annotations/hg19/UCSC_chromFa/` sub-directory. 
  
**2)** **HiCNAtra** also uses the unique mappability tracks for computing thr mappability scores that are used for correcting the Hi-C/3C-seq contact map  [Unique mappability tracks for several species](https://sites.google.com/site/anshulkundaje/projects/mappability). This includes per-base unique mappability tracks for a large range of read lengths for several key species [Umap and Bismap: quantifying genome and methylome mappability](https://academic.oup.com/nar/article/46/20/e120/5086676). Please download and extract the mappability tracks [globalmap_k101tok101](https://personal.broadinstitute.org/anshul/projects/umap/encodeHg19Male/globalmap_k101tok101.tgz) and [globalmap_k101tok101](https://personal.broadinstitute.org/anshul/projects/umap/encodeHg19Female/globalmap_k20tok81.tgz) in `...HiCNAtra/HiCNAtraTool/Annotations/hg19/Anshul_UniqueMappability/` sub-directory.

**3)** **HiCNAtra** uses the GC tracks for computing the GC scores that are used for correcting the Hi-C/3C-seq contact maps [GC tracks](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/index.html). This includes the GC tracks for a large range of read lengths for human genome [ReadDepth](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016327). Please download and extract the GC tracks based on the read length [gcWinds.readLength100.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength100.hg19.tar), [gcWinds.readLength200.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength200.hg19.tar), [gcWinds.readLength76.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength76.hg19.tar), [gcWinds.readLength50.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength50.hg19.tar), [gcWinds.readLength36.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength36.hg19.tar), and [gcWinds.readLength27.hg19](https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength27.hg19.tar) in `...HiCNAtra/HiCNAtraTool/Annotations/hg19/ChrisaMiller_GCContents/` sub-directory.
  

### <a name="parameters"></a>Parameters
The main analysis parameters of **HiCNAtra**:

    HDF5Files                - the input HDF5 files {'/input/file1.hdf5','/input/file2.hdf5', ..}. 
                               This HDF5 files are the output of Hiclib tool after running the iterative alignning module.
   
    readLength               - the short sequencing read length.
 
    restrictionEnzyme        - the name of the restriction-enzyme that is used for the Hi-C/3C-seq experiment.

    maximumMoleculeLength    - the maximum molecule length (in bps). 

    referenceGenome          - the reference genome (e.g. hg19).
      
    binSize                  - the bin size of the RD signal(default = 5Kb).

    contactMapBinSize        - the bin size of the contact map (default = 100 Kb).

    outputDirectory          - the directory that is used for save all CNV information, raw contact map, and corrected contact map.

    RDmethod                 - the method to be used for computing the RD signal. 1) "entire restriction fragment" counting 
                               (best for Hi-C data), 2) Paired-end method (best for 3C-seq), 3) Exact-cut position, 4) Midpoint 
                               approach (default = 1).
    
    cisOnly                  - a flag to compute and normalize only the cis interaction frequencies (default = 1).
 
### <a name="input_preparation"></a>Input Preparation
**HiCNAtra** input is HDF5 files that are generated by [hiclib](https://mirnylab.bitbucket.io/hiclib/index.html?) after applying the [iterative mapping](https://mirnylab.bitbucket.io/hiclib/tutorial/01_iterative_mapping.html) module only. They include the information needed for **HiCNAtra** in a dict-like structure with the main keys `'chrms1', 'chrms2', 'cuts1', 'cuts2', 'rfragIdxs1', 'rfragIdxs2', 'strands1', 'strands2','rsites1','rsites2'`.

   **(1)** install the [hiclib](https://mirnylab.bitbucket.io/hiclib/index.html?).
   **(2)** edit the iterative mapping module  based on the read length and restriction enzyme information.
   and run the iterative mapping [Mapping.py](./Scripts/Mapping.py) <!-- (https://mirnylab.bitbucket.io/hiclib/tutorial/01_iterative_mapping.html) -->


after applying the iterative-mapping technique
. . .
### <a name="usage"></a>Usage
. . .

 
## HiCNAtra software... Coming soon!  
