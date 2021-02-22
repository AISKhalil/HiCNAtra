classdef HiCNAtra < matlab.mixin.Copyable
	properties (GetAccess=public, SetAccess=public)

		%%%% -------------- Input-data/ Main-parameters -----------------%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%HDF5Files: Input HDF5 files (a string array ex: {'/input/file1.hdf5','/input/file2.hdf5'}) that include the mapped-reads (chromosomes, cut-sites,
        %restriction-fragments IDs, strands). The expected files are the hiclib mapped-reads after running parsing-sam method.
		HDF5Files;

		%readLength: The length of side/end read.
		readLength;

		%restrictionEnzyme: An enzyme that cleaves DNA into fragments at or near specific recognition sites within the molecule known as restriction sites.
		restrictionEnzyme;

		%maximumMoleculeLength: Maximum length of molecules in the HiC library, used as a cutoff for dangling ends filter.
		maximumMoleculeLength = 500;
		
		%referenceGenome: The name of the reference genome {'hg19','hg18','hg38','mm9','mm8'}. It is used for selecting the annotations folder and defining the chromosomes names.
		referenceGenome = 'hg19';
		
		%binSize: Number of bases per bin (5Kb is the default bin size which is of the same order of magnitude as the experimental resolution of 6-bp cutter, average 4096 bases for 6-bp cutter.
		binSize = 5000;
	
		%contactMapBinSize: bin-size for computing the contact-map.
		contactMapBinSize = 100000;
	

		%%%% -------------------- Annotation Files ---------------------%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		%referenceGenomeFolder: The folder that contains UCSC sequences files of the reference genome (.Fa files) that is used to calculate the effective-length feature.
		%A default hg19 folder: ".../gcWinds.readLength.100/"		
		referenceGenomeFolder;

		%mappabilityFolder: The folder that contains Anshul's uniquely mappability files of the reference genome (.unqiue files)  which can be used to filter highly repeated or unmappable regions. 
		%A default hg19-based folder ".../globalmap_k20tok81/" for read-length < 100bp.
		%A default hg19-based folder ".../globalmap_k101tok101/" for read-length >= 100bp.
		mappabilityFolder;

		%gcWindsFolder: Folder that contains the Christopher A. Miller's pre-calculated gc-contents per 100 bps at different read lengths.
		%A default hg19 folder: ".../gcWinds.readLength100/"
		gcWindsFolder;

		%blackListFile: File that includes black-listed regions. 
		%A default hg19-based file "...Anshul_wgEncodeHg19ConsensusSignalArtifactRegions.bed".
		blackListFile;

		%gapFile: File that includes unmappable hg19 regions. A default hg19-based file "...UCSC_gapRegions.txt".		
		gapFile;

		%centromeresFile: File that includes locations of centromeres. 
		%A default hg19 file "...UCSC_Centromeres.txt".
		centromeresFile;

		%telomeresFile: File that includes locations of telomeres. 
		%A default hg19 file "...UCSC_Telomeres.txt".
		telomeresFile;		

		%outputDirectory: Output directory to save results and figures at.
		outputDirectory = './HiCNAtraOutput';


		%%%% ----------------------- Parameters -------------------------%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%
		%% RD-calculation %%
		%%%%%%%%%%%%%%%%%%%%
		%RDmethod: Method that is used for computing the RD signal. 1) "entire restriction fragment" counting (best for Hi-C data), 2) Paired-end method(best for 3C-seq), 3) Exact-cut position, 4) Mid restriction-fragment mapping.
		RDmethod = 1;

		%gcCalculationMethod: (1) the Christopher A. Miller's pre-calculated gc-contents, (2) reference genome sequence.
		gcCalculationMethod = 2;

		%memoryFootPrint: Number of reads that can be loaded per iteration for informative and non-informative reads. User can choose it based on RAM size.
		memoryFootPrint = 1000000;
		

		%%%%%%%%%%%%%%%%%%%%
		%% Bin-filtering  %%
		%%%%%%%%%%%%%%%%%%%%		
		%minMappabilityThreshold: Minimum mappability-score of a bin to be kept for HiCNAtra analysis.
		minMappabilityThreshold = 0.5;
	
		%minGCThreshold: Minimum GC-score of a bin to be kept for HiCNAtra analysis.
		minGCThreshold = 0.25;

		%minEffectiveLength: minimum effective-length fraction of the bin to be kept for normalization.
		minEffectiveLength = 0.1;		

				
   		%%%% -------------------- Annotation Data ---------------------%%
   		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%rsites: Dictionary of restriction-sites per chromosome.
		rsites;

		%chrEffectiveLengthTracks: Dictionary of effective-length tracks per chromosome.
		chrEffectiveLengthTracks;

		%mappabilityWindowScores: Dictionary of mappability-scores for each GC-window per chromosome.
		mappabilityWindowScores;

		%chrMappabilityTracks: Dictionary of mappability tracks per chromosome (window-based for Hi-C, whole-bin for 3C-seq).
		chrMappabilityTracks;		

		%GCWindsWindowScores: Dictionary of GC-scores for each GC-window per chromosome.
		GCWindsWindowScores;

		%chrGCTracks: Dictionary of GC-tracks per chromosome.
		chrGCTracks;

		
		%%%% -------------------- ReadDepth Data ---------------------%%	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
		%chrDictionary: ReadDepth signal per chromosome before filtering and normalization.
		chrDictionary;

		%chrFDictionary: ReadDepth signal per chromosome after filtering gap-regions, black-listed, centromeres, telomeres, low-mappability regions.
		chrFDictionary;

		%chrFIndex: Filtered indices that can be used to filter-out RD signals.
		chrFIndex;

		%chrFNDictionary: ReadDepth signal per chromosome after normalization.
		chrFNDictionary;

		
		%%%% --------------------- CNA Computing ---------------------%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%ploidyLevel: Expected number of chromosome sets (default = 2), user can set it based on input-cell ploidy {'free', 'diploid', 'triploid', 'tetraploid'}.
		ploidyLevel = 'free'; 

		%minimumFAsize: A minimum segment size to be considered as focal alterations, default size = 100Kb.
		minimumFAsize = 100000;
        		
		%minimumIBsize: A minimum segment size to be considered as Iso-copy numeric block, default size = 1Mb.
		minimumIBsize = 1000000;

		%resolution: Minimum alteration-region that can be detected by HiCNAtra (interms of bins). It is computed using regression model based on the data coverage. We used it as a frame length for Savitzky Golay filter.
		resolution;

		%amplificationThreshold: Threshold for identifying the class2 amplifications.
		amplificationThreshold;

		%deletionThreshold: Threshold for identifying the class2 deletions.
		deletionThreshold;

		%%%%%%%%%%%%%%%%%
		%- CNV results -%
        %%%%%%%%%%%%%%%%%		
		%CNReference: The copy number reference of the genome.
		CNReference;

		%segmentsInfoDic: Dictionary of IBs per chromosome.
		segmentsInfoDic;

		%regionsInfoDic: Dictionary of CNVs per chromosome.
		regionsInfoDic;

		%copy-number tracks.
		chrCNVsTracks;
		
		%%%% ---------------- HiCNAtra Normalization -------------------%%
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    %cisOnly: flag to choose which matrices to compute & normalize
	    cisOnly = 1;
	    
		%cisTransNormaMethod: flag to choose the normalization method for cis-data: 1: genome-wise, 2: chromosome-wise.
		cisTransNormMethod = 1;
		
		%minInteractionFreq: minimum interaction-frequency between two-bins to be used for GLM-normalization.
		minInteractionFreq = 1;

		%contactMapEffectiveLengthTracks: 
		contactMapEffectiveLengthTracks;

		%contactMapGCTracks:
		contactMapGCTracks;

		%contactMapMappabilityTracks:
		contactMapMappabilityTracks;

		%contactMapCNVsTracks:
		contactMapCNVsTracks;


		%%%% ---------------------- Data Statistics -------------------%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%readStatistics: Number of reads for each category [total-reads, double-side, informative, self-circle, dangling-end, extraDangling-end, single-side];
		readStatistics;

		%dataCoverage: Genome coverage of the dataset = (# base-pairs covered by the restriction-fragments/the dangling-end fragments /genome-length(bps)).
		dataCoverage;		

		%noValidReadPairs: number of valid read-pairs.
		noValidReadPairs;
		
	end	
	%%%
	%%%
	%%%
	%%%
	%properties (GetAccess=private, SetAccess=private)
	properties (GetAccess=private, SetAccess=private)	
   		%%%% ------------------- Annotation Data --------------------%%
   		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%readTypes: Types of reads that are included to compute the RD signal of Hi-C data. 1) All-reads types, 2) Info-reads, 3) Non-info reads, and 4) single-sided reads
		readTypes = 1;

		%gcWindow: Half length of the window surrounding the restriction-site that is used to compute the gc-content feature.
		gcWindow = 200;
		
		%effectiveLengthWindow: Half length of the window surrounding the restriction-site that is used to compute the effective-length feature (ELW).
		%Intially, effectiveLengthWindow = maximumMoleculeLength.
		effectiveLengthWindow;

		%mappabilityWindow: Half length of the window surrounding the restriction-site that is used to compute the mappability feature.
		%Intially, mappabilityWindow = maximumMoleculeLength.
		mappabilityWindow;

		%cutToRsiteThreshold: Distance to restriction-site that is used to filtering the reads.
		cutToRsiteThreshold;

		%targetChrs: Target chromosomes for the HiCNAtra analysis. 0: All chromosomes, 23: chromosome X.
		targetChrs = 0;
		
		%removeBlackBins: Flag to remove black bins.		
		removeBlackBins = 1;	

		%removeLowMappabilityBins: Flag to remove low-mappability bins.		
		removeLowMappabilityBins = 1;

		%removeLowGCBins: Flag to remove low-GC bins.		
		removeLowGCBins = 1;	

		%largeRfragmentWidth: Minimum width to define large restriction fragments.		
		largeRfragmentWidth = 100000;	

		%uncertaintyDistanceToCentro: Area surrounding centromeres , with normal RD-value to extend the centromeres and telomeres boundaries (default size: 100Kb). 
		uncertaintyDistanceToCentro = 100000;
		
		%uncertaintyDistanceToTelo: Area surrounding telomeres, with normal RD-value to extend the centromeres and telomeres boundaries (default size: 10Kb). 
		uncertaintyDistanceToTelo = 10000;

		%chrLargeFragments: Dictionary of large-fragments per chromosome.
		chrLargeFragments;
		
		%chrCentroTeloBoundaries: Dictionary of extended telomeres and centromeres per chromosome.
		chrCentroTeloBoundaries;

		%chrBlackListedBoundaries: Dictionary of extended black-listed per chromosome.
		chrBlackListedBoundaries;

		%chrGapBoundaries: Dictionary of gap regions per chromosome.
		chrGapBoundaries;	
		
		%CNAtraRdNormalization: 1) use CNAtra interval-normalization for RD correction, 0) use GLM regression-model normalization.
		CNAtraRdNormalization = 1;
		
		%toolID: for testing contact-map normalization approaches 1) HiCNAtra, 2) HiCnorm, 3) OneD, 4) oneD+CN 
		toolID = 1;

		%%%% -------------------- ReadDepth Data ---------------------%%	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%chrNames: chromosome names.
		chrNames;
		
		%chrLengths: chromosome lengths.
		chrLengths;

		%genomeNormRD: ReadDepth signal of the genome after normalization and multiplied by 2/CNReference.
		genomeNormRD;

		%chrStartIndex: Indices of chromosome start bins for each chromosome. They are used for extract chromosome normalized RD-signal from genome normalized RD-signal.
		chrStartIndex;


		%%%% ---------------------- CNV calling ----------------------%%	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%maximumFalseBinsAllowed: Maximum allowed region percentage to keep alteration regions
		maximumFalseBinsAllowed = 0.5;	

		%%%% ---------------------- Data Statistics -------------------%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		%NDupReadStatistics: Number of reads after removing non-duplicated reads [total-reads, informative, non-informative, single-side];
		NDupReadStatistics;

		%nearRSitesSides: Number of reads after removing non-duplicated reads [total sides, informative, non-informative, single-side];
		nearRSitesSides;
	end
	%%%
	%%%
	%%%
	%%%
	methods
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		function obj = HiCNAtra(hiclibFiles, HiC_Directory, readLength, restrictionEnzyme, maximumMoleculeLength, referenceGenome)
		%% A pipeline for CNV detection and contact-map normalization of Hi-C/3C-Seq data.

			if nargin == 6
				obj.HDF5Files         = hiclibFiles;
				obj.readLength        = readLength;
				obj.restrictionEnzyme = restrictionEnzyme;
				obj.referenceGenome   = referenceGenome;
				%
				obj.referenceGenomeFolder    = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/UCSC_chromFa/');
				%
				if(readLength >= 100)
					obj.mappabilityFolder    = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/Anshul_UniqueMappability/globalmap_k101tok101/');
				else
					obj.mappabilityFolder    = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/Anshul_UniqueMappability/globalmap_k20tok81/');
				end
				%
				if(readLength >= 200)
					obj.gcWindsFolder        = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/ChrisaMiller_GCContents/gcWinds.readLength200/');				
				elseif(readLength >= 100) 
					obj.gcWindsFolder        = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/ChrisaMiller_GCContents/gcWinds.readLength100/');	
				elseif(readLength >= 76) 
					obj.gcWindsFolder        = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/ChrisaMiller_GCContents/gcWinds.readLength76/');	
				elseif(readLength >= 50) 
					obj.gcWindsFolder        = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/ChrisaMiller_GCContents/gcWinds.readLength50/');						
				elseif(readLength >= 36) 
					obj.gcWindsFolder        = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/ChrisaMiller_GCContents/gcWinds.readLength36/');		
				else
					obj.gcWindsFolder        = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/ChrisaMiller_GCContents/gcWinds.readLength27/');						
				end
				
				%				
				obj.blackListFile	         = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/Anshul_wgEncodeHg19ConsensusSignalArtifactRegions.bed');
				obj.centromeresFile          = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/UCSC_Centromeres.txt');	
				obj.telomeresFile    	     = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/UCSC_Telomeres.txt');
				obj.gapFile                  = strcat(HiC_Directory, '/Annotations/', referenceGenome, '/UCSC_gapRegions.txt');

				%
				%
				obj.maximumMoleculeLength = maximumMoleculeLength;
				obj.effectiveLengthWindow = obj.maximumMoleculeLength;
				obj.mappabilityWindow     = obj.maximumMoleculeLength;
				obj.cutToRsiteThreshold   = obj.maximumMoleculeLength;

				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%%%%%%%%%%%%%%  Create Empty-dictionaries   %%%%%%%%%%
				%% prasing the bias-sources %%
				obj.chrLengths = containers.Map({1},{[]});
				remove(obj.chrLengths,1);
				%
				obj.rsites = containers.Map({1},{[]});
				remove(obj.rsites,1);        	 
				%	
				obj.GCWindsWindowScores = containers.Map({1},{[]});
				remove(obj.GCWindsWindowScores,1);
				%
				obj.mappabilityWindowScores = containers.Map({1},{[]});
				remove(obj.mappabilityWindowScores,1);
				%%%		
				%%%
				%%%
				%% computing the bias-tracks based on the bin-size %%
				obj.chrEffectiveLengthTracks = containers.Map({1},{[]});
				remove(obj.chrEffectiveLengthTracks,1);   
				%
				obj.chrMappabilityTracks = containers.Map({1},{[]});
				remove(obj.chrMappabilityTracks,1);   
				%
				obj.chrGCTracks = containers.Map({1},{[]});
				remove(obj.chrGCTracks,1);				
				%%%
				%%%
				%%%
				%% extracting the dark-regions %%
				obj.chrCentroTeloBoundaries = containers.Map({1},{[]});
				remove(obj.chrCentroTeloBoundaries,1);		
				%
				obj.chrBlackListedBoundaries = containers.Map({1},{[]});
				remove(obj.chrBlackListedBoundaries,1);		
				%
				obj.chrGapBoundaries = containers.Map({1},{[]});
				remove(obj.chrGapBoundaries,1);		
				%
				obj.chrLargeFragments = containers.Map({1},{[]});
				remove(obj.chrLargeFragments,1);				
				%%%
				%%%
				%%%
				%% computing and normalization the RD-signal %%
				obj.chrDictionary = containers.Map({1},{[]});
				remove(obj.chrDictionary,1);		
				%
				obj.chrFDictionary = containers.Map({1},{[]});
				remove(obj.chrFDictionary,1);	
				%
				obj.chrFIndex = containers.Map({1},{[]});
				remove(obj.chrFIndex,1);
				%
				obj.chrFNDictionary = containers.Map({1},{[]});
				remove(obj.chrFNDictionary,1);		
				%%%
				%%%
				%%%
				%% detecting the LCVs/FAs of the RD-signal %%
				obj.segmentsInfoDic = containers.Map({1},{[]});
				remove(obj.segmentsInfoDic,1);
				%
				obj.regionsInfoDic = containers.Map({1},{[]});
				remove(obj.regionsInfoDic,1);
				%
				obj.chrCNVsTracks = containers.Map({1},{[]});
				remove(obj.chrCNVsTracks,1);  								
				%%%
				%%%
				%%%
				%%% normalization of contact map %%
				obj.contactMapEffectiveLengthTracks = containers.Map({1},{[]});
				remove(obj.contactMapEffectiveLengthTracks,1);
				%
				obj.contactMapGCTracks = containers.Map({1},{[]});
				remove(obj.contactMapGCTracks,1);
				%
				obj.contactMapMappabilityTracks = containers.Map({1},{[]});
				remove(obj.contactMapMappabilityTracks,1);
				%
				obj.contactMapCNVsTracks = containers.Map({1},{[]});
				remove(obj.contactMapCNVsTracks,1);
				%%%
				%%%
				%%%
			else
				error('No enough inputs')
			end
			%%% End-if %%%
		end   
		%%% End-constructor %%%





		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%% -------- Input Handling ------- %%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 1) Parsing the GC/Mappability/reference genome files.
		% 2) computing the read-depth signal from informative and non-informative reads, 
		% 3) bias-correction of th RD-signal.
		biasExtraction(obj)
		RDcalculator(obj)

		%--Biases-tracks--%
		restrictionSites(obj)
		effectiveLengthFeatureCalculater(obj)
		%
		mappabilityWindowsScores(obj)		
		mappabilityFeatureCalculater(obj)
		mappabilityRfragFeatureCalculater(obj)
		mappabilityBinFeatureCalculater(obj)
		%
		gcWindowsScores(obj)
		gcFeatureCalculater(obj)
		gcRfragFeatureCalculater(obj)
		gcBinFeatureCalculater(obj)
		%
		%
		%--Read-Depth Calculation--%
		readsExtractionExactCutPositions(obj)
		readsExtractionRestrictionFragmentsMid(obj)	
		readsExtractionRestrictionFragments(obj)
		readsExtractionDanglingFragments(obj)
		%
		darkRegionsExtraction(obj)
		%
		%
		%--Read-Depth Normalization--%
		RDfiltering(obj)
		RDnormalization(obj)
		pipelineParameters(obj)

		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%% ---- Copy-number Computing ----%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 1) copy-number reference computing, 
		% 2) Savitzky-Golay segmentation, 
		% 3) segments merging, 
		% 4) CNA calling.
		CNVcaller(obj)

		copyNumberReference(obj)
		[shortSegmentStart, shortSegmentEnd, longSegmentStart, longSegmentEnd] = rdSavitzkySegmentation(obj, targetChrIndex)
		[finalSegmentStart, finalSegmentEnd] = shortSegmentsMerging(obj, targetChrIndex, orgSegmentStart, orgSegmentEnd)
		[segmentsInfo, regionsPerSegment] = CNRegionCalling(obj, targetChrIndex, orgSegmentStart, orgSegmentEnd, finalSegmentStart, finalSegmentEnd)
		%
		[finalSegmentInfo, impRegionsPerSegment] = CNRegionFiltering(obj, chrIndex, segmentsInfo, regionsPerSegment)
		%
		[segmentsData, regionsData] = CNResultsUpdate(obj, chrIndex, segmentsInfo, impRegionsPerSegment)
		CNResultsUpdateGenome(obj)
		computeEventsTracks(obj)
		%
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%% -- contact-map normalization --%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 1) computing the 3D contact-map,
		% 2) normalizing the contact-map. 
		contactMapCorrector(obj)
		%
		contactMapComputing(obj)
		contactMapBiasExtraction(obj)
		contactMapNormalization(obj)
		cisGenomeWiseNormalization(obj)
		transGenomeWiseNormalization(obj)
		cisTransChromosomeWiseNormalization(obj)
		%
		HiCNormNormalization(obj)
		OneDNormalization(obj)
		OneDplusCnNormalization(obj)
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%% ------- Visualization ------- %%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% Visual test of the assumed plooidy-level
		ploidyTest(obj)
		% Plot chromosome or a specific area.
		CNAPlot(obj, saveResult, varargin)
		% Plot amplification or deletion region by its number.
		CNARegionPlot(obj, saveResult, varargin)
		% Plot the CNVs tracks.
		CNVsTrackPlot(obj, saveResult, varargin) 
		%
		plotGenome(obj, saveResult, nBinSize)
		
		% Plot raw interaction-matrix.
		rawContactMapPlot(obj, saveResult, varargin) 
		% Plot normalized interaction-matrix.
		normContactMapPlot(obj, saveResult, varargin)
		% Plot Genome-wide interaction-matrix
		rawGenomeContactMapPlot(obj)
		normGenomeContactMapPlot(obj)

		% Reads analysis.	
		plotFragmentLengthHistogram(obj)
		[correlationValues] = RDvsBiasesCorrelation(obj)
		[distancePercentages] = rsitesDistancesAnalysis(obj,plotCond)
		[gainInfo] = nonInfoGain(obj)
		
	
		end 
		%%%
		
end   
%%%
%%%
%%%