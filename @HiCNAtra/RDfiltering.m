function obj = RDfiltering(obj)
%Removing the low-mappability, large-fragments, black-listed, centromeres, telomeres, and gap regions.


%%---------------------------- Input Data -----------------------------%%
chrs = obj.targetChrs;
% chrs == 0: all chromosomes
% For human, chrs == 23: chrX
% For mouse, chrs == 19: chrX

if (chrs == 0) 
    chromosomes = [1:1:length(obj.chrNames)];
else
    chromosomes = chrs;
end
noChrs = length(chromosomes);
%
chrSelectedDictionary = obj.chrDictionary;


%%------------------------ Regions & Parameters -----------------------%%
chrLengths = obj.chrLengths;
binSize = obj.binSize;
%
chrMappabilityTracks = obj.chrMappabilityTracks;% dictionary of mappabiltiy per chromosome: mappability values per bins.
removeLowMappabilityBins = obj.removeLowMappabilityBins;
mappabilityThreshold = obj.minMappabilityThreshold;
%
chrGCTracks = obj.chrGCTracks;% dictionary of GC-ratios per chromosome: GC-scores per bins.
removeLowGCBins = obj.removeLowGCBins;
GCThreshold = obj.minGCThreshold; 
%
chrBlackListedBoundaries = obj.chrBlackListedBoundaries;% dictionary of black-listed regions per chromosomes.
removeBlackBins = obj.removeBlackBins;
%
chrCentroTeloBoundaries = obj.chrCentroTeloBoundaries;% dictionary of centromeres/telomeres per chromosomes.
chrGapBoundaries = obj.chrGapBoundaries;% dictionary of gap regions per chromosome.
%
rsites = obj.rsites;%dictionary of restriction-sites.
largeRfrags = obj.chrLargeFragments;%dictionary of large restriction-fragments.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------------------------- Filtering --------------------------%%
for i  = 1:1:noChrs 
	chrIndex = chromosomes(i);
	chrLengthBps  = chrLengths(chrIndex);
	chrLengthBins = ceil(chrLengthBps/binSize);
	selectedBins  = ones(chrLengthBins,1);

	
	%%%%%% centromere/telomeres %%%%%%%   
	centroTeloBoundaries = chrCentroTeloBoundaries(chrIndex);
	selectedBins(centroTeloBoundaries(1,1):centroTeloBoundaries(1,2))= 0;
	selectedBins(centroTeloBoundaries(2,1):centroTeloBoundaries(2,2))= 0;
	selectedBins(centroTeloBoundaries(3,1):centroTeloBoundaries(3,2))= 0;

	
	%%%%%%% gap regions %%%%%%%   
	gapBoundaries = chrGapBoundaries(chrIndex);
	[noGapRegions, ~]  = size(gapBoundaries);
	for j=1:noGapRegions
		selectedBins(gapBoundaries(j,1):gapBoundaries(j,2))= 0;    	
	end

	
	%%%%%%% Large fragments %%%%%%%%%
	chrRfragStart = [1; rsites(chrIndex)];
	chrRfragEnd   = [rsites(chrIndex)-1; chrLengthBps];
	%
	RfragStartBin = max(floor(chrRfragStart/binSize),1);
	RfragEndBin   = min(ceil(chrRfragStart/binSize),chrLengthBins);		
	%
	chrLargeRfrag = largeRfrags(chrIndex);
	noLargeRfrag  = length(chrLargeRfrag);
	for j=1:noLargeRfrag
		chrLargeRfragIndex = chrLargeRfrag(j);
		selectedBins(RfragStartBin(chrLargeRfragIndex):RfragEndBin(chrLargeRfragIndex))= 0;    	
	end

	
	%%%%% black-listed regions %%%%%%
	if(removeBlackBins)
	    blackListedBoundaries = chrBlackListedBoundaries(chrIndex);
	    [noBlackRegions,~] = size(blackListedBoundaries);
	    for j=1:noBlackRegions
			selectedBins(blackListedBoundaries(j,1):blackListedBoundaries(j,2))= 0;    	
		end	
	end   

	
	%%%%% low-mappability regions %%%%%
	if(removeLowMappabilityBins)
		mappabilityTracks = chrMappabilityTracks(chrIndex);
		mappabilityBins = filterLargeBins(obj.RDmethod, mappabilityTracks, mappabilityThreshold);
		selectedBins = bitand(selectedBins, mappabilityBins);
	end 	

	
	%%%%%%%%%% Low GC-ratios %%%%%%%%%%
	if(removeLowGCBins)
		GCTracks = chrGCTracks(chrIndex);
		GCBins   = filterLargeBins(obj.RDmethod, GCTracks,GCThreshold);
		selectedBins = bitand(selectedBins, GCBins);
	end	

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Filtering White/Black Regions %%%
	whiteRegions = find(selectedBins==1);%filtered indices that can be used to filter-out the RD signals.
	obj.chrFIndex(chrIndex) = whiteRegions;
	%
	chrData = chrSelectedDictionary(chrIndex);
	obj.chrFDictionary(chrIndex) = chrData(whiteRegions);
	
end
%%%

end
%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trueBins] = filterLargeBins(RDmethod, track, threshold)
%%%%%%%%
m = length(track);
trueBins = (track > threshold) & (~isnan(track));

%%% find continous-regions with low scores. 
if(RDmethod == 1 || RDmethod == 2)
	n = 15;%Minimum size of CNVs (in Bins)
	a = track;
	newtrack = arrayfun(@(k) nanmean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	newTrueBins = (newtrack > threshold) & (~isnan(newtrack));
	%%
	tmpTrueBins = repelem(newTrueBins,n);
	trueBins(1:m) = trueBins & tmpTrueBins(1:m);
end
%%%

end
%%%
