function obj = contactMapBiasExtraction(obj)
%%%%%%%%
contactMapEffLenFeature(obj);
contactMapGCFeature(obj);
contactMapMappFeature(obj);
contactMapCNVsFeature(obj);

end
%%%----- Contact-map bias extraction -----%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Sub-routines %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = contactMapEffLenFeature(obj)
%Computing the effective length feature per bin.


%%----------- chromosomes for Analysis -----------%%
chrNames = obj.chrNames;
chrs = obj.targetChrs;
if (chrs == 0) 
    chromosomes = [1:1:length(obj.chrNames)];
else
    chromosomes = chrs;
end
noChrs = length(chromosomes);

%%----------- Effective length feature -----------%%
effectiveLengthWindow = obj.effectiveLengthWindow;
%
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
rsites = obj.rsites;
%
binSize = obj.contactMapBinSize;

for i  = 1:1:noChrs 
	chrIndex = chromosomes(i);
	chromosome = cell2mat(chrNames(i));

	chrLengthBps = chrLengths(chrIndex);
	cut_sites2   = rsites(chrIndex);

	%%effectiveLengthWindow starts from "effectiveLengthWindow before Restriction_sites" to  "effectiveLengthWindow after Restriction_sites"
	effectiveLengthWindowStart = max(cut_sites2 - effectiveLengthWindow, 1);
	effectiveLengthWindowStop  = min(cut_sites2 + effectiveLengthWindow - 1,chrLengthBps);

	effectiveLengthBps  = zeros(chrLengthBps,1);
	for j = 1:length(cut_sites2)
	effectiveLengthBps(effectiveLengthWindowStart(j):effectiveLengthWindowStop(j)) = 1;
	end

	a = effectiveLengthBps;
	n = binSize;
	effectiveLengthBins = arrayfun(@(k) sum(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	%%
	obj.contactMapEffectiveLengthTracks(chrIndex) = effectiveLengthBins;
end    
%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = contactMapGCFeature(obj)
%Computing the gc feature per bin using the gc-window.

%%----------- chromosomes for Analysis ------------%%
chrs = obj.targetChrs;

if (chrs == 0) 
    chromosomes = [1:1:length(obj.chrNames)];
else
    chromosomes = chrs;
end
noChrs = length(chromosomes);


%%----------- gc feature (window-based) -----------%%
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
%
rsites = obj.rsites;
gcWindow = obj.gcWindow;
%
binSize = obj.contactMapBinSize;


%%%-------- checking GC-windows ----------%%%
if(length(obj.GCWindsWindowScores.keys) == 0)
	gcWindowsScores(obj);
end


%%%-------- computing GC-score/bin -------%%%
for i  = 1:1:noChrs 
	chrIndex = chromosomes(i);
	chromosome = cell2mat(chrNames(i));
	chrLengthBps = chrLengths(chrIndex);
	chrRsites = rsites(chrIndex);
	%%
	gcWindsWindowScores = obj.GCWindsWindowScores(chrIndex);
	%%
	chrLengthBins = ceil(chrLengthBps/binSize);
	gcBins = zeros(chrLengthBins,1);
	gcBinsStart = 0:binSize:chrLengthBps;
	gcBinsStop  = [binSize-1:binSize:chrLengthBps,chrLengthBps];
	for j = 1:chrLengthBins
	    rsitesPerBin = find(chrRsites >= gcBinsStart(j) & chrRsites < gcBinsStop(j));
	    rsitesFragments = gcWindsWindowScores(rsitesPerBin,:);
	    rsitesFragments = rsitesFragments(:);
	    gcBins(j) = nanmean(rsitesFragments);
	end
	%%
	obj.contactMapGCTracks(chrIndex) = gcBins;
	%%
end
%%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = contactMapMappFeature(obj)
%Computing the mappability feature per bin using the mappability window.

%%----------- chromosomes for Analysis -----------%%
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


%%----------- Mappability feature (window-based) -----%%
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
%
rsites = obj.rsites;
%
binSize = obj.contactMapBinSize;


%%%---------- checking Mappability-windows ----------%%%
if(length(obj.mappabilityWindowScores.keys) == 0)
	mappabilityWindowsScores(obj);
end


%%%-------- computing GC-score/bin -------%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
	chromosome = cell2mat(chrNames(i));
	chrLengthBps = chrLengths(chrIndex);
	chrLengthBins = ceil(chrLengthBps/binSize);
	chrRsites = rsites(chrIndex);
	%%
	mappabilityWindowScores = obj.mappabilityWindowScores(chrIndex);
	%%
	mappabilityBins = zeros(chrLengthBins,1);
	mappabilityBinsStart = 0:binSize:chrLengthBps;
	mappabilityBinsStop  = [binSize-1:binSize:chrLengthBps,chrLengthBps];

	for j = 1:length(mappabilityBins)
	    rsitesPerBin = find(mappabilityBinsStart(j) <= chrRsites & chrRsites < mappabilityBinsStop(j));
	    rsitesFragments = mappabilityWindowScores(rsitesPerBin,:);
	    if(~isempty(rsitesFragments))
		mappabilityBins(j) = nanmean(rsitesFragments(:));
	    end
	end
	%%%
	obj.contactMapMappabilityTracks(chrIndex) = mappabilityBins;    
end    
%%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = contactMapCNVsFeature(obj)  
% Compute CNV tracks from IBs/focal amplifcations and deletions.

%% ---------------- Data --------------- %%
chrNames = obj.chrNames;
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference  = obj.CNReference;
%
binSize = obj.binSize;
contactMapBinSize = obj.contactMapBinSize;
scalingRatio = contactMapBinSize/binSize;

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



%% ----------- Processing --------------%%
for i=1:1:noChrs
	targetChrIndex = chromosomes(i);
	chrLength = obj.chrLengths(targetChrIndex);
	segmentsCNVector = obj.chrCNVsTracks(targetChrIndex);

	%%%%%%%%%%%%%%%%%%%%%
	if(scalingRatio == 1)
		obj.contactMapCNVsTracks(targetChrIndex) = segmentsCNVector;
	else
		x = scalingCNVsTracks(segmentsCNVector, chrLength, binSize, contactMapBinSize);
		obj.contactMapCNVsTracks(targetChrIndex) = x;
	end
end
%%%


end
%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nCNVTrack] = scalingCNVsTracks(CNVTrack, chrLengthBps, RdBinSize, contactMapBinSize)
%%%%%%%%

%%%%%---------------- 1Kb Bin-size --------------%%%%%
targetChr1KbLength = ceil(chrLengthBps/1000);
scalingRatio = RdBinSize/1000;
CNVTrack1Kb = repelem(CNVTrack, scalingRatio);
CNVTrack1Kb = CNVTrack1Kb(1:targetChr1KbLength);

%%%%%---------- Contact-map Bin-size ------------%%%%%
a = CNVTrack1Kb;
n = contactMapBinSize/1000;
nCNVTrack = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';% the averaged vector


%%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nCNVTrack] = scalingCNVsTracksToLess1Kb(segmentsStart, segmentsCN, chrLength, binSize, contactMapBinSize)
%%%%%%%%

[segmentsStart, I] = sort(segmentsStart);
segmentsCN = segmentsCN(I);

%% segment-starts @bp
segmentsStart = ((segmentsStart-1)*binSize) + 1;
segmentsStart = [segmentsStart; chrLength];
noSegments = length(segmentsStart);

%% bp-level
a = zeros(chrLength,1);
for i = 1:noSegments - 1
	a(segmentsStart(i):segmentsStart(i+1)-1) = segmentsCN;
end

%% binning
n = contactMapBinSize;
nCNVTrack = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';% the averaged vector


%%%
end
