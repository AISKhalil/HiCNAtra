function obj = CNVcaller(obj)
% Identifying the copy-number events and CNV tracks of each chromosome.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------- Copy-Number (2N) Computing ----------------%%
% Computing CNReference (the copy-number reference of the genome)

% chrs == 0: all chromosomes
% For human, chrs == 23: chrX
% For mouse, chrs == 19: chrX
chrs = obj.targetChrs;
if (chrs == 0) 
    chromosomes = [1:1:length(obj.chrNames)];
else
    chromosomes = chrs;
end
noChrs = length(chromosomes);
%
copyNumberReference(obj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------- Copy-Number Events -----------------------%%
% Detecting both LCVs and FAs from the RD signal for each chromosome. 

for i = 1:noChrs
	targetChrIndex = chromosomes(i);
	
	%++Segmentation
	[shortSegmentStart, shortSegmentEnd, longSegmentStart, longSegmentEnd] = rdSavitzkySegmentation(obj, targetChrIndex);
	[finalSegmentStart, finalSegmentEnd] = shortSegmentsMerging(obj, targetChrIndex, longSegmentStart, longSegmentEnd);

	%++Detection of CNA regions
	[segmentsInfo, regionsPerSegment] = CNRegionCalling(obj, targetChrIndex, shortSegmentStart, shortSegmentEnd, finalSegmentStart, finalSegmentEnd);

	%++Filtering of CNA regions
	[finalSegmentInfo, impRegionsPerSegment] = CNRegionFiltering(obj, targetChrIndex, segmentsInfo, regionsPerSegment);

	%++Output Updating (individual chromosomes)
	[segmentsData, regionsData] = CNResultsUpdate(obj, targetChrIndex, finalSegmentInfo, impRegionsPerSegment);
	obj.segmentsInfoDic(targetChrIndex) = segmentsData;
	obj.regionsInfoDic(targetChrIndex) = regionsData;
	
end

%++Output Updating (genome)
CNResultsUpdateGenome(obj);
computeEventsTracks(obj);

end
%%%
