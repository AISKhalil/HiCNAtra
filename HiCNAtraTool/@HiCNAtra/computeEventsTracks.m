function obj = computeEventsTracks(obj)  
% Compute CNV tracks from IBs/focal amplifcations and deletions.

%% ---------------- Data --------------- %%
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference  = obj.CNReference;
%
binSize = obj.binSize;


%%%%%%%%%%%%%%%%%%%%%%
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
for i=1:noChrs
	targetChrIndex = chromosomes(i);
	chrLength = obj.chrLengths(targetChrIndex);

	segmentsCNVector = zeros(length(obj.chrDictionary(targetChrIndex)),1);%length of binned data
	%%------ Segment Data  ------------%%
	segmentsArray  = segmentsInfoDic(targetChrIndex);
	[noSegments,~] = size(segmentsArray);
	segmentsStart  = segmentsArray(:,2);
	segmentsEnd    = segmentsArray(:,3);
	segmentsCN     = segmentsArray(:,5);	
	for j=1:noSegments
		segmentsCNVector(segmentsStart(j):segmentsEnd(j)) = segmentsCN(j); 
	end

	%%------ Regions Data  ------------%%	
	regionsArray = regionsInfoDic(targetChrIndex);
	[noRegions,~] = size(regionsArray);	
	if(noRegions > 0)	
		regionsStart    = regionsArray(:,2);
		regionsEnd      = regionsArray(:,3);
		regionsCN       = regionsArray(:,5);
		regionsCategory = regionsArray(:,6);
		for j=1:noRegions
				segmentsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	obj.chrCNVsTracks(targetChrIndex) = segmentsCNVector;
end
%%%


end
%%%