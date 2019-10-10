function [segmentsInfo, regionsPerSegment] = CNRegionCalling(obj, targetChrIndex, orgSegmentStart, orgSegmentEnd, finalSegmentStart, finalSegmentEnd)
%identifying the significant alteration (amplified and deletion regions) using the coverage-based thresholds.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- Parameters ---------------------%%
% Copy-number reference
CNReference = obj.CNReference;
resolution    = obj.resolution;

% Coverage threshold for practically siginificant regions
coverageThAmp = obj.amplificationThreshold;
coverageThDel = obj.deletionThreshold;

% Dictionaries
segmentsInfo = containers.Map({1},{[]});
regionsPerSegment = containers.Map({1},{[]});
remove(segmentsInfo, 1);
remove(regionsPerSegment, 1);
%
regionsPerSegmentInfo0 = containers.Map({1},{[]});
regionsPerSegmentInfo1 = containers.Map({1},{[]});
regionsPerSegmentInfo2 = containers.Map({1},{[]});
remove(regionsPerSegmentInfo0, 1);
remove(regionsPerSegmentInfo1, 1);
remove(regionsPerSegmentInfo2, 1);
%
warning('off','all')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------- Region Classification ----------------%%
orgSignal = obj.chrFNDictionary(targetChrIndex)*2/CNReference;
noFinalSegments = length(finalSegmentStart);
%
for j =1:noFinalSegments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Large Segment %%%%%%%%%%%%%%%%%%%%%%%%%%%
    LSegmentStart = finalSegmentStart(j);
    LSegmentEnd   = finalSegmentEnd(j);
    LSegmentWidth = LSegmentEnd - LSegmentStart +1;
    
    %%%----------Segment Parameters -----------%%%
    sample = orgSignal(LSegmentStart:LSegmentEnd);
    LSegmentRD = median(sample);
    segmentsInfo(j) = [LSegmentStart, LSegmentEnd, LSegmentWidth, LSegmentRD];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% regions/Segment %%%%%%%%%%%%%%%%%%%%%%%%
    sSegStartNo = min(find(orgSegmentStart >= LSegmentStart));
    sSegEndNo   = max(find(orgSegmentEnd <= LSegmentEnd));
    %%
    if(sSegEndNo >= sSegStartNo)
        IBsegmentsStart = orgSegmentStart(sSegStartNo: sSegEndNo);
        IBsegmentsEnd   = orgSegmentEnd(sSegStartNo: sSegEndNo);
        if((IBsegmentsStart(1) - LSegmentStart) > resolution)
            IBsegmentsEnd   = [IBsegmentsStart(1)-1; IBsegmentsEnd];
            IBsegmentsStart = [LSegmentStart; IBsegmentsStart];
        end
        %
        if((LSegmentEnd - IBsegmentsEnd(end)) > resolution)
            IBsegmentsStart = [IBsegmentsStart; IBsegmentsEnd(end)+1];
            IBsegmentsEnd   = [IBsegmentsEnd; LSegmentEnd];
        end
        %
        noIBsegments = length(IBsegmentsStart);

        regionArray = [];%[Start, End, Width, RD, Interval-Category].
        for k = 1:noIBsegments
            regionStart = IBsegmentsStart(k);
            regionEnd   = IBsegmentsEnd(k);
            regionWidth = regionEnd - regionStart +1;
            %
            sample = orgSignal(regionStart:regionEnd);
    	    regionRD = median(sample);
        	if(LSegmentRD < 0.5 && regionRD == 0)
        		regionCategory = 1;
        	else
        	    regionCategory = regionClassify(sample, LSegmentRD, LSegmentWidth, coverageThAmp, coverageThDel);
        	end
            %
            regionParameters = [regionStart, regionEnd, regionWidth, regionRD,regionCategory];
            regionArray = [regionArray; regionParameters];
        end
    else
        regionArray = [];%[Start, End, Width, RD, Interval-Category].
    end
    regionsPerSegmentInfo0(j) = regionArray;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------- Significant Regions Merging ----------%%
%% Level#1 Merging: Same Level
for j = 1:noFinalSegments
    %%%Read
    segmentInfo = segmentsInfo(j);
    %%%Processing (Merging)
    LSegmentRD = segmentInfo(4);
    [mergedRegionArray] = L1RegionMerging(orgSignal, regionsPerSegmentInfo0(j), LSegmentRD);
    %%%Write-Back (After Merging)
    regionsPerSegmentInfo1(j) = mergedRegionArray;
end

regionsPerSegment = regionsPerSegmentInfo1;



%%%
%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------- Sub-routines -----------------------------------------------------------%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------ Category Assignment ---------------------------%%
function [regionCategory] = regionClassify(sample, LSegmentRD, LSegmentWidth, coverageThAmp, coverageThDel)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%% Threshold for top regions %%
	ampThreshold = coverageThAmp;
	delThreshold = coverageThDel;		
	if(abs(LSegmentRD - delThreshold) <= 0.1)
		delThreshold = LSegmentRD - 0.1;
	end

	regionRD    = median(sample);
	ampRegionRD = regionRD;
	delRegionRD = regionRD;


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Category %%%%%%%%%%%%
    cond1 = (regionRD >= LSegmentRD);
	if(cond1 ==1)
	    cond2 = abs(ampRegionRD - LSegmentRD) >= ampThreshold;
	else
	    cond2 = abs(LSegmentRD  - delRegionRD) >= delThreshold;
	end
	%
	if(cond1 == 1)
		%%Big difference (> threshold)
		if(cond2 == 1)
			regionCategory = 5;
		else
			regionCategory = 1;
		end
	else
		%%Big difference (< threshold)
		if(cond2 == 1)
			regionCategory = -5; 
		else
			regionCategory = 1;
		end           
	end
end
%%%
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------------------- Level#1 Merging  -------------------------------%%
function [mergedRegionArray] = L1RegionMerging(signal, regionArray, LSegmentRD)
%%%%%%%%
	[noRegions,noParameters] = size(regionArray);
	if(noRegions > 0)
		regionStart = regionArray(:,1);
		regionEnd   = regionArray(:,2);
		regionRD    = regionArray(:,4);
		regionCategory = regionArray(:,5);
		%%%
		mergedRegionIndex = zeros(noRegions,1);
		regionNumber = 1;
		mergedRegionIndex(1) = regionNumber;
		for j=2:noRegions
				if((regionCategory(j) == regionCategory(j-1)) && ((abs(regionRD(j) - regionRD(j-1)) < 1) || (min(regionRD(j),regionRD(j-1)) > 10)))
					mergedRegionIndex(j) = regionNumber;
				else	
					regionNumber = regionNumber + 1;
					mergedRegionIndex(j) = regionNumber;
				end	
		end
		%%%
		lastRegion = regionNumber;
		mergedRegionArray = zeros(lastRegion,noParameters);
		j = 1;
		while(j<=lastRegion)
			regionIndices = find(mergedRegionIndex == j);
			regionSL = min(regionIndices);
			regionEL  = max(regionIndices);
			%%%%%%%%%
			  LRegionStart = regionStart(regionSL);
			  LRegionEnd   = regionEnd(regionEL);
			  LRegionWidth = LRegionEnd-LRegionStart+1;
			  LRegionCategory = regionCategory(regionSL);
			  %%%%%%%%%
			  LRegionSignal = signal(LRegionStart:LRegionEnd);
			  LRegionRD = median(LRegionSignal);
			  mergedRegionArray(j,:) = [LRegionStart, LRegionEnd, LRegionWidth, LRegionRD, LRegionCategory];	
			%%%%%%%%%
			j = j+1;
		end
	else
		mergedRegionArray = [];
	end
end
%%%
%%%



