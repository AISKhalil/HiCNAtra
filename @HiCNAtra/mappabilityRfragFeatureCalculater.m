function obj = mappabilityRfragFeatureCalculater(obj)
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


%%----------- Mappability feature (rFrag-based) -----------%%
readLength = obj.readLength;
%% For using Anshul's mappability tracks
if(readLength == 100)
	readLength = readLength+1;
end

binSize = obj.binSize;
mappabilityWindow = obj.mappabilityWindow;
mappabilityFolder = obj.mappabilityFolder;
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
rsites = obj.rsites;
largeRfrags = obj.chrLargeFragments;

%%% mappability score per restriction-fragment.
rfragMappScoreDict = containers.Map({1},{[]});
remove(rfragMappScoreDict,1);

for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(i));
    chrLengthBps = chrLengths(chrIndex);
    chrLengthBins = ceil(chrLengthBps/binSize);
    chrRsites = rsites(chrIndex);

    %%%%%%%%%%%%%%%%%%%%%%
    % Anshul_Mappability %

    mappFileName = strcat(mappabilityFolder, chromosome, '.uint8.unique')
    tmp_uMap = fopen(mappFileName,'r');
    uMapdata = fread(tmp_uMap,'*uint8');
    fclose(tmp_uMap);
    mappedBps = (uMapdata > 0 & uMapdata<=readLength);

    %mappability score per restriction-fragment.
    mappabilityWindowLeft = max([chrRsites - mappabilityWindow, chrRsites-1],1);
    mappabilityWindowRight  = min([chrRsites, chrRsites + mappabilityWindow - 1],chrLengthBps);
    [rfragMappScore] = rfragScoring(mappedBps, chrRsites, mappabilityWindowLeft, mappabilityWindowRight, mappabilityWindow, chrLengthBps);
    rfragMappScoreDict(chrIndex) = rfragMappScore;

end

%%% mappability score per bin.
obj.chrMappabilityTracks = RDsignalfromRfrag(rfragMappScoreDict, rsites, chrLengths, binSize, noChrs, chromosomes,largeRfrags);


end
%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rfragScore] = rfragScoring(data, chrRsites, rsitesLeftWindow, rsitesRightWindow, windowSize, chrLengthBps)
%%%%%%%%
chrRfragStart  = [1; chrRsites];
chrRfragEnd    = [chrRsites-1; chrLengthBps]; 
chrRfragLength = chrRfragEnd - chrRfragStart + 1;
twoCompleteWindows = 2*windowSize;
%
noRfragments = length(chrRsites) + 1;
rfragScore = zeros(noRfragments,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First restriction-fragment,
% has only one window.
w1Start  = rsitesLeftWindow(1,1);
w1End    = rsitesLeftWindow(1,2);
w1Length = w1End - w1Start + 1;
w1Data   = data(w1Start:w1End);
rfragScore(1) = sum(w1Data)/w1Length;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last restriction-fragment,
% has only one window.
w1Start  = rsitesRightWindow(end,1);
w1End    = rsitesRightWindow(end,2);
w1Length = w1End - w1Start + 1;
w1Data   = data(w1Start:w1End);
rfragScore(end) = sum(w1Data)/w1Length;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other restriction-fragments,
% each have two restriction fragments
w1Array = rsitesRightWindow(1:end-1,:); 
w2Array = rsitesLeftWindow(2:end,:); 

for i = 2:noRfragments-1
	wIndex = i-1;
	rfragLen = chrRfragLength(i);
	%
	if(rfragLen > twoCompleteWindows)
		%rFrag contains two non-overlapped complete windows
		w1Start  = w1Array(wIndex,1);
		w1End    = w1Array(wIndex,2);
		w2Start  = w2Array(wIndex,1);
		w2End    = w2Array(wIndex,2);
		%
		wData    = [data(w1Start:w1End);data(w2Start:w2End)];
		rfragScore(i) = sum(wData)/twoCompleteWindows; 
	else
		%rFrag contains two overlapped complete windows 
		wData    = data(chrRfragStart(i):chrRfragEnd(i));
		rfragScore(i) = sum(wData)/rfragLen; 		 
	end

end

end
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = RDsignalfromRfrag(rfragCount, rsites, chrLengths, binSize, noChrs, chromosomes, largeRfrags)
%Each Restriction-fragment is counted once

chrDict = containers.Map({1},{[]});
remove(chrDict,1);
%%%%%%%%%%%%%%%%%%%
for i  = 1:1:noChrs 
    chrIndex      = chromosomes(i);
    chrLengthBps  = chrLengths(chrIndex);
    readCountsBps = zeros(chrLengthBps,1);
    %%
    chrRfragCount = rfragCount(chrIndex);
    chrRsites     = rsites(chrIndex);
    chrRfragStart = [1; chrRsites];
    chrRfragEnd   = [chrRsites-1; chrLengthBps]; 
    %% Remove counts for large restriction-fragments.
    chrRfragCount(largeRfrags(chrIndex)) = 0;

    %%
    noReads = length(chrRfragCount);
    for j = 1:noReads
	rfCount = chrRfragCount(j);
	rfStart = chrRfragStart(j);
	rfEnd   = chrRfragEnd(j);
	readCountsBps(rfStart:rfEnd) = rfCount;
    end
    %%
    chrLengthBins = ceil(chrLengthBps/binSize);
    n = binSize;
    a = readCountsBps;
    readCountsBins = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
    %%
    chrDict(chrIndex) = readCountsBins;
end

end
