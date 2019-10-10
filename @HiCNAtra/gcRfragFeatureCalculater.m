function obj = gcRfragFeatureCalculater(obj)
%Computing the gc feature per bin using the gc-window.


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



%%----------- GC feature (window-based) -----------%%
binSize = obj.binSize;
rsites = obj.rsites;
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
%
gcCalculationMethod = obj.gcCalculationMethod;
gcWindow = obj.gcWindow;
gcWindsFolder = obj.gcWindsFolder;
referenceGenomeFolder = obj.referenceGenomeFolder;
largeRfrags = obj.chrLargeFragments;

%%% GC score per restriction-fragment.
rfragGCScoreDict = containers.Map({1},{[]});
remove(rfragGCScoreDict,1);

%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(i));
    chrLengthBps = chrLengths(chrIndex);
    chrRsites = rsites(chrIndex);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % ChrisMiller_GCContents %
    if(gcCalculationMethod == 1)

        gcWindsFileName = strcat(gcWindsFolder, chromosome, '.gc')
        if exist(gcWindsFileName, 'file') == 0
            gunzip(strcat(gcWindsFileName,'.gz'));
        end

        fileID      = fopen(gcWindsFileName,'r');
        GC_Data     = textscan(fileID, '%s');
        GC_Contents = GC_Data{1};
        %%%%
        gcBps = nan(chrLengthBps,1);
        for j = 1:length(GC_Contents)
	    %Each value represent the GC-ratio of 100 bps.
            if(strcmp(GC_Contents(j),'NA')==0)
                segmentStart = (j-1)*100+1;
                segmentStop  = min(j*100,chrLengthBps);
                gcBps(segmentStart:segmentStop)= str2num(cell2mat(GC_Contents(j)));
            end
        end

	%%% GC score per restriction-fragment.
    gcWindsWindowLeft  = max([chrRsites - gcWindow, chrRsites-1],1);
    gcWindsWindowRight = min([chrRsites, chrRsites + gcWindow - 1], chrLengthBps);
	[rfragGCScore] = rfragScoring(gcBps, chrRsites, gcWindsWindowLeft, gcWindsWindowRight, gcWindow, chrLengthBps);
	rfragGCScoreDict(chrIndex) = rfragGCScore;



    %%%%%%%%%%%%%%%%%%%%%%%%%
    % UCSC reference genome %
    elseif(gcCalculationMethod == 2)

        FASTAfileName = strcat(referenceGenomeFolder, chromosome, '.fa')
        fidIn = fopen(FASTAfileName,'r');
        rawData = fread(fidIn,'*char')';
        fclose(fidIn);
        %
        firstBlock = rawData(1:100);
        NL = regexp(firstBlock, '[\n]');
        firstCharIndex = NL(1)+1;% to remove first line
        filteredData = rawData(firstCharIndex:end);
        chrSequence  = strrep(filteredData,sprintf('\n'),'');

	%%% GC score per restriction-fragment.
        gcWindsWindowLeft  = max([chrRsites - gcWindow, chrRsites-1],1);
        gcWindsWindowRight = min([chrRsites, chrRsites + gcWindow - 1], chrLengthBps);       
        [rfragGCScore] = rfragScoringFromSequence(chrSequence, chrRsites, gcWindsWindowLeft, gcWindsWindowRight, gcWindow, chrLengthBps);
        rfragGCScoreDict(chrIndex) = rfragGCScore;

    end
end   

%%% GC-score per bin.
obj.chrGCTracks = RDsignalfromRfrag(rfragGCScoreDict, rsites, chrLengths, binSize, noChrs, chromosomes,largeRfrags);


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
rfragScore(1) = nanmean(w1Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last restriction-fragment,
% has only one window.
w1Start  = rsitesRightWindow(end,1);
w1End    = rsitesRightWindow(end,2);
w1Length = w1End - w1Start + 1;
w1Data   = data(w1Start:w1End);
rfragScore(end) = nanmean(w1Data);

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
		rfragScore(i) = nanmean(wData); 
	else
		%rFrag contains two overlapped complete windows 
		wData    = data(chrRfragStart(i):chrRfragEnd(i));
		rfragScore(i) = nanmean(wData); 		 
	end

end

end
%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rfragScore] = rfragScoringFromSequence(data, chrRsites, rsitesLeftWindow, rsitesRightWindow, windowSize, chrLengthBps)
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
rfragScore(1) = rfragGCRatio(w1Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last restriction-fragment,
% has only one window.
w1Start  = rsitesRightWindow(end,1);
w1End    = rsitesRightWindow(end,2);
w1Length = w1End - w1Start + 1;
w1Data   = data(w1Start:w1End);
rfragScore(end) = rfragGCRatio(w1Data);

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
		rfragScore(i) = rfragGCRatio(wData);  
	else
		%rFrag contains two overlapped complete windows 
		wData    = data(chrRfragStart(i):chrRfragEnd(i));
		rfragScore(i) = rfragGCRatio(wData); 		 
	end

end

end
%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gcRatio] = rfragGCRatio(data)

% find the GC and AT content
gc = (sum(data == 'G' | data == 'C' | data == 'g' | data == 'c'));
at = (sum(data == 'A' | data == 'T' | data == 'a' | data == 't'));
% calculate the ratio of GC to the total known nucleotides
gcRatio = gc/(gc+at);   

end
%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    readCountsBins = arrayfun(@(k) nanmean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
    %%
    chrDict(chrIndex) = readCountsBins;
end

end







