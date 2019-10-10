function obj = readsExtractionExactCutPosition(obj)
%Extracting the info & non-info RD signal.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------- chromosomes for Analysis ------------%%
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------- Input-data -------------------%%
HDF5Files = obj.HDF5Files;
maximumMoleculeLength = obj.maximumMoleculeLength;
memoryFootPrint = obj.memoryFootPrint;
binSize = obj.binSize;
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
rsites = obj.rsites;
readLength = obj.readLength;
readTypes = obj.readTypes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chrLengthsVector = cell2mat(values(chrLengths));
maxChrmsLength = max(chrLengthsVector);
genomeLengthBps = sum(chrLengthsVector);
clear chrLengthsVector;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readsStatTot = [0,0,0,0,0,0];%Double-side, Informative, Self-circle, Dangling-end, extraDangling-end, single-side
totalReads = 0;
%%%
chrms1Info = [];
chrms2Info = [];
cuts1Info  = [];
cuts2Info  = [];
%%%
chrms1NonInfo = [];
chrms2NonInfo = [];
cuts1NonInfo  = [];
cuts2NonInfo  = [];
%%%
chrmsSingleSided = [];
cutsSingleSided = [];
%%%
chrInfoDictionary = containers.Map({1},{[]});
remove(chrInfoDictionary,1);
chrNonInfoDictionary = containers.Map({1},{[]});
remove(chrNonInfoDictionary, 1);
chrSSDictionary = containers.Map({1},{[]});
remove(chrSSDictionary, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- Reads Extraction & Duplicated Removal & Computing RD-signal--------------%
disp('Computing the RD-signal using the exact-cut positions approach ........')


%%%%%%%%%%%%%%%%
% -------------%
% (1) Info reads
% -------------%
if(readTypes == 1 || readTypes == 2)
    disp('Extracting informative Reads....')

    %-- Reads Classification --%
    noFiles = length(HDF5Files);

    for fileNo = 1:noFiles

        %----- File Details ----%
        hiclibDict = cell2mat(HDF5Files(fileNo))
        inputInfo = h5info(hiclibDict);
        inputDimensions = inputInfo.Datasets.Dataspace;
        noReads = inputDimensions.Size;    
        totalReads = totalReads + noReads;
        
        for readOrder = 1:memoryFootPrint:noReads
            startRead = readOrder;
            stopRead  = min(startRead+memoryFootPrint-1,noReads);
            countReads = stopRead - startRead + 1;

            %%
            chrms1Data = h5read(hiclibDict,'/chrms1',startRead,countReads);
            chrms2Data = h5read(hiclibDict,'/chrms2',startRead,countReads);
            cuts1Data  = h5read(hiclibDict,'/cuts1',startRead,countReads);
            cuts2Data  = h5read(hiclibDict,'/cuts2',startRead,countReads);
            rfragIdxs1Data = h5read(hiclibDict,'/rfragIdxs1',startRead,countReads);
            rfragIdxs2Data = h5read(hiclibDict,'/rfragIdxs2',startRead,countReads);
            strands1Data  = h5read(hiclibDict,'/strands1',startRead,countReads);
            strands2Data  = h5read(hiclibDict,'/strands2',startRead,countReads);

            %%
            [infoMask, ~, ~, ~, ~, ~, readsStat] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);
            readsStatTot = readsStatTot + readsStat;

            %% Informative Reads
            chrms1Info = [chrms1Info; uint8(chrms1Data(infoMask))];
            chrms2Info = [chrms2Info; uint8(chrms2Data(infoMask))];
            cuts1Info  = [cuts1Info ; uint32(cuts1Data(infoMask))];
            cuts2Info  = [cuts2Info ; uint32(cuts2Data(infoMask))];

        end
    end
    clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data;
    clear infoMask; 


    disp('Removing Duplicated Info Reads....')
    [uniqueMaskInfo] = removeDuplicatedMemoryEfficient(1,chrms1Info, chrms2Info, cuts1Info, cuts2Info, noChrs, chrLengths, memoryFootPrint);
    chrms1Info = chrms1Info(uniqueMaskInfo);
    chrms2Info = chrms2Info(uniqueMaskInfo);
    cuts1Info = cuts1Info(uniqueMaskInfo);
    cuts2Info = cuts2Info(uniqueMaskInfo);
    NDInfoReads = length(uniqueMaskInfo);
    clear uniqueMaskInfo;

    disp('Computing Info RD-signal....')
    chrInfoDictionary = RDsignalReadCountInfo(chrms1Info, chrms2Info, cuts1Info, cuts2Info, binSize, chrLengths, noChrs, chromosomes);
    clear chrms1Info chrms2Info cuts1Info cuts2Info;
end






%%%%%%%%%%%%%%%%%%%%
% ---------------- %
% (2) Non-Info reads
% -----------------%
if(readTypes == 1 || readTypes == 3)
    disp('Extracting Non-informative Reads....')

    %-- Reads Classification --%
    noFiles = length(HDF5Files);

    for fileNo = 1:noFiles

        %----- File Details ----%
        hiclibDict = cell2mat(HDF5Files(fileNo))
        inputInfo = h5info(hiclibDict);
        inputDimensions = inputInfo.Datasets.Dataspace;
        noReads = inputDimensions.Size;    
        
        for readOrder = 1:memoryFootPrint:noReads
            startRead = readOrder;
            stopRead  = min(startRead+memoryFootPrint-1,noReads);
            countReads = stopRead - startRead + 1;

            %%
            chrms1Data = h5read(hiclibDict,'/chrms1',startRead,countReads);
            chrms2Data = h5read(hiclibDict,'/chrms2',startRead,countReads);
            cuts1Data  = h5read(hiclibDict,'/cuts1',startRead,countReads);
            cuts2Data  = h5read(hiclibDict,'/cuts2',startRead,countReads);
            rfragIdxs1Data = h5read(hiclibDict,'/rfragIdxs1',startRead,countReads);
            rfragIdxs2Data = h5read(hiclibDict,'/rfragIdxs2',startRead,countReads);
            strands1Data  = h5read(hiclibDict,'/strands1',startRead,countReads);
            strands2Data  = h5read(hiclibDict,'/strands2',startRead,countReads);

            %%
            [~, selfCircleMask, danglingEndMask, extraDanglingEndMask, ~, ~, ~] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);

            %% Non-informative Reads
            nonInfoMask = selfCircleMask | danglingEndMask | extraDanglingEndMask;
            chrms1NonInfo = [chrms1NonInfo;uint8(chrms1Data(nonInfoMask))];
            chrms2NonInfo = [chrms2NonInfo;uint8(chrms2Data(nonInfoMask))];
            cuts1NonInfo  = [cuts1NonInfo;uint32(cuts1Data(nonInfoMask))];
            cuts2NonInfo  = [cuts2NonInfo;uint32(cuts2Data(nonInfoMask))];

        end
    end
    clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data;
    clear nonInfoMask selfCircleMask danglingEndMask extraDanglingEndMask; 


    disp('Removing Duplicated Non-Info Reads....')
    [uniqueMaskNonInfo] = removeDuplicatedMemoryEfficient(1,chrms1NonInfo, chrms2NonInfo, cuts1NonInfo, cuts2NonInfo, noChrs, chrLengths, memoryFootPrint);
    chrms1NonInfo = chrms1NonInfo(uniqueMaskNonInfo);
    chrms2NonInfo = chrms2NonInfo(uniqueMaskNonInfo);
    cuts1NonInfo = cuts1NonInfo(uniqueMaskNonInfo);
    cuts2NonInfo = cuts2NonInfo(uniqueMaskNonInfo);
    NDNonInfoReads = length(uniqueMaskNonInfo);
    clear uniqueMaskNonInfo;

    disp('Computing Non-Info RD-signal....')
    chrNonInfoDictionary = RDsignalReadCountNonInfo(chrms1NonInfo, chrms2NonInfo, cuts1NonInfo, cuts2NonInfo, binSize, chrLengths, noChrs, chromosomes);
    clear chrms1NonInfo chrms2NonInfo cuts1NonInfo cuts2NonInfo;
end








%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------%
% (3) Single-sided reads
% ---------------------%
if(readTypes == 1 || readTypes == 4)
    disp('Extracting Single-sided Reads....')

    %-- Reads Classification --%
    noFiles = length(HDF5Files);

    for fileNo = 1:noFiles

        %----- File Details ----%
        hiclibDict = cell2mat(HDF5Files(fileNo))
        inputInfo = h5info(hiclibDict);
        inputDimensions = inputInfo.Datasets.Dataspace;
        noReads = inputDimensions.Size;    
        
        for readOrder = 1:memoryFootPrint:noReads
            startRead = readOrder;
            stopRead  = min(startRead+memoryFootPrint-1,noReads);
            countReads = stopRead - startRead + 1;

            %%
            chrms1Data = h5read(hiclibDict,'/chrms1',startRead,countReads);
            chrms2Data = h5read(hiclibDict,'/chrms2',startRead,countReads);
            cuts1Data  = h5read(hiclibDict,'/cuts1',startRead,countReads);
            cuts2Data  = h5read(hiclibDict,'/cuts2',startRead,countReads);
            rfragIdxs1Data = h5read(hiclibDict,'/rfragIdxs1',startRead,countReads);
            rfragIdxs2Data = h5read(hiclibDict,'/rfragIdxs2',startRead,countReads);
            strands1Data  = h5read(hiclibDict,'/strands1',startRead,countReads);
            strands2Data  = h5read(hiclibDict,'/strands2',startRead,countReads);

            %%
            [~, ~, ~, ~, singleSide1Mask, singleSide2Mask, ~] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);

            %% Single-sided Reads
            chrmsSingleSided = [chrmsSingleSided; uint8(chrms1Data(singleSide1Mask)); uint8(chrms2Data(singleSide2Mask))];
            cutsSingleSided  = [cutsSingleSided;  uint32(cuts1Data(singleSide1Mask));  uint32(cuts2Data(singleSide2Mask))];

        end
    end
    clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data;
    clear singleSide1Mask singleSide2Mask; 

    disp('Removing Duplicated Single-sided Reads....')
    [uniqueMaskSingleSided] = removeDuplicatedMemoryEfficient(0,chrmsSingleSided, [], cutsSingleSided, [], noChrs, chrLengths, memoryFootPrint); 
    chrmsSingleSided = chrmsSingleSided(uniqueMaskSingleSided);
    cutsSingleSided  = cutsSingleSided(uniqueMaskSingleSided);
    NDSingleSidedReads = length(uniqueMaskSingleSided);


    clear uniqueMaskSingleSided;
    disp('Computing Single-sided RD-signal....')
    chrSSDictionary = RDsignalReadCountSS(chrmsSingleSided, cutsSingleSided, binSize, chrLengths, noChrs, chromosomes);
    clear chrmsSingleSided cutsSingleSided;
end










%%%%%%%%%%%%%%%
% ------------%
% All reads
% ------------%
disp('Computing exact-cut position RD-signal....')

for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    %
    if(readTypes == 1)
        obj.chrDictionary(chrIndex) = chrInfoDictionary(chrIndex) + chrNonInfoDictionary(chrIndex) + chrSSDictionary(chrIndex);
    elseif(readTypes == 2)
        obj.chrDictionary(chrIndex) = chrInfoDictionary(chrIndex);
    elseif(readTypes == 3)
        obj.chrDictionary(chrIndex) = chrNonInfoDictionary(chrIndex);        
    elseif(readTypes == 4)
        obj.chrDictionary(chrIndex) = chrSSDictionary(chrIndex);        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------- Coverage Calculation ---------------%%

% readStatistics = [total-reads, double-side, informative, self-circle, dangling-end, extraDangling-end, single-side];
readStatistics = [totalReads, readsStatTot];
obj.readStatistics = readStatistics;

% Number-of-sides
if(readTypes == 1)
    obj.NDupReadStatistics = [NDInfoReads, NDNonInfoReads, NDSingleSidedReads];
    obj.nearRSitesSides = [];
    totalSides = obj.NDupReadStatistics(1)*2 + obj.NDupReadStatistics(2)*2 + obj.NDupReadStatistics(3);
elseif(readTypes == 2)
    obj.NDupReadStatistics = [NDInfoReads];
    obj.nearRSitesSides = [];
    totalSides = obj.NDupReadStatistics(1)*2;
elseif(readTypes == 3)
    obj.NDupReadStatistics = [NDNonInfoReads];
    obj.nearRSitesSides = [];
    totalSides = obj.NDupReadStatistics(1)*2;
elseif(readTypes == 4)
    obj.NDupReadStatistics = [NDSingleSidedReads];
    obj.nearRSitesSides = [];
    totalSides = obj.NDupReadStatistics(1);
end


% total coverage
obj.dataCoverage = totalSides*readLength/genomeLengthBps; 








end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------- End "readsExtraction" ---------------------------------------------------%%




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis Sub-routines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------------------------------------------------------------------%%
%% Sub-routines for extracting the different types of reads.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [infoMask, selfCircleMask, danglingEndMask, extraDanglingEndMask, singleSide1Mask, singleSide2Mask, readsStat] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength)
%%hiclib.fragments.parseInputData():

DSMask = (chrms1Data >=0) & (chrms2Data >= 0);
singleSide1Mask = (chrms1Data >=0) & (chrms2Data <0);
singleSide2Mask = (chrms2Data >=0) & (chrms1Data <0);

%%
sameChrMask = (chrms1Data == chrms2Data);
sameFragmentMask = (rfragIdxs1Data == rfragIdxs2Data);
cutDifsMask = (cuts2Data > cuts1Data);
sameStrandMask = (strands1Data == strands2Data);
cutDiffEqualStrand2 = (cutDifsMask == strands2Data);

%% Converting strands to int64 for calculating the distance
strands1Data = int64(strands1Data);
strands2Data = int64(strands2Data);
dist = (-cuts1Data .* (2.*strands1Data - 1)) - (cuts2Data.*(2.*strands2Data - 1));
extraDEnotInfoMask = (chrms1Data == chrms2Data) & (strands1Data ~= strands2Data) & (dist >= 0) & (dist <= maximumMoleculeLength);

%%
selfCircleMask = DSMask & sameChrMask & sameFragmentMask & (~sameStrandMask) & cutDiffEqualStrand2;
danglingEndMask = DSMask & sameChrMask & sameFragmentMask & (~sameStrandMask) & (~cutDiffEqualStrand2);
extraDanglingEndMask = DSMask & sameChrMask & (~sameFragmentMask) & extraDEnotInfoMask;
infoMask = DSMask & ((sameChrMask &(~sameFragmentMask) & (~extraDEnotInfoMask)) | (~sameChrMask));

readsStat = [sum(DSMask), sum(infoMask), sum(selfCircleMask), sum(danglingEndMask), sum(extraDanglingEndMask), sum(singleSide1Mask) + sum(singleSide2Mask)];

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uniqueIndices] = removeDuplicated(doubleSidedFlag, chrms1Data, chrms2Data, cuts1Data, cuts2Data, maxChrmsLength)
%%fragments.filterDuplicates():

if(doubleSidedFlag == 1)
    %% Converting chrms to int64 for calculating the distance
    chrms1Data = int64(chrms1Data);
    chrms2Data = int64(chrms2Data);
    cuts1Data = int64(cuts1Data);
    cuts2Data = int64(cuts2Data);

    fragIDmult = maxChrmsLength + 1000;
    dups1 = chrms1Data*fragIDmult + cuts1Data;
    clear chrms1Data cuts1Data;
    dups2 = chrms2Data*fragIDmult + cuts2Data;
    clear chrms2Data cuts2Data;

    dups  = [dups1,dups2];
    dups  = sort(dups,2);        
    clear dups1 dups2;
    [~, uniqueIndices, ~] = unique(dups, 'rows');
else
    %% Converting chrms to int64 for calculating the distance
    chrms1Data = int64(chrms1Data);
    cuts1Data = int64(cuts1Data);

    fragIDmult = maxChrmsLength + 1000;

    dups = chrms1Data*fragIDmult + cuts1Data;
    clear chrms1Data cuts1Data;
    [~, uniqueIndices, ~] = unique(dups);

end 
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uniqueIndices] = removeDuplicatedMemoryEfficient(doubleSidedFlag, chrms1Data, chrms2Data, cuts1Data, cuts2Data, noChrs, chrLengths, memoryFootPrint)
%%fragments.filterDuplicates():

if(doubleSidedFlag == 1)
    %% Converting chrms to int32 for calculating the distance
    dups1 = cutAtGenomeCoordinates(chrms1Data, cuts1Data, noChrs, chrLengths);
    clear chrms1Data cuts1Data;
    dups2 = cutAtGenomeCoordinates(chrms2Data, cuts2Data, noChrs, chrLengths);
    clear chrms2Data cuts2Data;
	%
    dups  = [dups1,dups2];
    dups  = sort(dups,2);    
    clear dups1 dups2;
	%
    [uniqueIndices] = pairedUnique(dups, memoryFootPrint);
	clear dups;
else

    dups = cutAtGenomeCoordinates(chrms1Data, cuts1Data, noChrs, chrLengths);
    clear chrms1Data cuts1Data;
	%
    [uniqueIndices] = singleUnique(dups, memoryFootPrint);
	clear dups;
end
end

%%%%%%%%
function [genomeCuts] = cutAtGenomeCoordinates(chrms, cuts, noChrs, chrLengths)
    %%%chrms: int8, cuts: int64
    genomeCuts = zeros(length(chrms),1,'uint32');
    cuts = uint32(cuts);
    chrStartBps = uint32(0);
    for i = 1:noChrs
        chrLength = chrLengths(i);
        chrIndexHiCLib = i - 1;
        genomeCuts(chrms == chrIndexHiCLib) = chrStartBps + cuts(chrms == chrIndexHiCLib);
        chrStartBps = chrStartBps + uint32(chrLength);
    end
end

%%%%%%%% 
function [uniqueIndices] = pairedUnique(dups, memoryFootPrint)
	dups1 = dups(:,1);
	[~,sortedIndex]  = sort(dups1);
	dupsSorted = dups(sortedIndex,:);
	clear dups dups1;
	%
	[noReads, ~] = size(dupsSorted);
	newUniqueIndex = [];
	for i = 1:memoryFootPrint:noReads
		subsetEnd = min(i+memoryFootPrint-1,noReads);
		subset    = dupsSorted(i:subsetEnd,:);
		[~, setUniqueIndices, ~] = unique(subset, 'rows');
		newUniqueIndex = [newUniqueIndex; i-1+setUniqueIndices];
	end
	%
	uniqueIndices = sortedIndex(newUniqueIndex);
	clear subset dupsSorted sortedIndex newUniqueIndex;
end

%%%%%%%%
function [uniqueIndices] = singleUnique(dups, memoryFootPrint)
	[~,sortedIndex]  = sort(dups);
	dupsSorted = dups(sortedIndex);
	clear dups;
	%
	[noReads, ~] = size(dupsSorted);
	newUniqueIndex = [];
	for i = 1:memoryFootPrint:noReads
		subsetEnd = min(i+memoryFootPrint-1,noReads);
		subset    = dupsSorted(i:subsetEnd);
		[~, setUniqueIndices, ~] = unique(subset);
		newUniqueIndex = [newUniqueIndex; i-1+setUniqueIndices];
	end
	%
	uniqueIndices = sortedIndex(newUniqueIndex);
	clear subset dupsSorted sortedIndex newUniqueIndex;
end   






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = RDsignalReadCountInfo(chrms1, chrms2, cuts1, cuts2, binSize, chrLengths, noChrs, chromosomes)
%
% Each side is counted separately. So, if the two-sides are in the same bin, readCount = readCount + 2.

chrDict = containers.Map({1},{[]});
remove(chrDict,1);
%%%%%%%%%%%%%%%%%%
% Coordinates starts from 0 bp: bin#1 has bps from 0 to 999
cuts1Bins = floor(cuts1/binSize) + 1;
cuts2Bins = floor(cuts2/binSize) + 1;
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chrLengthBps = chrLengths(chrIndex);
    chrLengthBins = ceil(chrLengthBps/binSize);
    readCountsBins = zeros(chrLengthBins,1);
    %%
    chrIndexHiCLib = chrIndex - 1;
    chrms1Index = (chrms1 == chrIndexHiCLib);
    chrms2Index = (chrms2 == chrIndexHiCLib);
    cutsBins = [cuts1Bins(chrms1Index); cuts2Bins(chrms2Index)];
    cutsBins = min(cutsBins, chrLengthBins);
    noReads = length(cutsBins);
    %%
    for j = 1:noReads
        cutBin = cutsBins(j);
        readCountsBins(cutBin) = readCountsBins(cutBin) + 1;
    end
    %%
    chrDict(chrIndex) = readCountsBins;
end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = RDsignalReadCountNonInfo(chrms1, chrms2, cuts1, cuts2, binSize, chrLengths, noChrs, chromosomes)
%
% The two-sides are counted together. So, if the two-sides are in the same bin, readCount = readCount + 1. Else, count of each side is incremented 
% by 1: readCount1 = readCount1 + 1 & readCount2 = readCount2 + 1

chrDict = containers.Map({1},{[]});
remove(chrDict,1);
%%%%%%%%%%%%%%%%%%
% Coordinates starts from 0 bp: bin#1 has bps from 0 to 999
cuts1Bins = floor(cuts1/binSize) + 1;
cuts2Bins = floor(cuts2/binSize) + 1;

for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chrLengthBps = chrLengths(chrIndex);
    chrLengthBins = ceil(chrLengthBps/binSize);
    readCountsBins = zeros(chrLengthBins,1);
    %%
    chrIndexHiCLib = chrIndex - 1;
    chrmsIndex = (chrms1 ==chrIndexHiCLib);
    noReads = sum(chrmsIndex);
    %%
    sCuts1Bins  = cuts1Bins(chrmsIndex);
    sCuts2Bins  = cuts2Bins(chrmsIndex);
    sCuts1Bins = min(sCuts1Bins, chrLengthBins);
    sCuts2Bins = min(sCuts2Bins, chrLengthBins);    
    %% conditions %%
    twoSidesSameBin = (sCuts1Bins == sCuts2Bins);
    twoSidesDiffBin = (sCuts1Bins ~= sCuts2Bins);
    %%
    for j = 1:noReads
        if(twoSidesSameBin(j))%Counted once
            readCountsBins(sCuts1Bins(j)) = readCountsBins(sCuts1Bins(j)) + 1;
        elseif(twoSidesDiffBin(j))%Counted twice
            readCountsBins(sCuts1Bins(j)) = readCountsBins(sCuts1Bins(j)) + 1;
            readCountsBins(sCuts2Bins(j)) = readCountsBins(sCuts2Bins(j)) + 1;
        end
    end
    %%
    chrDict(chrIndex) = readCountsBins;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = RDsignalReadCountSS(chrms1, cuts1, binSize, chrLengths, noChrs, chromosomes)
%Binning the Single-sided reads

chrDict = containers.Map({1},{[]});
remove(chrDict,1);
%%%%%%%%%%%%%%%%%%
% Coordinates starts from 0 bp: bin#1 has bps from 0 to 999
cuts1Bins = floor(cuts1/binSize) + 1;

for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chrLengthBps = chrLengths(chrIndex);
    chrLengthBins = ceil(chrLengthBps/binSize);
    readCountsBins = zeros(chrLengthBins,1);
    %%
    chrIndexHiCLib = chrIndex - 1;
    chrms1Index = (chrms1 ==chrIndexHiCLib);
    cutsBins = cuts1Bins(chrms1Index);
    cutsBins = min(cutsBins, chrLengthBins);
    noReads = length(cutsBins);
    %%
    for j = 1:noReads
        cutBin = cutsBins(j);
        readCountsBins(cutBin) = readCountsBins(cutBin) + 1;
    end
    %%
    chrDict(chrIndex) = readCountsBins;
end

end
