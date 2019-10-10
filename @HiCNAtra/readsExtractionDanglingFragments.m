function obj = readsExtractionDanglingFragments(obj)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chrLengthsVector = cell2mat(values(chrLengths));
maxChrmsLength = max(chrLengthsVector);
genomeLengthBps = sum(chrLengthsVector);
clear chrLengthsVector;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readsStatTot = [0,0,0,0,0,0];%Double-side, Informative, Self-circle, Dangling-end, extraDangling-end, single-side
totalReads = 0;
%%%
chrms1Dangling = [];
chrms2Dangling = [];
cuts1Dangling  = [];
cuts2Dangling  = [];
%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- Reads Extraction & Duplicated Removal & Computing RD-signal--------------%
disp('Computing the RD-signal using the DE/extra DE fragments approach ........')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------%
% Dangling-end reads (fragments)
% -----------------------------%
disp('Extracting DE/Extra DE Reads....')

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
        [~, ~, danglingEndMask, extraDanglingEndMask, ~, ~, readsStat] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);
        readsStatTot = readsStatTot + readsStat;

        %% Dangling Reads
        danglingMask = danglingEndMask | extraDanglingEndMask;
        chrms1Dangling = [chrms1Dangling; uint8(chrms1Data(danglingMask))];
        chrms2Dangling = [chrms2Dangling; uint8(chrms2Data(danglingMask))];
        cuts1Dangling  = [cuts1Dangling;  uint32(cuts1Data(danglingMask))];
        cuts2Dangling  = [cuts2Dangling;  uint32(cuts2Data(danglingMask))];


    end
end
clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data;
clear infoMask nonInfoMask danglingMask selfCircleMask danglingEndMask extraDanglingEndMask singleSide1Mask singleSide2Mask; 


disp('Removing duplicated DE/Extra DE reads')
[uniqueMaskDangling] = removeDuplicatedMemoryEfficient(1,chrms1Dangling, chrms2Dangling, cuts1Dangling, cuts2Dangling, noChrs, chrLengths,memoryFootPrint);
chrms1Dangling = chrms1Dangling(uniqueMaskDangling);
chrms2Dangling = chrms2Dangling(uniqueMaskDangling);
cuts1Dangling = cuts1Dangling(uniqueMaskDangling);
cuts2Dangling = cuts2Dangling(uniqueMaskDangling);
clear uniqueMaskDangling;

disp('Computing DE/Extra DE RD-signal....')
[obj.chrDictionary, totalDanglingBpsCounts]= RDsignalBpsCountDangling(chrms1Dangling, chrms2Dangling, cuts1Dangling, cuts2Dangling, binSize, chrLengths, noChrs, chromosomes);
clear chrms1Dangling chrms2Dangling cuts1Dangling cuts2Dangling;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Data Statistics %%%%%%%%%%%%%%%%%%%%%

% readStatistics = [total-reads, double-side, informative, self-circle, dangling-end, extraDangling-end, single-side];
readStatistics = [totalReads, readsStatTot];
obj.readStatistics = readStatistics;

% total coverage
obj.dataCoverage = totalDanglingBpsCounts/ genomeLengthBps; 








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
function [chrDict,totalDanglingBpsCounts] = RDsignalBpsCountDangling(chrms1, chrms2, cuts1, cuts2, binSize, chrLengths, noChrs, chromosomes)
%
% Counting Bps after forming dangling-fragments.

chrDict = containers.Map({1},{[]});
remove(chrDict,1);
totalDanglingBpsCounts = 0;
%%%%%%%%%%%%%%%%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chrLengthBps = chrLengths(chrIndex);
    readCountsBps = zeros(chrLengthBps,1);
    %%
    chrIndexHiCLib = chrIndex - 1;
    chrmsIndex = (chrms1 == chrIndexHiCLib);% & (chrms2 == chrIndexHiCLib);
    noReads = sum(chrmsIndex);
    %%
    sCuts1Bps  = cuts1(chrmsIndex);
    sCuts2Bps  = cuts2(chrmsIndex);

    %%
    readFragmentStart = min(sCuts1Bps,sCuts2Bps);
    readFragmentStop   = max(sCuts1Bps,sCuts2Bps);
    for j = 1:noReads
        readCountsBps(readFragmentStart(j):readFragmentStop(j)) = readCountsBps(readFragmentStart(j):readFragmentStop(j)) + 1;
    end
    %%
    %% Binning the readCountsBps.
    %% Coordinates starts from 0 bp: bin#1 has bps from 0 to 999
    a = [0;readCountsBps(1:end-1)];
    n = binSize;
    readCountsBins = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
    %%
    totalDanglingBpsCounts = totalDanglingBpsCounts + sum(readCountsBps);
    chrDict(chrIndex) = readCountsBins;
end

end



