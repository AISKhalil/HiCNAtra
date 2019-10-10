function obj = readsExtractionRestrictionFragments(obj)
%Extracting the info & non-info RD signal by counting the restriction fragments of the cuts assuming that each cut means a restriction fragment.


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
readTypes  = obj.readTypes;
cutToRsiteTh = obj.cutToRsiteThreshold;
largeRfrags = obj.chrLargeFragments;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chrLengthsVector = cell2mat(values(chrLengths));
maxChrmsLength = max(chrLengthsVector);
genomeLengthBps = sum(chrLengthsVector);
clear chrLengthsVector;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readsStatTot = [0,0,0,0,0,0];%Double-side, Informative, Self-circle, Dangling-end, ExtraDangling-end, Single-side
totalReads = 0;
%%%
chrms1Info = [];
chrms2Info = [];
cuts1Info  = [];
cuts2Info  = [];
rfrag1Info = [];
rfrag2Info = [];
rsites1Info = [];
rsites2Info = [];
%%%
chrms1NonInfo = [];
chrms2NonInfo = [];
cuts1NonInfo  = [];
cuts2NonInfo  = [];
rfrag1NonInfo = [];
rfrag2NonInfo = [];
rsites1NonInfo = [];
rsites2NonInfo = [];
%%%
chrmsSingleSided  = [];
cutsSingleSided   = [];
rfragsSingleSided = [];
rsitesSingleSided = [];
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- Reads Extraction & Duplicated Removal & Computing RD-signal--------------%
disp('Computing the RD-signal using the restriction-fragments counting approach ........')


%%%%%%%%%%%%%%%%
% -------------%
% (1) Info reads
% -------------%
if(readTypes == 1 || readTypes == 2)
    disp('Extracting informative Reads....')
    %%%

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
            rsites1Data  = h5read(hiclibDict,'/rsites1',startRead,countReads);
            rsites2Data  = h5read(hiclibDict,'/rsites2',startRead,countReads);        

            %%
            [infoMask, ~, ~, ~, ~, ~, readsStat] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);
            readsStatTot = readsStatTot + readsStat;

            %% Informative Reads (before distance-filtering)
            chrms1Info = [chrms1Info; uint8(chrms1Data(infoMask))];
            chrms2Info = [chrms2Info; uint8(chrms2Data(infoMask))];
            cuts1Info  = [cuts1Info ; uint32(cuts1Data(infoMask))];
            cuts2Info  = [cuts2Info ; uint32(cuts2Data(infoMask))];
            rfrag1Info = [rfrag1Info; uint32(rfragIdxs1Data(infoMask))];
            rfrag2Info = [rfrag2Info; uint32(rfragIdxs2Data(infoMask))];
            rsites1Info  = [rsites1Info ; uint32(rsites1Data(infoMask))];
            rsites2Info  = [rsites2Info ; uint32(rsites2Data(infoMask))];

            %% 

        end
    end

    clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data rsites1Data rsites2Data;
    clear infoMask;  

    disp('Removing Duplicated Info Reads....')
    [uniqueMaskInfo] = removeDuplicatedMemoryEfficient(1,chrms1Info, chrms2Info, cuts1Info, cuts2Info, noChrs, chrLengths, memoryFootPrint);
    chrms1Info = [chrms1Info(uniqueMaskInfo)];
    chrms2Info = [chrms2Info(uniqueMaskInfo)];
    cuts1Info  = [cuts1Info(uniqueMaskInfo)];
    cuts2Info  = [cuts2Info(uniqueMaskInfo)];
    rfrag1Info = [rfrag1Info(uniqueMaskInfo)];
    rfrag2Info = [rfrag2Info(uniqueMaskInfo)];
    rsites1Info  = [rsites1Info(uniqueMaskInfo)];
    rsites2Info  = [rsites2Info(uniqueMaskInfo)];
    NDInfoReads = length(uniqueMaskInfo);
    clear uniqueMaskInfo;

    disp('Removing Info Reads with large cut-rsite distance....')
    InfoDistance1 = abs(int64(cuts1Info) - int64(rsites1Info));
    InfoDistance2 = abs(int64(cuts2Info) - int64(rsites2Info));
    InfoDistance1Mask = (InfoDistance1 < cutToRsiteTh);
    InfoDistance2Mask = (InfoDistance2 < cutToRsiteTh);
    clear cuts1Info cuts2Info rsites1Info rsites2Info InfoDistance1 InfoDistance2;

    chrms1Info  = chrms1Info(InfoDistance1Mask);
    chrms2Info  = chrms2Info(InfoDistance2Mask);
    rfrag1Info  = rfrag1Info(InfoDistance1Mask);
    rfrag2Info  = rfrag2Info(InfoDistance2Mask);
    noFilteredInfo = sum(InfoDistance1Mask) + sum(InfoDistance2Mask);
    clear InfoDistance1Mask InfoDistance2Mask;

    %%
    disp('Counting Info restriction-fragments....')
    rfragCounts1 = countRsitesInfo(chrms1Info, chrms2Info, rfrag1Info, rfrag2Info, rsites, noChrs, chromosomes);

    %chrInfoDictionary = RDsignalfromRfrag(rfragCounts1, rsites, chrLengths, binSize, noChrs, chromosomes);
    clear chrms1Info chrms2Info rfrag1Info rfrag2Info;
end





%%%%%%%%%%%%%%%%%%%%
% ---------------- %
% (2) Non-Info reads
% ------------------
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
            rsites1Data  = h5read(hiclibDict,'/rsites1',startRead,countReads);
            rsites2Data  = h5read(hiclibDict,'/rsites2',startRead,countReads);

            %%
            [~, selfCircleMask, danglingEndMask, extraDanglingEndMask, ~, ~, ~] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);

            %% Non-informative Reads
            nonInfoMask = selfCircleMask | danglingEndMask | extraDanglingEndMask;
            chrms1NonInfo = [chrms1NonInfo; uint8(chrms1Data(nonInfoMask))];
            chrms2NonInfo = [chrms2NonInfo; uint8(chrms2Data(nonInfoMask))];
            cuts1NonInfo  = [cuts1NonInfo;  uint32(cuts1Data(nonInfoMask))];
            cuts2NonInfo  = [cuts2NonInfo;  uint32(cuts2Data(nonInfoMask))];
            rfrag1NonInfo = [rfrag1NonInfo; uint32(rfragIdxs1Data(nonInfoMask))];
            rfrag2NonInfo = [rfrag2NonInfo; uint32(rfragIdxs2Data(nonInfoMask))];
            rsites1NonInfo  = [rsites1NonInfo;uint32(rsites1Data(nonInfoMask))];
            rsites2NonInfo  = [rsites2NonInfo;uint32(rsites2Data(nonInfoMask))];        

        end
    end
    clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data rsites1Data rsites2Data;
    clear nonInfoMask selfCircleMask danglingEndMask extraDanglingEndMask; 


    disp('Removing Duplicated Non-Info Reads....')
    [uniqueMaskNonInfo] = removeDuplicatedMemoryEfficient(1,chrms1NonInfo, chrms2NonInfo, cuts1NonInfo, cuts2NonInfo, noChrs, chrLengths, memoryFootPrint);
    chrms1NonInfo = [chrms1NonInfo(uniqueMaskNonInfo)];
    chrms2NonInfo = [chrms2NonInfo(uniqueMaskNonInfo)];
    cuts1NonInfo  = [cuts1NonInfo(uniqueMaskNonInfo)];
    cuts2NonInfo  = [cuts2NonInfo(uniqueMaskNonInfo)];
    rfrag1NonInfo = [rfrag1NonInfo(uniqueMaskNonInfo)];
    rfrag2NonInfo = [rfrag2NonInfo(uniqueMaskNonInfo)];
    rsites1NonInfo  = [rsites1NonInfo(uniqueMaskNonInfo)];
    rsites2NonInfo  = [rsites2NonInfo(uniqueMaskNonInfo)];
    NDNonInfoReads = length(uniqueMaskNonInfo);
    clear uniqueMaskNonInfo;


    disp('Removing Non-Info Reads with large cut-rsite distance....')
    NonInfoDistance1 = min([abs(int64(cuts1NonInfo) - int64(rsites1NonInfo))], [abs(int64(cuts1NonInfo) - int64(rsites2NonInfo))]);
    NonInfoDistance2 = min([abs(int64(cuts2NonInfo) - int64(rsites1NonInfo))], [abs(int64(cuts2NonInfo) - int64(rsites2NonInfo))]);
    NonInfoDistance  = min(NonInfoDistance1, NonInfoDistance2);
    NonInfoDistanceMask = (NonInfoDistance < cutToRsiteTh);
    clear cuts1NonInfo cuts2NonInfo rsites1NonInfo rsites2NonInfo NonInfoDistance1 NonInfoDistance2 NonInfoDistance;

    chrms1NonInfo = chrms1NonInfo(NonInfoDistanceMask);
    chrms2NonInfo = chrms2NonInfo(NonInfoDistanceMask);
    rfrag1NonInfo = rfrag1NonInfo(NonInfoDistanceMask);
    rfrag2NonInfo = rfrag2NonInfo(NonInfoDistanceMask);
    noFilteredNonInfo = sum(NonInfoDistanceMask);
    clear NonInfoDistanceMask;

    disp('Counting Non-Info restriction-fragments....')
    rfragCounts2 = countRsitesNonInfo(chrms1NonInfo, chrms2NonInfo, rfrag1NonInfo, rfrag2NonInfo, rsites, noChrs, chromosomes);
    %chrNonInfoDictionary = RDsignalfromRfrag(rfragCounts2, rsites, chrLengths, binSize, noChrs, chromosomes);
    clear chrms1NonInfo chrms2NonInfo rfrag1NonInfo rfrag1NonInfo;
end




%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------%
% (3) Single-sided reads
% ----------------------
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
            rsites1Data  = h5read(hiclibDict,'/rsites1',startRead,countReads);
            rsites2Data  = h5read(hiclibDict,'/rsites2',startRead,countReads);

            %%
            [~, ~, ~, ~, singleSide1Mask, singleSide2Mask, ~] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);

            %% Single-sided Reads
            chrmsSingleSided  = [chrmsSingleSided; uint8(chrms1Data(singleSide1Mask)); uint8(chrms2Data(singleSide2Mask))];
            cutsSingleSided   = [cutsSingleSided;  uint32(cuts1Data(singleSide1Mask));  uint32(cuts2Data(singleSide2Mask))];
            rfragsSingleSided = [rfragsSingleSided; uint32(rfragIdxs1Data(singleSide1Mask));  uint32(rfragIdxs2Data(singleSide2Mask))];
    	    rsitesSingleSided = [rsitesSingleSided; uint32(rsites1Data(singleSide1Mask));  uint32(rsites2Data(singleSide2Mask))];
        end
    end
    clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data rsites1Data rsites2Data;
    clear singleSide1Mask singleSide2Mask; 

    disp('Removing Duplicated Single-sided Reads....')
    [uniqueMaskSingleSided] = removeDuplicatedMemoryEfficient(0,chrmsSingleSided, [], cutsSingleSided, [], noChrs, chrLengths, memoryFootPrint); 
    chrmsSingleSided   = [chrmsSingleSided(uniqueMaskSingleSided)];
    cutsSingleSided    = [cutsSingleSided(uniqueMaskSingleSided)];
    rfragsSingleSided  = [rfragsSingleSided(uniqueMaskSingleSided)];
    rsitesSingleSided  = [rsitesSingleSided(uniqueMaskSingleSided)];
    NDSingleSidedReads = length(uniqueMaskSingleSided);
    clear uniqueMaskSingleSided;

    disp('Removing Single-sided Reads with large cut-rsite distance....')
    SingleSidedDistance = abs(int64(cutsSingleSided) - int64(rsitesSingleSided));
    SingleSidedDistanceMask = (SingleSidedDistance < cutToRsiteTh);
    clear cutsSingleSided rsitesSingleSided SingleSidedDistance;

    chrmsSingleSided  = chrmsSingleSided(SingleSidedDistanceMask);
    rfragsSingleSided = rfragsSingleSided(SingleSidedDistanceMask);
    noFilteredSingleSided = sum(SingleSidedDistanceMask);

    clear SingleSidedDistanceMask;
    disp('Counting Single-sided restriction-fragments....')
    rfragCounts3 = countRsitesSingleSided(chrmsSingleSided, rfragsSingleSided, rsites, noChrs, chromosomes);
    %chrSSDictionary = RDsignalfromRfrag(rfragCounts3, rsites, chrLengths, binSize, noChrs, chromosomes);
    clear chrmsSingleSided rsitesSingleSided;
end




%%%%%%%%%%%%%%%
% -------------
% All reads
% -------------
disp('Computing restriction-fragments RD-signal....')

rfragCounts = containers.Map({1},{[]});
remove(rfragCounts, 1);
totalRfrags = 0;

for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    %
    if(readTypes == 1)
        rfragCounts(chrIndex) = rfragCounts1(chrIndex) + rfragCounts2(chrIndex) + rfragCounts3(chrIndex);
        chrRfrags = sum(rfragCounts(chrIndex));
    elseif(readTypes == 2)
        rfragCounts(chrIndex) = rfragCounts1(chrIndex);
        chrRfrags = sum(rfragCounts(chrIndex));
    elseif(readTypes == 3)
        rfragCounts(chrIndex) = rfragCounts2(chrIndex);
        chrRfrags = sum(rfragCounts(chrIndex));
    elseif(readTypes == 4)
        rfragCounts(chrIndex) = rfragCounts3(chrIndex);
        chrRfrags = sum(rfragCounts(chrIndex));
    end
    %
    totalRfrags = totalRfrags + chrRfrags;
end
obj.chrDictionary = RDsignalfromRfrag(rfragCounts, rsites, chrLengths, binSize, noChrs, chromosomes, largeRfrags);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------- Coverage Calculation ---------------%%
 
% readStatistics = [total-reads, double-side, informative, self-circle, dangling-end, extraDangling-end, single-side];
readStatistics = [totalReads, readsStatTot];
obj.readStatistics = readStatistics;

% Number-of-sides
if(readTypes == 1)
    obj.NDupReadStatistics = [NDInfoReads, NDNonInfoReads, NDSingleSidedReads];
    obj.nearRSitesSides = [noFilteredInfo, noFilteredNonInfo, noFilteredSingleSided];
    totalSides = noFilteredInfo + noFilteredNonInfo*2 + noFilteredSingleSided;
elseif(readTypes == 2)
    obj.NDupReadStatistics = [NDInfoReads];
    obj.nearRSitesSides = [noFilteredInfo];
    totalSides = noFilteredInfo;
elseif(readTypes == 3)
    obj.NDupReadStatistics = [NDNonInfoReads];
    obj.nearRSitesSides = [noFilteredNonInfo];    
    totalSides = noFilteredNonInfo*2;    
elseif(readTypes == 4)
    obj.NDupReadStatistics = [NDSingleSidedReads];
    obj.nearRSitesSides = [noFilteredSingleSided];
    totalSides = noFilteredSingleSided;
end

% total-coverage %
%%%%%%%%%%%%%%%%%%
% 1) Counting the sides of the reads
obj.dataCoverage = totalSides*readLength/genomeLengthBps; 

% 2) Counting restriction-fragments
%avgRfragLength = 4^obj.restrictionEnzymeLength;
%obj.dataCoverage = totalRfrags*avgRfragLength/genomeLengthBps; 






end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------- End "readsExtraction" ---------------------------------------------------%%



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis Sub-routines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------------------------------------------------------------------%%
%% Sub-routines for extracting the different types of reads.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = countRsitesInfo(chrms1, chrms2, rfrag1, rfrag2, rsites, noChrs, chromosomes)
%Each cut-site is translated to a restriction-fragment: we count a restriction-fragment for each cut. 
%Hiclib rfrag is started from 0. For example, rfrag = 0 means region from 0 to rsites(1). Our indices is from 1, so all rfrags are incremented first.

%%%%%%%%
chrDict = containers.Map({1},{[]});
remove(chrDict,1);
%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chrRsites = rsites(chrIndex);
    chrLengthRFragments = length(chrRsites) + 1;
    readCountsRsites = zeros(chrLengthRFragments,1);
    
    %%%
    chrIndexHiCLib = chrIndex - 1;
    chrms1Index = (chrms1 == chrIndexHiCLib);
    chrms2Index = (chrms2 == chrIndexHiCLib);
    rfrags = [rfrag1(chrms1Index); rfrag2(chrms2Index)];
    %
    rfrags = rfrags + 1;%since hiclib rfrag starts from 0.
    rfrags = rfrags(rfrags <= chrLengthRFragments);%to remove false-reads
    
    %%% Counting the restriction-fragments %%%
    rfrags = double(rfrags);
    [a,b] = hist(rfrags,unique(rfrags));
    readCountsRsites(b) = a;

    %noReads = length(rfrags);
    %for j = 1:noReads
    %    rfragment = rfrags(j);
    %    readCountsRsites(rfragment) = readCountsRsites(rfragment)+ 1;
    %end
    %%
    chrDict(chrIndex) = readCountsRsites;
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = countRsitesNonInfo(chrms1, chrms2, rfrag1, rfrag2, rsites, noChrs,chromosomes)
%Each restriction-fragment, intersected by the read fragment, is counted once regardless its type (DE, Extra DE, Self-Circle). 
%%%%%%%%
chrDict = containers.Map({1},{[]});
remove(chrDict,1);
%%%%%%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chrRsites = rsites(chrIndex);
    chrLengthRFragments = length(chrRsites) + 1;
    readCountsRsites = zeros(chrLengthRFragments,1);
    
    %%
    chrIndexHiCLib = chrIndex - 1;
    chrmsIndex = (chrms1 == chrIndexHiCLib);
    rfrags1 = rfrag1(chrmsIndex); 
    rfrags2 = rfrag2(chrmsIndex);
    %
    rfrags1 = rfrags1 + 1;%since hiclib rfrag starts from 0.
    rfrags2 = rfrags2 + 1;%since hiclib rfrag starts from 0.
    correctFragIndices = (rfrags1 <= chrLengthRFragments) & (rfrags2  <= chrLengthRFragments);
    rfrags1 = rfrags1(correctFragIndices);%to remove false-reads    
    rfrags2 = rfrags2(correctFragIndices);%to remove false-reads    
    %%

    noReads = length(rfrags1);
    for j = 1:noReads
        rfragment1 = rfrags1(j);
	    rfragment2 = rfrags2(j);
        readCountsRsites(rfragment1:rfragment2) = readCountsRsites(rfragment1:rfragment2)+ 1;
    end
    %%
    chrDict(chrIndex) = readCountsRsites;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = countRsitesSingleSided(chrms1, rfrag1, rsites, noChrs,chromosomes)
%Each restriction-fragment, intersected by the read fragment, is counted once regardless its type (DE, Extra DE, Self-Circle). 
%%%%%%%%
chrDict = containers.Map({1},{[]});
remove(chrDict,1);
%%%%%%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chrRsites = rsites(chrIndex);
    chrLengthRFragments = length(chrRsites) + 1;
    readCountsRsites = zeros(chrLengthRFragments,1);

    %%
    chrIndexHiCLib = chrIndex - 1;
    chrmsIndex = (chrms1 == chrIndexHiCLib);
    rfrags1 = rfrag1(chrmsIndex); 
    %
    rfrags1 = rfrags1 + 1;%since hiclib rfrag starts from 0.
    rfrags1 = rfrags1(rfrags1 <= chrLengthRFragments);%to remove false-reads    
    %%

    rfrags1 = double(rfrags1);
    [a,b] = hist(rfrags1,unique(rfrags1));
    readCountsRsites(b) = a;

    %noReads = length(rfrags1);
    %for j = 1:noReads
    %    rfragment1 = rfrags1(j);
    %    readCountsRsites(rfragment1) = readCountsRsites(rfragment1) + 1;
    %end
    %%
    chrDict(chrIndex) = readCountsRsites;
end

end


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
