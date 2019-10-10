function [gainInfo] = nonInfoGain(obj)
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
%%%%%%%%%%%%%%%%%%%%
% ---------------- %
%  Non-Info reads  %
% -----------------%

disp('Extracting Non-informative Reads....')

%-- Reads Classification --%
noFiles = length(HDF5Files);

for fileNo = 1:noFiles

    %----- File Details ----%
    hiclibDict = cell2mat(HDF5Files(fileNo));
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
%
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
%
NonInfoDistance1 = min([abs(int64(cuts1NonInfo) - int64(rsites1NonInfo))], [abs(int64(cuts1NonInfo) - int64(rsites2NonInfo))]);
NonInfoDistance2 = min([abs(int64(cuts2NonInfo) - int64(rsites1NonInfo))], [abs(int64(cuts2NonInfo) - int64(rsites2NonInfo))]);
NonInfoDistance  = min(NonInfoDistance1, NonInfoDistance2);
NonInfoDistanceMask = (NonInfoDistance < cutToRsiteTh);
clear cuts1NonInfo cuts2NonInfo rsites1NonInfo rsites2NonInfo NonInfoDistance1 NonInfoDistance2 NonInfoDistance;
%
chrms1NonInfo = chrms1NonInfo(NonInfoDistanceMask);
chrms2NonInfo = chrms2NonInfo(NonInfoDistanceMask);
rfrag1NonInfo = rfrag1NonInfo(NonInfoDistanceMask);
rfrag2NonInfo = rfrag2NonInfo(NonInfoDistanceMask);
noFilteredNonInfo = sum(NonInfoDistanceMask);
clear NonInfoDistanceMask;
%
rfragNonInfo = countRsitesNonInfo(chrms1NonInfo, chrms2NonInfo, rfrag1NonInfo, rfrag2NonInfo, rsites, noChrs, chromosomes);
%%%%
noNonInfoRestFrag = sum(rfragNonInfo);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extracting informative Reads....')
%%%

%-- Reads Classification --%
noFiles = length(HDF5Files);

for fileNo = 1:noFiles

    %----- File Details ----%
    hiclibDict = cell2mat(HDF5Files(fileNo));
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

InfoDistance1 = abs(int64(cuts1Info) - int64(rsites1Info));
InfoDistance2 = abs(int64(cuts2Info) - int64(rsites2Info));
InfoDistance1Mask = (InfoDistance1 < cutToRsiteTh);
InfoDistance2Mask = (InfoDistance2 < cutToRsiteTh);
InfoDistanceMask  = (InfoDistance1Mask) & (InfoDistance2Mask);
clear cuts1Info cuts2Info rsites1Info rsites2Info InfoDistance1 InfoDistance2;

chrms1Info  = chrms1Info(InfoDistanceMask);
chrms2Info  = chrms2Info(InfoDistanceMask);
rfrag1Info  = rfrag1Info(InfoDistanceMask);
rfrag2Info  = rfrag2Info(InfoDistanceMask);
noFilteredInfo = sum(InfoDistanceMask);
clear InfoDistance1Mask InfoDistance2Mask;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Sampling info-reads %%%%%%%%%%%%%%%%%%%%%
%% 1st %%
%%%%%%%%%

noValidPairs1    = round(noNonInfoRestFrag);
selectedIndices = randi(noFilteredInfo,noValidPairs1,1); 
chrms1Info1  = chrms1Info(selectedIndices);
chrms2Info1  = chrms2Info(selectedIndices);
rfrag1Info1  = rfrag1Info(selectedIndices);
rfrag2Info1  = rfrag2Info(selectedIndices);
%%
rfragInfo1 = countRsitesInfo(chrms1Info1, chrms2Info1, rfrag1Info1, rfrag2Info1, rsites, noChrs, chromosomes);

%%%%%%%%%
%% 2nd %%
%%%%%%%%%

noValidPairs2    = round(noNonInfoRestFrag/2);
selectedIndices = randi(noFilteredInfo,noValidPairs2,1); 
chrms1Info2  = chrms1Info(selectedIndices);
chrms2Info2  = chrms2Info(selectedIndices);
rfrag1Info2  = rfrag1Info(selectedIndices);
rfrag2Info2  = rfrag2Info(selectedIndices);
%%
rfragInfo2 = countRsitesInfo(chrms1Info2, chrms2Info2, rfrag1Info2, rfrag2Info2, rsites, noChrs, chromosomes);

clear chrms1Info chrms2Info rfrag1Info rfrag2Info;
clear chrms1Info1 chrms2Info1 rfrag1Info1 rfrag2Info1;
clear chrms1Info2 chrms2Info2 rfrag1Info2 rfrag2Info2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  Analysis  %%%%%%%%%%%%%%%%%%%
rfragInfoA = rfragInfo1;
rfragInfoB = rfragNonInfo + rfragInfo2;
%%

noReadsA = sum(rfragInfoA);
noReadsB = sum(rfragInfoB);
%
meanA = mean(rfragInfoA);
meanB = mean(rfragInfoB);
%
stdA = std(rfragInfoA);
stdB = std(rfragInfoB);
%
noZerosA = sum(rfragInfoA == 0);
noZerosB = sum(rfragInfoB == 0);
%
[corrValue, ~]  = corr(rfragInfoA, rfragInfoB,'Type','Spearman');
%
gainInfo = [noReadsA, noReadsB, meanA, meanB, stdA, stdB, noZerosA, noZerosB, corrValue];
%%


%%%%%%%%%%%%%
%% Figures %%
%%%%%%%%%%%%%
format long;
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
    mkdir(dir);
end
%
runDirectory = strcat(dir,'/','figures');
if(exist(runDirectory,'dir') ~= 7)
    mkdir(runDirectory);
end
%%%
ThA = max(rfragInfoA)*0.90;
ThB = max(rfragInfoB)*0.90;
rfragInfoA = rfragInfoA(rfragInfoA < ThA);
rfragInfoB = rfragInfoB(rfragInfoB < ThB);

%%%%%%
%%%%%%
%f1 = randi(10000,1,1);
%figure;
%subplot(1,2,1);
%h1 = histogram(rfragInfoA, 200);
%h1.FaceColor = [0 0 0.6];
%h1.EdgeColor = [0 0 0.6];
%title('Info only');
%xlabel('Read per restriction-fragment');
%ylabel('Frequency');

%ax = gca;
%ax.FontSize = 10;
%ax.FontWeight = 'bold';
%%
%subplot(1,2,2);
%h2 = histogram(rfragInfoB, 200);
%h2.FaceColor = [0 0 0.6];
%h2.EdgeColor = [0 0 0.6];
%title('Info + Non-info');
%xlabel('Read per restriction-fragment');
%ylabel('Frequency');
%ax = gca;
%ax.FontSize = 10;
%ax.FontWeight = 'bold';

%ff = strcat('-f',num2str(f1));
%hh = strcat(runDirectory, '/nonInfoGain');
%print(ff,hh,'-dpng');
%savefig(hh);
%close All;







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
chrDict = [];
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
    chrDict = [chrDict; readCountsRsites];
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = countRsitesNonInfo(chrms1, chrms2, rfrag1, rfrag2, rsites, noChrs,chromosomes)
%Each restriction-fragment, intersected by the read fragment, is counted once regardless its type (DE, Extra DE, Self-Circle). 
%%%%%%%%
chrDict = [];
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
    chrDict = [chrDict; readCountsRsites];
end

end
