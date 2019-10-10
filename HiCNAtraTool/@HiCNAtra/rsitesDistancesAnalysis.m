function [distancePercentages] = rsitesDistancesAnalysis(obj,plotCond)
%Analysis of the distances between rsites and cuts for different read-types


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chrLengthsVector = cell2mat(values(chrLengths));
maxChrmsLength = max(chrLengthsVector);
genomeLengthBps = sum(chrLengthsVector);
clear chrLengthsVector;
noFiles = length(HDF5Files);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% -------------%
% (1) Info reads
% -------------%
cuts1Info  = [];
cuts2Info  = [];
rsites1Info = [];
rsites2Info = [];

%%%
disp('Extracting informative Reads....')
%%%

%-- Reads Classification --%

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
        [infoMask, ~, ~, ~, ~, ~, readsStat] = readsInfoVsNonInfo(chrms1Data, chrms2Data, cuts1Data, cuts2Data, rfragIdxs1Data, rfragIdxs2Data, strands1Data, strands2Data, maximumMoleculeLength);

        %% Informative Reads (before distance-filtering)
        cuts1Info  = [cuts1Info ; uint32(cuts1Data(infoMask))];
        cuts2Info  = [cuts2Info ; uint32(cuts2Data(infoMask))];
        rsites1Info  = [rsites1Info ; uint32(rsites1Data(infoMask))];
        rsites2Info  = [rsites2Info ; uint32(rsites2Data(infoMask))];

        %% 

    end
end
%%%

clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data rsites1Data rsites2Data;
clear infoMask;  
%%%

% Distance Calculation %
InfoDistance = [abs(int32(cuts1Info) - int32(rsites1Info)); abs(int32(cuts2Info) - int32(rsites2Info))];
clear cuts1Info cuts2Info rsites1Info rsites2Info;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------- %
% (2) Non-Info reads: Dangling-end %
% -------------------------------- %
cuts1Dangling   = [];
cuts2Dangling   = [];
rsites1Dangling = [];
rsites2Dangling = [];
%
cuts1ExDangling   = [];
cuts2ExDangling   = [];
rsites1ExDangling = [];
rsites2ExDangling = [];
%
cuts1SelfCircle   = [];
cuts2SelfCircle   = [];
rsites1SelfCircle = [];
rsites2SelfCircle = [];

%%%
disp('Extracting Non-informative Reads....')

%-- Reads Classification --%

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
        %%
        cuts1Dangling   = [cuts1Dangling; uint32(cuts1Data(danglingEndMask))];
        cuts2Dangling   = [cuts2Dangling; uint32(cuts2Data(danglingEndMask))];
        rsites1Dangling = [rsites1Dangling; uint32(rsites1Data(danglingEndMask))];
        rsites2Dangling = [rsites2Dangling; uint32(rsites2Data(danglingEndMask))];
        %
        cuts1ExDangling   = [cuts1ExDangling; uint32(cuts1Data(extraDanglingEndMask))];
        cuts2ExDangling   = [cuts2ExDangling; uint32(cuts2Data(extraDanglingEndMask))];
        rsites1ExDangling = [rsites1ExDangling; uint32(rsites1Data(extraDanglingEndMask))];
        rsites2ExDangling = [rsites2ExDangling; uint32(rsites2Data(extraDanglingEndMask))];
        %
        cuts1SelfCircle   = [cuts1SelfCircle; uint32(cuts1Data(selfCircleMask))];
        cuts2SelfCircle   = [cuts2SelfCircle; uint32(cuts2Data(selfCircleMask))];
        rsites1SelfCircle = [rsites1SelfCircle; uint32(rsites1Data(selfCircleMask))];
        rsites2SelfCircle = [rsites2SelfCircle; uint32(rsites2Data(selfCircleMask))];
    end
end
clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data rsites1Data rsites2Data;
clear nonInfoMask selfCircleMask danglingEndMask extraDanglingEndMask; 

%%% Dangling %%%
DanglingDistance1 = min([abs(int64(cuts1Dangling) - int64(rsites1Dangling))], [abs(int64(cuts1Dangling) - int64(rsites2Dangling))]);
DanglingDistance2 = min([abs(int64(cuts2Dangling) - int64(rsites1Dangling))], [abs(int64(cuts2Dangling) - int64(rsites2Dangling))]);
DanglingDistance  = [DanglingDistance1; DanglingDistance2];
DanglingFragmentLength = [abs(int32(cuts1Dangling) - int32(cuts2Dangling))];

%%% Extra-Dangling %%%
ExtraDanglingDistance = [abs(int32(cuts1ExDangling) - int32(rsites1ExDangling)); abs(int32(cuts2ExDangling) - int32(rsites2ExDangling))];
ExtraDanglingFragmentLength = [abs(int32(cuts1ExDangling) - int32(cuts2ExDangling))];

%%% Self-Circle %%%
SelfCircleDistance = [abs(int32(cuts1SelfCircle) - int32(rsites1SelfCircle)); abs(int32(cuts2SelfCircle) - int32(rsites2SelfCircle))];

%%%%%
clear cuts1Dangling cuts2Dangling rsites1Dangling rsites2Dangling;
clear cuts1ExDangling cuts2ExDangling rsites1ExDangling rsites2ExDangling;
clear cuts1SelfCircle cuts2SelfCircle rsites1SelfCircle rsites2SelfCircle;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prcInfo = 100*sum(InfoDistance < maximumMoleculeLength)/length(InfoDistance);
prcDangling = 100*sum(DanglingDistance < maximumMoleculeLength)/length(DanglingDistance);
prcExDangling = 100*sum(ExtraDanglingDistance < maximumMoleculeLength)/length(ExtraDanglingDistance);
prcSelfCircle = 100*sum(SelfCircleDistance < maximumMoleculeLength)/length(SelfCircleDistance);
distancePercentages = [prcInfo, prcDangling, prcExDangling, prcSelfCircle];
%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------------------------------------------------ Figures -----------------------------------------------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%----- Distance between cut-site & restriction-site ------%%%%%%%%%%%%%
infoTh = max(InfoDistance)*0.99;
extraDanglingTh = max(ExtraDanglingDistance)*0.99;
danglingTh = max(DanglingDistance)*0.99;
selfCircleTh = max(SelfCircleDistance)*0.99;
maxTh = max([infoTh, extraDanglingTh, danglingTh, selfCircleTh]);
maxTh = min(maxTh,3000);
%%
fInfoDistance = InfoDistance(InfoDistance < maxTh);
fExtraDanglingDistance = ExtraDanglingDistance(ExtraDanglingDistance < maxTh);
fDanglingDistance   = DanglingDistance(DanglingDistance < maxTh);
fSelfCircleDistance = SelfCircleDistance(SelfCircleDistance < maxTh);

%%%%%%%%%%%%%%%%%
if(plotCond == 1)
    f1 = randi(10000,1,1);
    figure;
    h1 = histogram(fInfoDistance,100,'Normalization','probability');
    h1.FaceColor = [1 0 0];
    h1.LineWidth = 0.2;
    hold on;
    h2 = histogram(fDanglingDistance,100, 'Normalization','probability');
    h2.FaceColor = [0 1 0];
    h2.LineWidth = 0.2;
    hold on;
    h3 = histogram(fExtraDanglingDistance,100, 'Normalization','probability');
    h3.FaceColor = [0 0 1];
    h3.LineWidth = 0.2;
    hold on;
    h4 = histogram(fSelfCircleDistance,100, 'Normalization','probability');
    h4.FaceColor = [0.6 0.24 0.6];
    h4.LineWidth = 0.2;

    xlabel('Cut-site to restriction-site');
    ylabel('Frequency');
    lgd = legend('Valid read-pairs','Dangling-end', 'Extra dangling-end', 'Self-circle');
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    lgd.FontSize = 8;

    %ff = strcat('-f',num2str(f1));
    %hh = strcat(runDirectory, '/readTypeAnalysis');
    %print(ff,hh,'-dpng');
    %savefig(hh);
end



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
