function obj = plotFragmentLengthHistogram(obj)
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
noFiles = length(HDF5Files);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------- %
%  Non-Info reads: Dangling-end    %
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

%%%
disp('Extracting Non-informative Reads....')


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
    end
end
clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data rsites1Data rsites2Data;
clear nonInfoMask selfCircleMask danglingEndMask extraDanglingEndMask; 

%%% Dangling %%%
DanglingFragmentLength = [abs(int32(cuts1Dangling) - int32(cuts2Dangling))];

%%% Extra-Dangling %%%
ExtraDanglingFragmentLength = [abs(int32(cuts1ExDangling) - int32(cuts2ExDangling))];

%%%%%
clear cuts1Dangling cuts2Dangling rsites1Dangling rsites2Dangling;
clear cuts1ExDangling cuts2ExDangling rsites1ExDangling rsites2ExDangling;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------------------------------------------------ Figures -----------------------------------------------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%----------------------- Fragment length -----------------%%%%%%%%%%%%%
extraDanglingTh2 = max(ExtraDanglingFragmentLength)*0.99;
danglingTh2 = max(DanglingFragmentLength)*0.99;
maxTh2 = max(danglingTh2, extraDanglingTh2);
maxTh2 = min(maxTh2, 1200);
%%
fExtraDanglingFragmentLength = ExtraDanglingFragmentLength(ExtraDanglingFragmentLength < maxTh2);
fDanglingFragmentLength   = DanglingFragmentLength(DanglingFragmentLength < maxTh2);

%%%%%%
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


%%%%%%
figure;
subplot(1,2,1);
h1 = histogram(fDanglingFragmentLength,200);
h1.FaceColor = [0 0 0.6];
h1.EdgeColor = [0 0 0.6];
title('Dangling-end');
xlabel('Fragment length');
ylabel('Frequency');
xlim([0 maxTh2])
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';

subplot(1,2,2);
h2 = histogram(fExtraDanglingFragmentLength,200);
h2.FaceColor = [0 0 0.6];
h2.EdgeColor = [0 0 0.6];
title('Extra dangling-end');
xlabel('Fragment length');
ylabel('Frequency');
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
xlim([0 maxTh2])

%ff = strcat('-f',num2str(f2));
%hh = strcat(runDirectory, '/fragmetLength');
%print(ff,hh,'-dpng');
%savefig(hh);
%close All;


end
%%%


















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


