function obj = contactMapComputing(obj)
%Computing the interaction-frequencies between all genomic loci.
%Interaction-matrices will be save at the output directory "outputDirectory". 

outputDirectory = obj.outputDirectory;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------- chromosomes for Analysis -------------%%
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
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
readLength = obj.readLength;
cutToRsiteTh = obj.cutToRsiteThreshold;
rsites = obj.rsites;
%%%
binSize = obj.contactMapBinSize;
%%%
chrLengthsVector = cell2mat(values(chrLengths));
maxChrmsLength = max(chrLengthsVector);
genomeLengthBps = sum(chrLengthsVector);
clear chrLengthsVector;
%%%
readsStatTot = [0,0,0,0,0,0];%Double-side, Informative, Self-circle, Dangling-end, ExtraDangling-end, Single-side
totalReads = 0;
%%%
chrms1Info  = [];
chrms2Info  = [];
cuts1Info   = [];
cuts2Info   = [];
rfrag1Info  = [];
rfrag2Info  = [];
rsites1Info = [];
rsites2Info = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Reads Extraction & Duplicated Removal & Computing RD-signal---%

disp('Extracting informative-reads....')

%-- Reads Classification --%
noFiles = length(HDF5Files)

for fileNo = 1:noFiles
	%----- File Details ----%
	hiclibDict = cell2mat(HDF5Files(fileNo))
	inputInfo = h5info(hiclibDict);
	inputDimensions = inputInfo.Datasets.Dataspace;
	noReads = inputDimensions.Size;    
	totalReads = totalReads + noReads;
	%%
	for readOrder = 1:memoryFootPrint:noReads
	    startRead = readOrder;
	    stopRead  = min(startRead+memoryFootPrint-1,noReads);
	    countReads = stopRead - startRead + 1;

	    %%
	    chrms1Data     = h5read(hiclibDict,'/chrms1',startRead,countReads);
	    chrms2Data     = h5read(hiclibDict,'/chrms2',startRead,countReads);
	    cuts1Data      = h5read(hiclibDict,'/cuts1',startRead,countReads);
	    cuts2Data      = h5read(hiclibDict,'/cuts2',startRead,countReads);
	    rfragIdxs1Data = h5read(hiclibDict,'/rfragIdxs1',startRead,countReads);
	    rfragIdxs2Data = h5read(hiclibDict,'/rfragIdxs2',startRead,countReads);
	    strands1Data   = h5read(hiclibDict,'/strands1',startRead,countReads);
	    strands2Data   = h5read(hiclibDict,'/strands2',startRead,countReads);
	    rsites1Data    = h5read(hiclibDict,'/rsites1',startRead,countReads);
	    rsites2Data    = h5read(hiclibDict,'/rsites2',startRead,countReads);        

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
	    rsites1Info= [rsites1Info ; uint32(rsites1Data(infoMask))];
	    rsites2Info= [rsites2Info ; uint32(rsites2Data(infoMask))];
	    %% 
	end
end

clear chrms1Data chrms2Data cuts1Data cuts2Data rfragIdxs1Data rfragIdxs2Data strands1Data strands2Data rsites1Data rsites2Data;
clear infoMask;  

disp('Removing Duplicated Info Reads....')
[uniqueMaskInfo] = removeDuplicatedMemoryEfficient(1,chrms1Info, chrms2Info, cuts1Info, cuts2Info, noChrs, chrLengths, memoryFootPrint);
chrms1Info   = [chrms1Info(uniqueMaskInfo)];
chrms2Info   = [chrms2Info(uniqueMaskInfo)];
cuts1Info    = [cuts1Info(uniqueMaskInfo)];
cuts2Info    = [cuts2Info(uniqueMaskInfo)];
rfrag1Info   = [rfrag1Info(uniqueMaskInfo)];
rfrag2Info   = [rfrag2Info(uniqueMaskInfo)];
rsites1Info  = [rsites1Info(uniqueMaskInfo)];
rsites2Info  = [rsites2Info(uniqueMaskInfo)];
clear uniqueMaskInfo;

disp('Removing Info-reads with large cut-rsite distance....')
InfoDistance1 = abs(int64(cuts1Info) - int64(rsites1Info));
InfoDistance2 = abs(int64(cuts2Info) - int64(rsites2Info));
InfoDistanceMask = (InfoDistance1 < cutToRsiteTh) & (InfoDistance2 < cutToRsiteTh);
clear cuts1Info cuts2Info rsites1Info rsites2Info InfoDistance1 InfoDistance2;

chrms1Info  = chrms1Info(InfoDistanceMask);
chrms2Info  = chrms2Info(InfoDistanceMask);
rfrag1Info  = rfrag1Info(InfoDistanceMask);
rfrag2Info  = rfrag2Info(InfoDistanceMask);
noValidReadPairs = sum(InfoDistanceMask);
clear InfoDistanceMask;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Computing the contact-maps & saving them in the outputDirectory---%
disp('Computing the contact maps ......')

computeInteractionFreq(obj, chrms1Info, chrms2Info, rfrag1Info, rfrag2Info, rsites, noChrs, chrLengths, chromosomes);
obj.noValidReadPairs = noValidReadPairs;
clear chrms1Info chrms2Info rfrag1Info rfrag2Info;




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------------- End "contactMap Computing" ------------------%%%




















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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrDict] = computeInteractionFreq(obj, chrms1, chrms2, rfrag1, rfrag2, rsites, noChrs, chrLengths, chromosomes)
%We count reads based on the middle-bp of restriction fragments.
%Hiclib rfrag is started from 0. For example, rfrag = 0 means region from 0 to rsites(1). Our indices is from 1, so all rfrags are incremented first.

binSize = obj.contactMapBinSize;
%%%%%%%%%%%%%%%%%%%
for i  = 1:1:noChrs
	%
	chr1Index = chromosomes(i)
	chr1Rsites = rsites(chr1Index);
	chr1LengthRFragments = length(chr1Rsites) + 1;
	chr1LengthBps   = chrLengths(chr1Index);
	chr1LengthBins  = ceil(chr1LengthBps/binSize);
	%
    chr1RfragStart = [1; chr1Rsites];
    chr1RfragEnd   = [chr1Rsites-1; chr1LengthBps]; 
	chr1RfragMid   = round((chr1RfragStart+chr1RfragEnd)/2);
	%
	% cis only condition
	if(obj.cisOnly == 1)
		lastChr = i;
	else
		lastChr = noChrs;
	end		
	%%%%%%%%%%%%%%%%%%%
	for j = i:1:lastChr
		%
		chr2Index = chromosomes(j)
		chr2Rsites = rsites(chr2Index);
		chr2LengthRFragments = length(chr2Rsites) + 1;
		chr2LengthBps  = chrLengths(chr2Index);
		chr2LengthBins  = ceil(chr2LengthBps/binSize);
		%
		chr2RfragStart = [1; chr2Rsites];
		chr2RfragEnd   = [chr2Rsites-1; chr2LengthBps]; 
		chr2RfragMid   = round((chr2RfragStart+chr2RfragEnd)/2);
		%
		chr1IndexHiCLib = chr1Index - 1;
		chr2IndexHiCLib = chr2Index - 1;


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%--------- Reads-selection --------%%
		if(chr1IndexHiCLib == chr2IndexHiCLib)
		% cis-matrix
			%%%
			chrmsIndex = (chrms1 == chr1IndexHiCLib & chrms2 == chr2IndexHiCLib);
			%
			rfrags1 = rfrag1(chrmsIndex);
			rfrags1 = rfrags1 + 1;%since hiclib rfrag starts from 0.
			rfrags2 = rfrag2(chrmsIndex);
			rfrags2 = rfrags2 + 1;%since hiclib rfrag starts from 0.
			%
			rfragsIndex = (rfrags1 <= chr1LengthRFragments & rfrags2 <= chr2LengthRFragments);%to remove false-reads	
			rfrags1 = rfrags1(rfragsIndex);		
			rfrags2 = rfrags2(rfragsIndex);
			%%%
			rfrags1Mid = chr1RfragMid(rfrags1);
			rfrags2Mid = chr2RfragMid(rfrags2);		
			%
			rfrags1Bin = floor(rfrags1Mid/binSize) + 1;%Chr coordinates start from 0.
			rfrags2Bin = floor(rfrags2Mid/binSize) + 1;%Chr coordinates start from 0.			
			%
			%sorting bins	
			rfragsBins  = [rfrags1Bin,rfrags2Bin];
			rfragsBinsSorted  = sort(rfragsBins,2);   
			rfrags1Bin = rfragsBinsSorted(:,1);
			rfrags2Bin = rfragsBinsSorted(:,2);
			clear rfragsBins rfragsBinsSorted;				
		else
		% trans-matrix
			%%% 
			chrmsIndexA = (chrms1 == chr1IndexHiCLib & chrms2 == chr2IndexHiCLib);
			chrmsIndexB = (chrms1 == chr2IndexHiCLib & chrms2 == chr1IndexHiCLib);
			%
			rfrags1 = [rfrag1(chrmsIndexA); rfrag2(chrmsIndexB)];
			rfrags1 = rfrags1 + 1;%since hiclib rfrag starts from 0.
			rfrags2 = [rfrag2(chrmsIndexA); rfrag1(chrmsIndexB)];
			rfrags2 = rfrags2 + 1;%since hiclib rfrag starts from 0.
			%
			rfragsIndex = (rfrags1 <= chr1LengthRFragments & rfrags2 <= chr2LengthRFragments);%to remove false-reads	
			rfrags1 = rfrags1(rfragsIndex);		
			rfrags2 = rfrags2(rfragsIndex);
			%%%
			rfrags1Mid = chr1RfragMid(rfrags1);
			rfrags2Mid = chr2RfragMid(rfrags2);		
			%
			rfrags1Bin = floor(rfrags1Mid/binSize) + 1;%Chr coordinates start from 0.
			rfrags2Bin = floor(rfrags2Mid/binSize) + 1;%Chr coordinates start from 0.						
		end
		%%%


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%% Matrix-computing %%%%%%%%	
		noInteractions = length(rfrags1Bin);
		intFreq = zeros(chr1LengthBins,chr2LengthBins);
		%
		for k = 1:1:noInteractions
			intFreq(rfrags1Bin(k), rfrags2Bin(k)) = intFreq(rfrags1Bin(k), rfrags2Bin(k)) + 1;
		end
		%%% 


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% Saving Interaction-matrices %%
		format long;
		dir = obj.outputDirectory;
		if(exist(dir,'dir') ~= 7)
			mkdir(dir);
		end
		%
		runDirectory = strcat(dir,'/','rawContactMap_binSize',int2str(obj.contactMapBinSize));
		if(exist(runDirectory,'dir') ~= 7)
			mkdir(runDirectory);
		end
		%
		cisDirectory = strcat(runDirectory,'/','cis_matrics');
		transDirectory = strcat(runDirectory,'/','trans_matrics');
		if(exist(cisDirectory,'dir') ~= 7)
			mkdir(cisDirectory);
		end
		if(exist(transDirectory,'dir') ~= 7)
			mkdir(transDirectory);
		end
		%%%
		if(chr1Index ~= 23)
			chr1Name = strcat('chr',int2str(chr1Index));
		else
			chr1Name = 'chrX';
		end
		%
		if(chr2Index ~= 23)
			chr2Name = strcat('chr',int2str(chr2Index));
		else
			chr2Name = 'chrX';
		end
		%%%
		if(chr1Index == chr2Index)
			filePath = strcat(cisDirectory, '/', chr1Name, '.txt');
		else
			filePath = strcat(transDirectory, '/', chr1Name, '_', chr2Name, '.txt');

		end
		%
		% remove diagonal reads
		M = intFreq;
		intFreq = M - diag(diag(M));
		%
		[row col v] = find(intFreq);
		dlmwrite(filePath, [row col v], 'delimiter','\t');
	%%%	
	end
end



%%%
end