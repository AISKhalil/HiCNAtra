function obj = darkRegionsExtraction(obj)
%Extracting the black-listed regions, centromeres, telomeres, and gap-regions.


%%---------------- chromosomes for Analysis ----------------%%
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


%%---------------------- Dark Region ----------------------%%
binSize = obj.binSize;
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
%
blackListFile = obj.blackListFile;
gapFile = obj.gapFile;
centromeresFile = obj.centromeresFile;
telomeresFile = obj.telomeresFile;
%
chrCentroTeloBoundaries  = obj.chrCentroTeloBoundaries;
chrBlackListedBoundaries = obj.chrBlackListedBoundaries;
chrGapBoundaries = obj.chrGapBoundaries;
%
uncertaintyDistanceToCentro = ceil(obj.uncertaintyDistanceToCentro/obj.binSize);%number of normal bins (lower-threshold < RD-values < upper-threshold), around the centromeres, that are considered as uncertainty regions.
uncertaintyDistanceToTelo = ceil(obj.uncertaintyDistanceToTelo/obj.binSize);  %number of normal bins (lower-threshold < RD-values < upper-threshold), around the telomeres, that are considered as uncertainty regions.
%

if(length(obj.chrEffectiveLengthTracks.keys) > 0)
    chrDict = obj.chrDictionary;
else
    disp('Error: compute RD-signal first');
end



%%------------------ Filters threshold -------------------%%
%++Telomeres Threshold (for normal bins)
thrLT = 3;
thrUT = 97;
%++Centromeres Threshold (for normal bins)
thrLC = 3;
thrUC = 97;


%%%%%%%%%%%%%%%%%%%%%%% Input Files %%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------%
%%%% -- Black-list regions -- %%%%
fileID = fopen(blackListFile,'r');
blacklistRegions = textscan(fileID, '%s%f%f%s%f');
chrBlackName = blacklistRegions{1};
chrBlackStart = blacklistRegions{2};
chrBlackEnd = blacklistRegions{3};
fclose(fileID);

%%%% ------- gapFile -------- %%%%
fileID = fopen(gapFile,'r');
gapRegions = textscan(fileID, '%f%s%f%f%f%s%f%s%s');
chrGapName  = gapRegions{2};
chrGapStart = gapRegions{3};
chrGapEnd   = gapRegions{4};
fclose(fileID);

%%%% Centromeres & Telomeres %%%%
fileID = fopen(centromeresFile,'r');
centromeresRegions = textscan(fileID, '%f%s%f%f%f%s%f%s%s');
chrCentroName = centromeresRegions{2};
chrCentroStart = centromeresRegions{3};
chrCentroEnd = centromeresRegions{4};
fclose(fileID);

%%%%
fileID = fopen(telomeresFile,'r');
telomeresRegions = textscan(fileID, '%f%s%f%f%f%s%f%s%s');
chrTeloName = telomeresRegions{2};
chrTeloStart = telomeresRegions{3};
chrTeloEnd = telomeresRegions{4};
fclose(fileID);


%%%%%%%%%%%%%%%%%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(chrIndex));
    chrLengthBps = chrLengths(chrIndex);
    targetChrLength = ceil(chrLengthBps/binSize);
    targetChrData = chrDict(chrIndex);
	%
	minRegionSizeToFilter = 0.25*obj.binSize;
	
	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gap Bins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    targetChrGapStart = [];
    targetChrGapEnd = [];
    for j=1:length(chrGapName)
        if strcmp(chromosome,chrGapName(j))== 1
			if(chrGapEnd(j) - chrGapStart(j) > minRegionSizeToFilter)
				chrGapStartRegion = max(floor(chrGapStart(j)/binSize)+1,1);
				chrGapEndRegion = ceil(chrGapEnd(j)/binSize);
				%
				targetChrGapStart = [targetChrGapStart; chrGapStartRegion];
				targetChrGapEnd   = [targetChrGapEnd  ; chrGapEndRegion];
			end
        end
    end
    chrGapBoundaries(chrIndex) = [targetChrGapStart, targetChrGapEnd];


    %%%%%%%%%%%%%%%%%%%%%%%%% Black-Listed Bins %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    %Filtering the black-listed regions based on the RD of surrounding bins
    targetChrBlStart = [];
    targetChrBlEnd = [];
    for j=1:length(chrBlackName)
        if strcmp(chromosome,chrBlackName(j))== 1
			if(chrBlackEnd(j) - chrBlackStart(j) > minRegionSizeToFilter)
				chrBlackStartRegion = max(floor(chrBlackStart(j)/binSize)+1,1);
				chrBlackEndRegion = ceil(chrBlackEnd(j)/binSize); 
				%
				targetChrBlStart = [targetChrBlStart; chrBlackStartRegion];
				targetChrBlEnd = [targetChrBlEnd; chrBlackEndRegion]; 
			end
        end
    end    
    chrBlackListedBoundaries(chrIndex) = [targetChrBlStart,targetChrBlEnd];


    %%%%%%%%%%%%%%%%%%%%%% Centromere/Telomere Bins %%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    centroTeloBins = ones(targetChrLength,1);
    %%%
    targetChrCentroStart = [];
    targetChrCentroEnd = [];
    for j=1:length(chrCentroName)
        if strcmp(chromosome,chrCentroName(j))== 1
            targetChrCentroStart = max(floor(chrCentroStart(j)/binSize)+1,1);
            targetChrCentroEnd = ceil(chrCentroEnd(j)/binSize);
        end
    end
    %%%
    targetChrTeloStart = [];
    targetChrTeloEnd = [];
    for j=1:length(chrTeloName)
        if strcmp(chromosome,chrTeloName(j))== 1
            TeloStartTemp = max(floor(chrTeloStart(j)/binSize)+1,1);
            if TeloStartTemp == 0
                TeloStartTemp = 1;
            end
            TeloEndTemp = ceil(chrTeloEnd(j)/binSize);
            if TeloEndTemp > targetChrLength
                TeloEndTemp = targetChrLength;%%%%
            end            
            targetChrTeloStart = [targetChrTeloStart; TeloStartTemp];
            targetChrTeloEnd = [targetChrTeloEnd; TeloEndTemp];
        end
    end
    %%% Chr17 Telomeres %%%
    if (length(targetChrTeloStart)~=2)
        targetChrTeloStart = [1,targetChrLength-1];
        targetChrTeloEnd = [2,targetChrLength];
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Extend the centromeres %%%%%%%%%
    countTh = uncertaintyDistanceToCentro;
    %%%%% Centromeres Extension %%%%%
    firstBin = targetChrTeloEnd(1)+1;
    lastBin  = targetChrCentroStart;
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample(sample ~= 0);   
    if(isempty(nonZeroData))
        targetChrCentroExStart= firstBin;
    else
        extensionTh = prctile(nonZeroData,thrLC);%%% threshold for extension
        extensionThUp = prctile(nonZeroData,thrUC);
        count = 0;
        for j = lastBin:-1:firstBin
            if(targetChrData(j)<= extensionTh ||targetChrData(j)>= extensionThUp)
                centroTeloBins(j)= 0;
            elseif (count == countTh)
                targetChrCentroExStart= j;
                break;
            else
                centroTeloBins(j)= 0;
                count = count +1;
            end
        end
        if(count ~=countTh)
           targetChrCentroExStart= firstBin;
        end
    end
    %%%%%%%%
    %%%%%%%%
    firstBin = targetChrCentroEnd;
    lastBin  = targetChrTeloStart(2)-1;
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample(sample ~= 0);
    if(isempty(nonZeroData))
       targetChrCentroExEnd= lastBin;
    else
        extensionTh = prctile(nonZeroData,thrLC);%%% threshold for extension     
        extensionThUp = prctile(nonZeroData,thrUC);   
        count = 0;    
        for j = firstBin:1:lastBin
            if(targetChrData(j) <= extensionTh ||targetChrData(j) >= extensionThUp)
                centroTeloBins(j)= 0;
            elseif (count == countTh)
                targetChrCentroExEnd= j;
                break;
            else
                centroTeloBins(j)= 0;
                count = count +1;           
            end
        end
        if(count ~=countTh)
           targetChrCentroExEnd= lastBin;
        end
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Telomere#1 Extension %%%%%
    countTh = uncertaintyDistanceToTelo;
    firstBin = targetChrTeloEnd(1);
    lastBin = targetChrCentroExStart-1;
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample(sample ~= 0);
    extensionTh = prctile(nonZeroData,thrLT);%%% threshold for extension     
    extensionThUp = prctile(nonZeroData,thrUT);  
    
    targetChrTelo1ExStart = targetChrTeloStart(1);
    count = 0;       
    for j = firstBin:1:lastBin
        if(targetChrData(j) <= extensionTh ||targetChrData(j)>= extensionThUp)
            centroTeloBins(j)= 0;
        elseif (count == countTh)
            targetChrTelo1ExEnd= j;
            break;
        else
            centroTeloBins(j)= 0;
            count = count + 1;           
        end
    end
    if(count ~=countTh)
       targetChrTelo1ExEnd= lastBin;
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Telomere#2 Extension %%%%%
    countTh = uncertaintyDistanceToTelo;
    firstBin = targetChrCentroExEnd+1;
    lastBin  = targetChrTeloStart(2);
    sample = targetChrData(firstBin:lastBin);
    nonZeroData = sample;%(sample ~= 0);
    extensionTh = prctile(nonZeroData,thrLT);%%% threshold for extension     
    extensionThUp = prctile(nonZeroData,thrUT);     
    
    count = 0;        
    for j = lastBin:-1:firstBin
        if(targetChrData(j)<= extensionTh ||targetChrData(j)>= extensionThUp)
            centroTeloBins(j)= 0;
        elseif (count == countTh)
            targetChrTelo2ExStart= j;
            break;
        else
            centroTeloBins(j)= 0;
            count = count + 1;           
        end
    end
    if(count ~=countTh)
       targetChrTelo2ExStart= firstBin;
    end
    targetChrTelo2ExEnd = targetChrTeloEnd(2); 
    %
    centroTeloBoundaries  = [targetChrCentroExStart, targetChrCentroExEnd; targetChrTelo1ExStart, targetChrTelo1ExEnd; targetChrTelo2ExStart, targetChrTelo2ExEnd];
    chrCentroTeloBoundaries(chrIndex) = centroTeloBoundaries;

end    



end
