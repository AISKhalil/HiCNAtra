function obj = gcBinFeatureCalculater(obj)
%Computing the gc feature per bin.


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



%%----------- gc feature (window-based) -----------%%
binSize    = obj.binSize;
rsites     = obj.rsites;
chrNames   = obj.chrNames;
chrLengths = obj.chrLengths;
%
gcCalculationMethod = obj.gcCalculationMethod;
gcWindsFolder = obj.gcWindsFolder;
referenceGenomeFolder = obj.referenceGenomeFolder;


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
        gc100bins = nan(length(GC_Contents),1);
        for j = 1:length(GC_Contents)
            if(strcmp(GC_Contents(j),'NA')==0)
                gc100bins(j)= str2num(cell2mat(GC_Contents(j)));
            end
        end

        n = binSize/100;
        a = gc100bins;
        gcBins = arrayfun(@(k) nanmean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
        obj.chrGCTracks(chrIndex) = gcBins;


    %%%%%%%%%%%%%%%%%%%%%%%%%
    % UCSC reference genome %
    elseif(gcCalculationMethod == 2)

        FASTAfileName = strcat(referenceGenomeFolder, chromosome, '.fa');
        fidIn = fopen(FASTAfileName,'r');
        rawData = fread(fidIn,'*char')';
        fclose(fidIn);
        %
        firstBlock = rawData(1:100);
        NL = regexp(firstBlock, '[\n]');
        firstCharIndex = NL(1)+1;% to remove first line
        filteredData = rawData(firstCharIndex:end);
        chrSequence  = strrep(filteredData,sprintf('\n'),'');
        %
        chrLengthBins = ceil(chrLengthBps/binSize);
        gcBins = zeros(chrLengthBins,1);
        gcBinsStart = 1:binSize:chrLengthBps;
        gcBinsStop  = [binSize:binSize:chrLengthBps,chrLengthBps];
        %

        for j = 1:length(gcBins)

            %% gc block
            gcBlock = chrSequence(gcBinsStart(j):gcBinsStop(j));
            % find the GC and AT content
            gc = (sum(gcBlock == 'G' | gcBlock == 'C' | gcBlock == 'g' | gcBlock == 'c'));
            at = (sum(gcBlock == 'A' | gcBlock == 'T' | gcBlock == 'a' | gcBlock == 't'));
            % calculate the ratio of GC to the total known nucleotides
            gcRatio = gc/(gc+at);   
            gcBins(j) = gcRatio;
        end
        obj.chrGCTracks(chrIndex) = gcBins;
    end
end    



end
