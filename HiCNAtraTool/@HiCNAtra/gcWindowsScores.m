function obj = gcWindowsScores(obj)
%Computing the gc scores for the gc-window around the restriction-sites.



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



%%----------------- GC windows ------------------%%
rsites = obj.rsites;
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
%
gcCalculationMethod = obj.gcCalculationMethod;
gcWindow = obj.gcWindow;
gcWindsFolder = obj.gcWindsFolder;
referenceGenomeFolder = obj.referenceGenomeFolder;
%

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
            if(strcmp(GC_Contents(j),'NA')==0)
                segmentStart = (j-1)*100+1;
                segmentStop  = min(j*100,chrLengthBps);
                gcBps(segmentStart:segmentStop)= str2num(cell2mat(GC_Contents(j)));
            end
        end


        gcWindsWindowLeft  = max([chrRsites - gcWindow, chrRsites-1],1);
        gcWindsWindowRight = min([chrRsites, chrRsites + gcWindow - 1], chrLengthBps);

        gcWindsWindowScores = zeros(length(gcWindsWindowLeft),2);
        for i = 1:length(gcWindsWindowLeft)

            %% Left block
            gcLeftBlock = gcBps(gcWindsWindowLeft(i,1):gcWindsWindowLeft(i,2));
            gcLeftRatio = nanmean(gcLeftBlock);   

            %% Right block
            gcRightBlock = gcBps(gcWindsWindowRight(i,1):gcWindsWindowRight(i,2));
            gcRightRatio = nanmean(gcRightBlock);        

            gcWindsWindowScores(i,1) = gcLeftRatio;
            gcWindsWindowScores(i,2) = gcRightRatio;
        end


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
        %
        gcWindsWindowLeft  = max([chrRsites - gcWindow, chrRsites-1],1);
        gcWindsWindowRight = min([chrRsites, chrRsites + gcWindow - 1], chrLengthBps);       

        gcWindsWindowScores = zeros(length(gcWindsWindowLeft),2);
        %
        for j = 1:length(gcWindsWindowLeft)

            %% Left block
            gcLeftBlock = chrSequence(gcWindsWindowLeft(j,1):gcWindsWindowLeft(j,2));
            % find the GC and AT content
            gc = (sum(gcLeftBlock == 'G' | gcLeftBlock == 'C' | gcLeftBlock == 'g' | gcLeftBlock == 'c'));
            at = (sum(gcLeftBlock == 'A' | gcLeftBlock == 'T' | gcLeftBlock == 'a' | gcLeftBlock == 't'));
            % calculate the ratio of GC to the total known nucleotides
            gcLeftRatio = gc/(gc+at);   

            %% Right block
            gcRightBlock = chrSequence(gcWindsWindowRight(j,1):gcWindsWindowRight(j,2));
            % find the GC and AT content
            gc = (sum(gcRightBlock == 'G' | gcRightBlock == 'C' | gcRightBlock == 'g' | gcRightBlock == 'c'));
            at = (sum(gcRightBlock == 'A' | gcRightBlock == 'T' | gcRightBlock == 'a' | gcRightBlock == 't'));
            % calculate the ratio of GC to the total known nucleotides
            gcRightRatio = gc/(gc+at);  

            gcWindsWindowScores(j,1) = gcLeftRatio;
            gcWindsWindowScores(j,2) = gcRightRatio;
        end
    end


    %%------------------ output dictionary ---------------%%
    obj.GCWindsWindowScores(chrIndex) = gcWindsWindowScores;
end    



end
%%%
