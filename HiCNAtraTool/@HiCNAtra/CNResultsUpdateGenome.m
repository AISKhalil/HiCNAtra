function [] = CNResultsUpdateGenome(obj)
%%Writing CNAtra results to output files including CNV regions and IBs for all chromosomes.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------- Output Directory ------------ %%
format long;
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
	mkdir(dir);
end

runDirectory = strcat(dir,'/','CNVs_binSize',int2str(obj.binSize));
if(exist(runDirectory,'dir') ~= 7)
	mkdir(runDirectory);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------ Output Files ------------------------%%
%%
statisticsFilePath = strcat(runDirectory, '/', 'genome_Statistics.txt');
writeStatistics(obj.regionsInfoDic, obj.segmentsInfoDic, statisticsFilePath);
%
compareFilePath_focalRegions = strcat(runDirectory, '/', 'genome_focalRegions_comparsionList.txt');
bedFilePath1_focalRegions = strcat(runDirectory, '/', 'genome_focalRegions_narrowPeak.bed');
bedFilePath2_focalRegions = strcat(runDirectory, '/', 'genome_focalRegions_RGB.bed');
focalCNAsForComparsion(obj, obj.regionsInfoDic, bedFilePath1_focalRegions, bedFilePath2_focalRegions, compareFilePath_focalRegions);
%
compareFilePath_IBs = strcat(runDirectory, '/', 'genome_IBs_comparsionList.txt');
bedFilePath1_IBs = strcat(runDirectory, '/', 'genome_IBs_narrowPeak.bed');
bedFilePath2_IBs = strcat(runDirectory, '/', 'genome_IBs_RGB.bed');
IBsForComparsion(obj, obj.segmentsInfoDic, bedFilePath1_IBs, bedFilePath2_IBs, compareFilePath_IBs);


end
%%%
%%%
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------------------------------------- Subroutines -------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%++ Writing statistics  ++%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parsedData] = writeStatistics(regionsDic, segmentsDic, fileName)

chromosomes = cell2mat(keys(regionsDic));

%Grey-Scale File
fid = fopen(fileName,'w');
fprintf(fid, 'Chromosome\t#IBs\t#CNVs\t#Amplifications\t#Deletions \n');

genomeAlterations = 0;
genomeAmplifications = 0;
genomeDeletions = 0;
genomeSegments = 0;

for k=1:length(chromosomes)
    chrNumber      = chromosomes(k);
    regionsPerChr  = regionsDic(k);
    [noRegions,~]  = size(regionsPerChr);
    segmentsPerChr = segmentsDic(k);
    [noSegments,~] = size(segmentsPerChr);
    
    if(noRegions > 0)
        regionsStart = regionsPerChr(:,2);
        regionsEnd   = regionsPerChr(:,3);
        copyNumber   = regionsPerChr(:,5);
        alteration   = regionsPerChr(:,6);

        %% ------- Chromosome Level -------%%
        noAlterarions    = length(alteration);
        noAmplifications = sum(alteration > 0);
        noDeletions      = sum(alteration < 0);

        %% --------- Genome Level ---------%%
        genomeAlterations = genomeAlterations + noAlterarions;
        genomeAmplifications = genomeAmplifications + noAmplifications;
        genomeDeletions = genomeDeletions + noDeletions;
        genomeSegments = genomeSegments + noSegments;

        if(chrNumber==23)
            chrName = 'chrX';
        else
            chrName = strcat('chr',int2str(chrNumber));
        end
        fprintf(fid, '%s\t%i\t%i\t%i\t%i \n', chrName, noSegments, noAlterarions, noAmplifications, noDeletions); 
    else
        genomeSegments = genomeSegments + 1;
        if(chrNumber==23)
            chrName = 'chrX';
        else
            chrName = strcat('chr',int2str(chrNumber));
        end	
        fprintf(fid, '%s\t%i\t%i\t%i\t%i \n', chrName, 1, 0, 0, 0);
    end
end

fprintf(fid,'%s\t%i\t%i\t%i\t%i \n', 'genome', genomeSegments, genomeAlterations, genomeAmplifications, genomeDeletions); 
fclose(fid);
%sprintf('%s\t%i\t%i\t%i\t%i \n', 'genome', genomeSegments, genomeAlterations, genomeAmplifications, genomeDeletions 
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%+++ Writing Bed Files contains all focal alterations +++%%%%%%%%%%%%%%%%%%%%%%%
function [parsedData] = focalCNAsForComparsion(obj, regionsDic, fileName1, fileName2, fileName3)

chromosomes = cell2mat(keys(regionsDic));
binSize = obj.binSize;

%NarrowPeak File
fid = fopen(fileName1,'w');
fprintf(fid, 'track type=narrowPeak visibility=3 description="CNAtra focal amplifications/deletions" \n');

%RGB File
%fid2 = fopen(fileName2,'w');
%fprintf(fid2, 'track name="RGBScale" visibility=2 description="CNAtra focal amplifications/deletions" itemRgb="On" \n');

%CNAtra Output to compare
fid3 = fopen(fileName3,'w');
%fprintf(fid3, 'Chromosome Region-start Region-end Rank Copy-number focal-Amp/Del Class ID \n');


ampNumber = 0;
delNumber = 0;
parsedData1 = [];
for k=1:length(chromosomes)
    chrNumber     = chromosomes(k);
    regionsPerChr = regionsDic(k);
    [noRegions,~] = size(regionsPerChr);
    
    if(noRegions > 0)
        regionsStart = regionsPerChr(:,2);
        regionsEnd   = regionsPerChr(:,3);
        copyNumber   = regionsPerChr(:,5);
        alteration   = regionsPerChr(:,6);

        for i=1:noRegions  
            
            %%Bed File
            if(chrNumber==23)
                chrName = 'chrX';
            else
                chrName = strcat('chr',int2str(chrNumber));
            end
            aRegionStart = ((regionsStart(i)-1)*binSize);
            aRegionEnd   = (regionsEnd(i)*binSize)-1;
            if(sign(alteration(i))==1)
                ampNumber = ampNumber+1;
                altName = strcat('amp',int2str(ampNumber));
            else
                delNumber = delNumber+1;
                altName = strcat('del',int2str(delNumber));
            end
	    %%%%%%
	    %Normalized copy number
	    %assume copy-number has a range [0,10]
	    CN = copyNumber(i);
	    nCN = min(CN/10,1);

        %Grey-score for IB.
	    score = max(min(floor(nCN*1000),999),1);

	    %%%%%%
            switch alteration(i)
                case -5
                    rValue = 0;
                    gValue = 0;
                    bValue = 50 - min(copyNumber(i)*10, 40);
                    regionClass = 'Deletion';
                case 5
                    rValue = 200+ min(copyNumber(i)*10, 40);
                    gValue = 0;
                    bValue = 0;                
                    regionClass = 'Amplification';
                otherwise
					rValue = 0;
                    gValue = 100;
                    bValue = 0;
                    regionClass = 'Neutral';   
            end
	    %
        if(alteration(i) > 0)
			strand = '+';
	    else
            strand = '-';
	    end

	    %%%%%%%%
            score = score + min(copyNumber(i)*20, 90);
            fprintf(fid, '%s\t%i\t%i\t%s\t%i\t%s\t%2.4f\t%i\t-1\t-1 \n', chrName, aRegionStart, aRegionEnd, regionClass, score, strand, CN, i);
            %fprintf(fid2, '%s\t%i\t%i\t%s\t%2.4f\t%s\t%i\t%i\t%i,%i,%i\t%i \n', chrName, aRegionStart, aRegionEnd, regionClass, CN, strand, aRegionStart, aRegionEnd, rValue, gValue, bValue, i);
            fprintf(fid3, '%s\t%i\t%i\t%i\t%2.4f\t%s\t%s\t%i\n', chrName, aRegionStart, aRegionEnd, alteration(i), copyNumber(i), altName, regionClass, i);
        end 
    end
end

fclose(fid);
%fclose(fid2);
fclose(fid3);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%+++++++++ Writing Bed Files contains all IBs ++++++++++%%%%%%%%%%%%%%%%%%%
function [parsedData] = IBsForComparsion(obj, regionsDic, fileName1, fileName2, fileName3)

chromosomes = cell2mat(keys(regionsDic));
binSize = obj.binSize;

%Grey-Scale File
fid = fopen(fileName1,'w');
fprintf(fid, 'track type=narrowPeak visibility=3 description="CNAtra IBs" \n');

%RGB File
%fid2 = fopen(fileName2,'w');
%fprintf(fid2, 'track name="RGBScale" visibility=2 description="CNAtra IBs" itemRgb="On" \n');

%CNAtra Output to compare
fid3 = fopen(fileName3,'w');
%fprintf(fid3, 'Chromosome IB-start IB-end Copy-number ID \n');


IBsNumber = 0;
parsedData1 = [];
for k=1:length(chromosomes)
    chrNumber     = chromosomes(k);
    regionsPerChr = regionsDic(k);
    [noRegions,~] = size(regionsPerChr);
    
    if(noRegions > 0)
        regionsStart = regionsPerChr(:,2);
        regionsEnd   = regionsPerChr(:,3);
        copyNumber   = regionsPerChr(:,5);

        for i=1:noRegions  
            
            %%Bed File
            if(chrNumber==23)
                chrName = 'chrX';
            else
                chrName = strcat('chr',int2str(chrNumber));
            end
            aRegionStart = ((regionsStart(i)-1)*binSize);
            aRegionEnd   = (regionsEnd(i)*binSize)-1;
	    %
            IBsNumber = IBsNumber+1;
            altName = strcat('IB',int2str(IBsNumber));
	    %
	    %Normalized copy number
	    %assume copy-number has a range [0,10]
	    CN = copyNumber(i);
	    nCN = min(CN/10,1);

        %Grey-score for IB.
	    score = max(1,min(floor(nCN*1000),999));

	    %RGB score for IB.
	    if (CN < 1.5)
	            rValue = 0;
        	    bValue = floor(nCN*250);
        	    gValue = 0;
	    elseif(CN < 2.5)
	            rValue = 0;
        	    bValue = 0;
        	    gValue = floor(nCN*250);
	    else
	            rValue = floor(nCN*250);
        	    bValue = 0;
        	    gValue = 0;
	    end

            fprintf(fid, '%s\t%i\t%i\t%s\t%i\t%s\t%2.4f\t%i\t-1\t-1 \n', chrName, aRegionStart, aRegionEnd, altName, score, '.', CN, i);
            %fprintf(fid2, '%s\t%i\t%i\t%s\t%2.4f\t%s\t%i\t%i\t%i,%i,%i\t%i \n', chrName, aRegionStart, aRegionEnd, altName, CN, '.', aRegionStart, aRegionEnd, rValue, gValue, bValue, i);
            fprintf(fid3, '%s\t%i\t%i\t%2.4f\t%s\t%i\n', chrName, aRegionStart, aRegionEnd, copyNumber(i), altName, i);
        end 
    end
end
fclose(fid);
%fclose(fid2);
fclose(fid3);

end