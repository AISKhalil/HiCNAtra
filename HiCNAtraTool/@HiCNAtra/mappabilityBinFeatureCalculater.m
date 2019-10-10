function obj = mappabilityBinFeatureCalculater(obj)
%Computing the mappability feature per bin.


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


%%----------- Mappability feature (window-based) -----------%%
readLength = obj.readLength;
%% For using Anshul's mappability tracks
if(readLength == 100)
	readLength = readLength+1;
end

mappabilityFolder = obj.mappabilityFolder;
%
chrNames = obj.chrNames;
%
binSize = obj.binSize;


%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(i));

    %%%%%%%%%%%%%%%%%%%%%%
    % Anshul_Mappability %
    mappFileName = strcat(mappabilityFolder, chromosome, '.uint8.unique')
    tmp_uMap = fopen(mappFileName,'r');
    uMapdata = fread(tmp_uMap,'*uint8');
    fclose(tmp_uMap);
    mappedBps = (uMapdata > 0 & uMapdata<=readLength);

    n = binSize;
    a = mappedBps;
    mappabilityBins = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
    %%%
    obj.chrMappabilityTracks(chrIndex) = mappabilityBins;
end    
%%%


end
