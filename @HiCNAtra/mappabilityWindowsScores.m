function obj = mappabilityWindowsScores(obj)
%Computing the mappability score for each mappability window.


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


%%----------- Mappability Scores ----------------%%
%% For using Anshul's mappability tracks
readLength = obj.readLength;
if(readLength == 100)
	readLength = readLength+1;
end
%
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
rsites = obj.rsites;
%
mappabilityWindow = obj.mappabilityWindow;
mappabilityFolder = obj.mappabilityFolder;


%%%
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(i));
    chrLengthBps = chrLengths(chrIndex);
    chrRsites = rsites(chrIndex);

    %%%%%%%%%%%%%%%%%%%%%%
    % Anshul_Mappability %

    mappFileName = strcat(mappabilityFolder, chromosome, '.uint8.unique')
    tmp_uMap = fopen(mappFileName,'r');
    uMapdata = fread(tmp_uMap,'*uint8');
    fclose(tmp_uMap);
    mappedBps = (uMapdata > 0 & uMapdata<=readLength);


    %mappabilityWindow is divided into two windows (left and right)
    mappabilityWindowLeft = max([chrRsites - mappabilityWindow, chrRsites-1],1);
    mappabilityWindowRight  = min([chrRsites, chrRsites + mappabilityWindow - 1],chrLengthBps);

    mappabilityWindowScores = zeros(length(mappabilityWindowLeft),2);
    for i = 1:length(mappabilityWindowLeft)
        mappabilityWindowScores(i,1) = mean(mappedBps(mappabilityWindowLeft(i,1):mappabilityWindowLeft(i,2)));
        mappabilityWindowScores(i,2) = mean(mappedBps(mappabilityWindowRight(i,1):mappabilityWindowRight(i,2)));
    end
    %%%

    %%------------------ output dictionary -------------------%%
    obj.mappabilityWindowScores(chrIndex) = mappabilityWindowScores;
end    


end
