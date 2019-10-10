function obj = mappabilityFeatureCalculater(obj)
%Computing the mappability feature per bin using the mappability window.


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
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
%
mappabilityWindow = obj.mappabilityWindow;
rsites = obj.rsites;
%
binSize = obj.binSize;


%%%
for i  = 1:1:noChrs 
	chrIndex = chromosomes(i);
	chromosome = cell2mat(chrNames(i));
	chrLengthBps = chrLengths(chrIndex);
	chrLengthBins = ceil(chrLengthBps/binSize);
	chrRsites = rsites(chrIndex);
	%%
	mappabilityWindowLeft = max([chrRsites - mappabilityWindow, chrRsites-1],1);
	mappabilityWindowRight  = min([chrRsites, chrRsites + mappabilityWindow - 1],chrLengthBps);
	%%
	mappabilityWindowScores = obj.mappabilityWindowScores(chrIndex);
	%%
	mappabilityBins = zeros(chrLengthBins,1);
	mappabilityBinsStart = 0:binSize:chrLengthBps;
	mappabilityBinsStop  = [binSize-1:binSize:chrLengthBps,chrLengthBps];

	for j = 1:length(mappabilityBins)
	    rsitesPerBin = find(mappabilityBinsStart(j) <= chrRsites & chrRsites < mappabilityBinsStop(j));
	    rsitesFragments = mappabilityWindowScores(rsitesPerBin,:);
	    if(~isempty(rsitesFragments))
		mappabilityBins(j) = mean(rsitesFragments(:));
	    end
	end
	%%%
	obj.chrMappabilityTracks(chrIndex) = mappabilityBins;
end    
%%%s


end
