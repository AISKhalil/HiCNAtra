function obj = gcFeatureCalculater(obj)
%Computing the gc feature per bin using the gc-window.


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
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
%
rsites = obj.rsites;
gcWindow = obj.gcWindow;
%
binSize = obj.binSize;


%%%
for i  = 1:1:noChrs 
	chrIndex = chromosomes(i);
	chromosome = cell2mat(chrNames(i));
	chrLengthBps = chrLengths(chrIndex);
	chrRsites = rsites(chrIndex);
	%%
	gcWindsWindowLeft  = max([chrRsites - gcWindow, chrRsites-1],1);
	gcWindsWindowRight = min([chrRsites, chrRsites + gcWindow - 1], chrLengthBps);
	%%
	gcWindsWindowScores = obj.GCWindsWindowScores(chrIndex);
	%%
	chrLengthBins = ceil(chrLengthBps/binSize);
	gcBins = zeros(chrLengthBins,1);
	gcBinsStart = 0:binSize:chrLengthBps;
	gcBinsStop  = [binSize-1:binSize:chrLengthBps,chrLengthBps];
	for j = 1:chrLengthBins
	    rsitesPerBin = find(chrRsites >= gcBinsStart(j) & chrRsites < gcBinsStop(j));
	    rsitesFragments = gcWindsWindowScores(rsitesPerBin,:);
	    rsitesFragments = rsitesFragments(:);
	    gcBins(j) = nanmean(rsitesFragments);
	end
	%%
	obj.chrGCTracks(chrIndex) = gcBins;
end    
%%%


end
