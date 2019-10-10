function obj = effectiveLengthFeatureCalculater(obj)
%Computing the effective length feature per bin.




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



%%----------- Effective length feature -----------%%
effectiveLengthWindow = obj.effectiveLengthWindow;
%
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
rsites = obj.rsites;
%
binSize = obj.binSize;


for i  = 1:1:noChrs 
	chrIndex = chromosomes(i);
	chromosome = cell2mat(chrNames(i));

	chrLengthBps = chrLengths(chrIndex);
	cut_sites2   = rsites(chrIndex);

	%%effectiveLengthWindow starts from "effectiveLengthWindow before Restriction_sites" to  "effectiveLengthWindow after Restriction_sites"
	effectiveLengthWindowStart = max(cut_sites2 - effectiveLengthWindow, 1);
	effectiveLengthWindowStop  = min(cut_sites2 + effectiveLengthWindow - 1,chrLengthBps);

	effectiveLengthBps  = zeros(chrLengthBps,1);
	for j = 1:length(cut_sites2)
	effectiveLengthBps(effectiveLengthWindowStart(j):effectiveLengthWindowStop(j)) = 1;
	end

	a = effectiveLengthBps;
	n = binSize;
	effectiveLengthBins = arrayfun(@(k) sum(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	%%
	obj.chrEffectiveLengthTracks(chrIndex) = effectiveLengthBins;
end    
%%%



end
