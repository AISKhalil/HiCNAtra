function obj = RDnormalization(obj)
%Normalization the RD signal for the bias sources.

%%-------------------- RD-signal --------------------%%
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

RDdata = [];
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    RDdata   = [RDdata; obj.chrFDictionary(chrIndex)];
end


%%---------------- Bias-sources flatten -------------%%
[e, g, m, chrStart] = flattenFilterBias(obj);
obj.chrStartIndex = chrStart;


%%--------------------Normalization -----------------%%
if (obj.CNAtraRdNormalization == 1)
	%%%% CNAtra-normalization %%%%
	[RDdataNormalized] = CNAtraNormalization(RDdata, m, g, e, obj.binSize);
	for i  = 1:1:noChrs 
	    chrIndex = chromosomes(i);
	    obj.chrFNDictionary(chrIndex) = RDdataNormalized(chrStart(i):chrStart(i+1)-1);
	end
else
	%%%% GLM-normalization %%%%
	[RDdataNormalized, normalizationParameters] = GLMDataNormalization(RDdata, m, g, e);
	for i  = 1:1:noChrs 
	    chrIndex = chromosomes(i);
	    obj.chrFNDictionary(chrIndex) = RDdataNormalized(chrStart(i):chrStart(i+1)-1)*median(RDdata);%To have the same range of numbers.
	end
end
%%%


end
%%%
%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------------------------- Sub-routines ------------------------------%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e1, g1, m1, chrStart1] = flattenFilterBias(obj)
%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
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
chrNames = obj.chrNames;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chrEffectiveLengthTracks = obj.chrEffectiveLengthTracks;
chrMappabilityTracks = obj.chrMappabilityTracks;
chrGCTracks = obj.chrGCTracks;
chrFIndex = obj.chrFIndex;



%%%%%%%%
e1 = [];
g1 = [];
m1 = [];
%%
chrStart1 = [];
firstBinPerChr1= 1;
chrStart1 = [chrStart1; firstBinPerChr1];



for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(chrIndex));

    %
    whiteRegions = chrFIndex(chrIndex);
	firstBinPerChr1= firstBinPerChr1 + length(whiteRegions);
	chrStart1 = [chrStart1; firstBinPerChr1];
    %
    e1ChrData = chrEffectiveLengthTracks(chrIndex);
    m1ChrData = chrMappabilityTracks(chrIndex);
    g1ChrData = chrGCTracks(chrIndex);
	%
	e1 = [e1; e1ChrData(whiteRegions)];
	m1 = [m1; m1ChrData(whiteRegions)];
	g1 = [g1; g1ChrData(whiteRegions)];
end
%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RDnormalized, iParam] = GLMDataNormalization(RDdata, m, g, e)
%%%%%%%%
	m1 = (m-mean(m))/std(m);
	e1 = (e-mean(e))/std(e);
	g1 = (g-mean(g))/std(g);
	%
	m1 = m1 - min(m1) + 1;
	e1 = e1 - min(e1) + 1;
	g1 = g1 - min(g1) + 1;
	%
	features = [e1, m1, g1];	
	%
	distr = 'poisson';
	link = 'log';
	X = log(features);
	y = RDdata;

	iParam = glmfit(X,y, distr,'link', link);
	yhat   = glmval(iParam,X,link);
	RDnormalized = y ./ yhat;
end
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RDnormalized] = CNAtraNormalization(RD, m, g, e, binSize)
%%%%%%%%

%Method#1: parallel-normalization
[coeff1] = intervalNormalization(RD, m);
[coeff2] = intervalNormalization(RD, g);
[coeff3] = intervalNormalization(RD, e/binSize);
RDnormalized = ((RD.*coeff1).*coeff2).*coeff3;

end
%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correctionRatioVec] = intervalNormalization(RD, gc)

%
correctionRatioVec = ones(length(RD),1);
rdGlobal = nanmedian(RD);
%
minGC = min(gc);
maxGC = max(gc);
rangeBoundaries = minGC:0.01:maxGC;
noIntervals = length(rangeBoundaries)-1;
%
for j = 1:noIntervals
	intervalIndex   = (rangeBoundaries(j) <= gc & gc < rangeBoundaries(j+1));
	if(sum(intervalIndex)>1)
	   intervalData    = RD(intervalIndex);
	   intervalLength  = length(intervalData);
	   intervalRDMedian  = nanmedian(intervalData);
	   correctionRatio = rdGlobal/intervalRDMedian;
	   if(intervalLength == 0 || intervalRDMedian == 0)
		intervalDataMod = 1;
	   else
		intervalDataMod = correctionRatio;
	   end
	   correctionRatioVec(intervalIndex) = intervalDataMod;
	end 
end 
%%%

end
