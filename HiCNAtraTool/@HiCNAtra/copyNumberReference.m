function [] = copyNumberReference(obj)
%Computing the copy number reference (2N) by fitting the input data to the multi-model distribution.

ploidyLevel   = obj.ploidyLevel;
%
chrs = obj.targetChrs;
if (chrs == 0) 
    chromosomes = [1:1:length(obj.chrNames)];
else
    chromosomes = chrs;
end
noChrs = length(chromosomes);

signal = [];
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    signal   = [signal; obj.chrFNDictionary(chrIndex)];
end


%%%% Fitting multi-model distribution %%%%%
RDDist = formDistribution(ploidyLevel);

%%%%%%%%%% Copy-Number Reference %%%%%%%%%%
n = min(round(100000/obj.binSize),1);
data = signal(signal ~= 0);
lowerLimit = prctile(data,1);
upperLimit = prctile(data,99);
a = data(data>lowerLimit & data<upperLimit);
binnedD = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';% the averaged vector
[f,l] = hist(binnedD,1000);

scanInterval = l;
noScanItems = length(scanInterval);
inIntervalSum = zeros(noScanItems,1);
for j =1:noScanItems
	CN = scanInterval(j);%Corresponding RD of CN = 2
	nScanInterval = scanInterval*2/CN;
	correlationData = RDDist(nScanInterval).*f;
	inIntervalSum(j) = sum(correlationData);   
end 
chrCNIndex = find(inIntervalSum == max(inIntervalSum));
CNReference = mean(scanInterval(chrCNIndex));    

%%%%%%%%%%%%%%%%%% RD-signal Update %%%%%%%%%%%%%%%%%%%
obj.genomeNormRD = signal *2/CNReference;
obj.CNReference  = CNReference;


%%%
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RDDist] = formDistribution(ploidyLevel)
%%%%%%%%

sigma = 0.5/3;%standard deviation
a1 = 1;
a2 = 0.5;
a3 = 0.25;
a4 = 0.125;
a5 = 0.0625;


if(strcmp(ploidyLevel,'free'))
	c = 1/(a1 + a1 + a1 + a1 + a1);
	RDDist = @(x)(c*(a1*normpdf(x,2,sigma)+a1*normpdf(x,3,sigma)+a1*normpdf(x,4,sigma)+a1*normpdf(x,5,sigma)+a1*normpdf(x,6,sigma)));
elseif(strcmp(ploidyLevel,'diploid'))
	c = 1/(a1 + a2 + a3 + a4 + a5);
	RDDist = @(x)(c*(a1*normpdf(x,2,sigma)+a2*normpdf(x,3,sigma)+a3*normpdf(x,4,sigma)+a4*normpdf(x,5,sigma)+a5*normpdf(x,6,sigma)));
elseif(strcmp(ploidyLevel,'triploid'))
	c = 1/(a1 + a2 + a3 + a4 + a2);
	RDDist = @(x)(c*(a2*normpdf(x,2,sigma)+a1*normpdf(x,3,sigma)+a2*normpdf(x,4,sigma)+a3*normpdf(x,5,sigma)+a4*normpdf(x,6,sigma)));
elseif(strcmp(ploidyLevel,'tetraploid'))
	c = 1/(a1 + a2 + a3 + a2 + a3);
	RDDist = @(x)(c*(a3*normpdf(x,2,sigma)+a2*normpdf(x,3,sigma)+a1*normpdf(x,4,sigma)+a2*normpdf(x,5,sigma)+a3*normpdf(x,6,sigma)));
else
	disp('Free-model is used for estimating the copy-number reference')	
	c = 1/(a1 + a1 + a1 + a1 + a1);
	RDDist = @(x)(c*(a1*normpdf(x,2,sigma)+a1*normpdf(x,3,sigma)+a1*normpdf(x,4,sigma)+a1*normpdf(x,5,sigma)+a1*normpdf(x,6,sigma)));
end
%%%
end
