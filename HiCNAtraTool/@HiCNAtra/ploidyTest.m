function obj = ploidyTest(obj)
%Plotting the histogram of copy numbers for all IBs and the multimodal distribution. This can be used to validate the assumption of the ploidy level. 


%% ---------------- Data --------------- %%
chrData = obj.chrDictionary ;
chrNames = obj.chrNames;
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
chrCentroTeloBoundaries  = obj.chrCentroTeloBoundaries;
chrBlackListedBoundaries = obj.chrBlackListedBoundaries;
chromosomes = cell2mat(keys(chrData));
noChrs = length(chromosomes);


selectedChrs = chromosomes;
colorCode = ['b';'b';'b'];
CNnormalization = 1;
[colorCodeLength, ~] = size(colorCode);
binSize = 1;%1Kb

%% ---------- Interval --------------%%
figure;
hold on; 
totalIBs = 0;
maxWidth = 0;
for i=1:length(selectedChrs)
	IBs = segmentsInfoDic(i);
	noIBs = length(IBs(:,1));
	totalIBs = totalIBs + noIBs;
	for j=1:noIBs
		IB_width = IBs(j,4)/binSize;
		maxWidth = max(maxWidth,IB_width);
		IB_CN = IBs(j,5);
		IB_CN_Range = int32(IB_CN);
                plot([IB_CN,IB_CN],[0,IB_width],'Color',colorCode(mod(IB_CN_Range,colorCodeLength)+1,:),'LineWidth',1);   
	end
end
xlabel('Estimated CN of the IB');
ylabel('Width of the IB (Kbps)');
xlim([0 8]);
%
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';

%% ------- Multimodal Distribution ----------%%
RDDist = formDistribution(obj.ploidyLevel);
x = 0:0.001:8;
y = RDDist(x);
yNormalized = y * maxWidth/max(y);
plot(x,yNormalized,'k','lineWidth',2);


%%%
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	c = 1/(a1 + a1 + a1 + a1 + a1 + a1);
	RDDist = @(x)(c*(a1*normpdf(x,1,sigma)+a1*normpdf(x,2,sigma)+a1*normpdf(x,3,sigma)+a1*normpdf(x,4,sigma)+a1*normpdf(x,5,sigma)+a1*normpdf(x,6,sigma)));
end

%%%
end
