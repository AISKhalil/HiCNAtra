function [correlationValues] = RDvsBiasesCorrelation(obj)
%Normalization the RD signal for the bias sources.


noIntervals = 20;
%%----------------------------- RD-signal ----------------------------%%
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
RDdataNormalized = [];

for i  = 1:1:noChrs 
	chrIndex = chromosomes(i);
	RDdata   = [RDdata; obj.chrFDictionary(chrIndex)];
	RDdataNormalized   = [RDdataNormalized; obj.chrFNDictionary(chrIndex)];
end


%%------------------------ Bias-sources flatten ----------------------%%
[e1, g1, m1, chrStart] = flattenFilterBias(obj);
features = [e1, m1, g1];	


%%---------------------------- Binning -------------------------------%%
n = 1;
%%%%%%%%%
if(n > 1)
	a = RDdata;
	RDdata = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	%
	a = RDdataNormalized;
	RDdataNormalized = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	%
	a = e1;
	e1 = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	%
	a = m1;
	m1 = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	%
	a = g1;
	g1 = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';
	%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------------------- Plotting ------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-normalization  %%
%%%%%%%%%%%%%%%%%%%%%%%%

[RDdataNormalizedEff,  eBarInterval] = getBarValues2(RDdata, e1, noIntervals);
[RDdataNormalizedMapp, mBarInterval] = getBarValues(RDdata, m1, noIntervals);
[RDdataNormalizedGC,   gBarInterval] = getBarValues2(RDdata, g1, noIntervals);


f1 = randi(10000,1,1);
figure(f1);
subplot(2,3,1);
barPlot(RDdataNormalizedEff, eBarInterval, noIntervals, 'Effective length', 'Read counts (mean)', [0,0,1]);
xlim([0 noIntervals+1])
subplot(2,3,2);
barPlot(RDdataNormalizedMapp, mBarInterval, noIntervals, 'Mappability', '', [0,0,1]);
xlim([0 noIntervals+1])
subplot(2,3,3);
barPlot(RDdataNormalizedGC, gBarInterval, noIntervals, 'GC-score', '', [0,0,1]);
xlim([0 noIntervals+1])
	
[preSpearmanCorrelationEff,~]  = corr(RDdata, e1,'Type','Spearman');
[preSpearmanCorrelationMapp,~] = corr(RDdata, m1,'Type','Spearman');
[preSpearmanCorrelationGC,~]   = corr(RDdata, g1,'Type','Spearman');

disp(strcat('Pre-normalization correlation'));
disp(strcat('Effective-length = ', num2str(preSpearmanCorrelationEff)));
disp(strcat('Mappability = ', num2str(preSpearmanCorrelationMapp)));
disp(strcat('GC-contents = ', num2str(preSpearmanCorrelationGC)));



%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-normalization  %%
%%%%%%%%%%%%%%%%%%%%%%%%%



[RDdataNormalizedEff,  eBarInterval] = getBarValues2(RDdataNormalized, e1, noIntervals);
[RDdataNormalizedMapp, mBarInterval] = getBarValues(RDdataNormalized, m1, noIntervals);
[RDdataNormalizedGC,   gBarInterval] = getBarValues2(RDdataNormalized, g1, noIntervals);

subplot(2,3,4);
barPlot(RDdataNormalizedEff, eBarInterval, noIntervals, 'Effective length', 'Read count (mean)', [1,0,0]);
xlim([0 noIntervals+1])
subplot(2,3,5);
barPlot(RDdataNormalizedMapp, mBarInterval, noIntervals, 'Mappability', '', [1,0,0]);
xlim([0 noIntervals+1])
subplot(2,3,6);
barPlot(RDdataNormalizedGC, gBarInterval, noIntervals, 'GC-score', '', [1,0,0]);
xlim([0 noIntervals+1])
	
[postSpearmanCorrelationEff,~]  = corr(RDdataNormalized, e1,'Type','Spearman');
[postSpearmanCorrelationMapp,~] = corr(RDdataNormalized, m1,'Type','Spearman');
[postSpearmanCorrelationGC,~]   = corr(RDdataNormalized, g1,'Type','Spearman');

disp(strcat('Pre-normalization correlation'))
disp(strcat('Effective length = ', num2str(postSpearmanCorrelationEff)));
disp(strcat('Mappability = ', num2str(postSpearmanCorrelationMapp)));
disp(strcat('GC-contents = ', num2str(postSpearmanCorrelationGC)));

%%%%%%%%%%%%%%%%%%%%%%%
%%  output           %%
%%%%%%%%%%%%%%%%%%%%%%%
correlationValues = [preSpearmanCorrelationEff, preSpearmanCorrelationMapp, preSpearmanCorrelationGC, postSpearmanCorrelationEff, postSpearmanCorrelationMapp, postSpearmanCorrelationGC];


%% Figures %%
%%%%%%%%%%%%%
format long;
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
    mkdir(dir);
end
%
runDirectory = strcat(dir,'/','figures');
if(exist(runDirectory,'dir') ~= 7)
    mkdir(runDirectory);
end

ff = strcat('-f',num2str(f1));
hh = strcat(runDirectory, '/Interval_normalization_',int2str(obj.CNAtraRdNormalization));
print(ff,hh,'-dpng');
savefig(hh);
close All;

%%%
end
%%%







%%----------------------------------- Sub-routines ------------------------------%%
function [e1, g1, m1, chrStart1] = flattenFilterBias(obj)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chrEffectiveLengthTracks = obj.chrEffectiveLengthTracks;
chrMappabilityTracks = obj.chrMappabilityTracks;
chrGCTracks = obj.chrGCTracks;
chrFIndex = obj.chrFIndex;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [barValues, eBarInterval] = getBarValues(data, counts, noIntervals)
%%%%%%%%
NotNanIndices = find(~isnan(counts));
data = data(NotNanIndices);
counts = counts(NotNanIndices);

minCount = min(counts);
maxCount = 0.95;
y = linspace(minCount, maxCount, noIntervals+1);

barValues = [];
for i = 1:noIntervals
	intervalStart = y(i);
	intervalStop  = y(i+1);
	%%
	countIndex   = find(counts > intervalStart & counts <= intervalStop);
	countData    = data(countIndex);
	if(length(countData)>0)
		countDataAvg = median(countData);
	else
		countDataAvg = 0;
	end
	barValues = [barValues, countDataAvg];
end

eBarInterval = y;
%%%
end
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [barValues, eBarInterval] = getBarValues2(data, counts, noIntervals)
%%%%%%%%
NotNanIndices = find(~isnan(counts));
data = data(NotNanIndices);
counts = counts(NotNanIndices);

minCount = min(counts);
maxCount = max(counts);
y = linspace(minCount, maxCount, noIntervals+1);

barValues = [];
for i = 1:noIntervals
	intervalStart = y(i);
	intervalStop  = y(i+1);
	%%
	countIndex   = find(counts > intervalStart & counts <= intervalStop);
	countData    = data(countIndex);
	if(length(countData)>0)
		countDataAvg = median(countData);
	else
		countDataAvg = 0;
	end
	barValues = [barValues, countDataAvg];
end

eBarInterval = y;
%%%
end
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = barPlot(barValues, intervalEdges, noIntervals, xlabelName, ylabelName, colorCode)

nBars = noIntervals;
barXticks = {};
for i = 1:nBars
	if(intervalEdges(i) < 1)
		barXticks(i) = {num2str(intervalEdges(i),2)};
	else
		barXticks(i) = {num2str(int32(intervalEdges(i)))};
	end
end

b = bar(barValues,0.5);
b(1).FaceColor = colorCode;
b.LineWidth = 0.01;
%set(gca,'xticklabel',barXticks);
ylim([0 max(barValues)*1.1]);
ax = gca;
ax.FontSize = 8;
ax.FontWeight = 'bold';
set(get(gca,'XLabel'), 'String', xlabelName);
set(get(gca,'YLabel'), 'String', ylabelName); 

end
