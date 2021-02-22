function obj = plotGenome(obj, saveResult, nBinSize)  

%Plotting the genome RD signal.
if(strcmp(saveResult,'save')== 1)
    saveResult = 1;
else
    saveResult = 0;
end
%%%
if(saveResult == 1)
    figure;
    set(gcf, 'Position', [1 1 1200 600]);
    set(gcf,'Renderer','painters');
    hold on;
else
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
end
%%%
format long;
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
    mkdir(dir);
end
runDirectory = strcat(dir,'/','figures');
if(exist(runDirectory,'dir') ~= 7)
    mkdir(runDirectory);
end

%% ---------------- Data --------------- %%
chrData = obj.chrDictionary;
segmentsInfoDic = obj.segmentsInfoDic;
CNReference  = obj.CNReference;
binSize = obj.binSize;
%
chrs = obj.targetChrs;
if (chrs == 0) 
    chromosomes = [1:1:length(obj.chrNames)];
else
    chromosomes = chrs;
end
noChrs = length(chromosomes);

%% ------------ Processing --------------%%
chrStart = 1;
globalMax = 0;
IBs = [];
xticksPositions = [];
%%%
for i=1:noChrs
    targetChrIndex = chromosomes(i);
    %%%%%%%%%% Data Reading %%%%%%%%

    chrFNDictionary = obj.chrFNDictionary(targetChrIndex);
    chrFIndex = obj.chrFIndex(targetChrIndex);
    chrLength = length(obj.chrDictionary(targetChrIndex));
    %
    orgData = zeros(chrLength,1);
    orgData(chrFIndex) = chrFNDictionary;
    %
    copyNumber = CNReference;
    normData = orgData*2/copyNumber;
    %
    clippedDataLen = length(normData);

    %%-------------- Segment Data  ---------------%%
    segmentsArray = segmentsInfoDic(targetChrIndex);
    [noSegments,~] = size(segmentsArray);
    segmentsStart = segmentsArray(:,2);
    segmentsEnd   = segmentsArray(:,3);
    segmentsCN     = segmentsArray(:,5);
    %
    segmentsCNVector = nan(clippedDataLen,1);
    for j=1:noSegments-1
    	segmentsCNVector(segmentsStart(j):segmentsEnd(j)) = segmentsCN(j);
    end
    segmentsCNVector(segmentsStart(noSegments):end) = segmentsCN(noSegments);
    %
    %centroData = obj.chrCentroTeloBoundaries(targetChrIndex);
    %segmentsCNVector(centroData(1,1):centroData(1,2)) = 0;

    %
    %%----- Clipped Data -------%%
    %Clip the Read-Depth data to a maximum copyNumber of 10.
    maxCN = max(segmentsCN)+2;
    globalMax = max(maxCN,globalMax);
    clippedData = min(normData, maxCN);


    %%%%%%%%%% Scaling %%%%%%%%%
    n = round(nBinSize/binSize);
    %
    a = clippedData;
    sClippedData = arrayfun(@(k) mean(a(k:k+n-1)),1:n:length(a)-n+1)';
    %
    a = segmentsCNVector;
    sSegmentsCNVector = arrayfun(@(k) mean(a(k:k+n-1)),1:n:length(a)-n+1)';


    %%%%%%%%% Plotting %%%%%%%%%
    plotIndices = chrStart:(chrStart+length(sClippedData)-1);
    chrMid = round(chrStart + length(sClippedData)/2);
    chrStart = chrStart + length(sClippedData);
    IBs = [IBs; sSegmentsCNVector];
    %
    if(mod(i,2)==1)
        plot(plotIndices, sClippedData,'LineStyle','none','Marker','.', 'MarkerSize', 4,'Color', [0.3 0.3 0.3]);
        xticksPositions = [xticksPositions; chrMid];
    else
        plot(plotIndices, sClippedData,'LineStyle','none','Marker','.', 'MarkerSize', 4,'Color', [0.6 0.6 0.6]);
    end
end
%%%

plot(IBs,'k','lineWidth',2);
%%%

h = zeros(1, 1);
h(1) = plot(NaN,NaN,'k','lineWidth',2);
%lgd = legend(h, 'IB copy number');
lgd.FontSize = 8;
lgd.FontWeight = 'bold';
%%%

xticks(xticksPositions)
xticklabels({'1','3','5','7','9','11','13','15','17','19','21','X'})
ylim([-0.2,globalMax]);
xlim([0, length(IBs)+10])
set(gca,'FontSize', 8);
set(gca, 'FontWeight', 'bold');
ylabel('Copy number', 'FontSize', 8,'FontWeight','bold');
xlabel('Chromosome', 'FontSize', 8,'FontWeight','bold');
title('Genome CNV profile', 'FontSize', 10,'FontWeight','bold');
hold off;
%%%
if(saveResult == 1)
    ff = '-f111';
    hh = strcat(runDirectory, '/','GenomeCNVsTrack');
    print(ff,hh,'-dpng');
    savefig(hh);
    close All;
end








end
%%%