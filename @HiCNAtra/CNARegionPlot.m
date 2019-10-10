function obj = CNARegionPlot(obj, saveResult, varargin)  
%plotting the CNAtra result of a CNV region using its order in the chromosome.  

%% ------------- Plot Mode -------------- %%
nargin = length(varargin);
switch nargin,
    case 0,
        error('Error: No enough arguments.');    
    case 1,
        error('Error: No enough arguments.');    		
    case 2,
		chrNo = varargin{1};
		regionNo = varargin{2};
    otherwise
        error('Error: Too many arguments.')
end

if(strcmp(saveResult,'save')== 1)
	saveResult = 1;
else
	saveResult = 0;
end

%%%%%%%%%%%%
format long;
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
	mkdir(dir);
end

runDirectory = strcat(dir,'/','figures');
if(exist(runDirectory,'dir') ~= 7)
	mkdir(runDirectory);
end
%
if(saveResult == 1)
	figure(targetChrIndex);
	set(gcf, 'Position', [1 1 800 400]);
	set(gcf,'Renderer','painters');
else
	figure;
	set(gcf, 'Position', get(0, 'Screensize'));
end



%% ---------------- Data --------------- %%
chrNames = obj.chrNames;
segmentsInfoDic = obj.segmentsInfoDic;
regionsInfoDic  = obj.regionsInfoDic;
CNReference  = obj.CNReference;
chrCentroTeloBoundaries  = obj.chrCentroTeloBoundaries;
binSize = obj.binSize;


%% ----------- Processing --------------%%
targetChrIndex = chrNo;
targetChr = cell2mat(chrNames(targetChrIndex));


%%%%%%%%%%%%%% Data Reading %%%%%%%%%%%%%%
chrFNDictionary = obj.chrFNDictionary(targetChrIndex);
chrFIndex = obj.chrFIndex(targetChrIndex);
chrLength = length(obj.chrDictionary(targetChrIndex));
%
orgData = zeros(chrLength,1);
orgData(chrFIndex) = chrFNDictionary;
%
copyNumber = CNReference;
normData = orgData*2/copyNumber;
clippedDataLen = length(normData);	


%%------ Segment Data  ------------%%
segmentsArray = segmentsInfoDic(targetChrIndex);
[noSegments,~] = size(segmentsArray);
segmentsStart = segmentsArray(:,2);
segmentsEnd   = segmentsArray(:,3);
segmentsCN     = segmentsArray(:,5);

segmentsCNVector = nan(clippedDataLen,1);
for j=1:noSegments
	segmentsCNVector(segmentsStart(j):segmentsEnd(j)) = segmentsCN(j);
end


%%----------- Clipped Data -----------%%
%Clip the Read-Depth data to a maximum copyNumber of 10.
maxCN = max(10, max(segmentsCN )+ 2);
%
clippedData = min(normData, maxCN);



%%------ Regions Data  ------------%%
regionsArray = regionsInfoDic(targetChrIndex);
[noRegions,~] = size(regionsArray);
regionsStart  = regionsArray(:,2);
regionsEnd    = regionsArray(:,3);
regionsWidth  = regionsArray(:,4);
regionsCN     = regionsArray(:,5);
regionsCategory = regionsArray(:,6);

regionsCNData   = nan(clippedDataLen,1);
regionsCNVector = segmentsCNVector;
regionsCategoryVector = zeros(clippedDataLen,1);
for j=1:noRegions
	regionsCNData(regionsStart(j):regionsEnd(j)) = clippedData(regionsStart(j):regionsEnd(j));
	regionsCNVector(regionsStart(j):regionsEnd(j)) = regionsCN(j);
	regionsCategoryVector(regionsStart(j):regionsEnd(j)) = regionsCategory(j);
end


%%---- BlackListed/Gap Regions ----%%
regionsBlackGapData = nan(clippedDataLen,1);
chrBinsNumbers = 1:clippedDataLen;
blackBinsNumbers = setdiff(chrBinsNumbers, chrFIndex);
regionsBlackGapData(blackBinsNumbers) = clippedData(blackBinsNumbers);


%%----- Telomeres/Centromeres -----%%
centroTeloLocations = chrCentroTeloBoundaries(targetChrIndex);



%%%%%%%%%%%%%%%%%%%%%%%%%% Area Selecting for plot %%%%%%%%%%%%%%%%%%%%%%%%%%
regionStart = regionsStart(regionNo);
regionEnd   = regionsEnd(regionNo);
regionWidth = regionsWidth(regionNo);

%
maxRegionCN = min(max(clippedData(regionStart:regionEnd)),maxCN);
minRegionCN = max(min(clippedData(regionStart:regionEnd)),0);
areaStart = regionStart-20*regionWidth;
areaEnd   = regionEnd + 20*regionWidth; 
plotStart = max(1, areaStart);
plotEnd   = min(areaEnd, clippedDataLen);

plotIndices = plotStart:plotEnd;

%%%%%%%%%%%%%%%%%%%%%% Ploting Copy-Number regions %%%%%%%%%%%%%%%%%%%%%%%%%%

hold on;
%Data
plot(plotIndices, clippedData(plotIndices),'LineStyle','none','Marker','.','Color', [0.4 0.4 0.4]);
%Amplification and deletion data
plot(plotIndices, regionsCNData(plotIndices),'b.','MarkerSize',6);
%BlackListed/gapRegions
plot(plotIndices, regionsBlackGapData(plotIndices),'k.','MarkerSize',4);
% Segments CN vector
plot(plotIndices, segmentsCNVector(plotIndices),'k','lineWidth',4);
%Amplification and deletion CN vector
plot(plotIndices, regionsCNVector(plotIndices),'r','lineWidth',1.5);
% Centromeres & Telomeres
centroTelo = centroTeloLocations(:);


legend('Read-Depth','Alteration bins', 'BlackListed/Gap/Low Mapp bins', 'Segments Copy-Number', 'Alteration Copy-Number');
for j =1:length(centroTelo)
   linePt = centroTelo(j); 
   if (linePt >= plotStart && linePt <= plotEnd);
          plot([linePt,linePt],[0,5],'Color','g','LineWidth',3);   
		  legend('Read-Depth','Alteration bins', 'BlackListed/Gap/Low-Mapp bins', 'Segments Copy-Number', 'Alteration Copy-Number','Centromeres-Telomeres');
   end
end    
  

ylim([-0.5,maxCN]);
ylabel('Copy number');
xlabel('bin Number');
title(targetChr);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Saving Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targetChr = cell2mat(chrNames(targetChrIndex));

if(saveResult == 1)
	ff = strcat('-f',num2str(targetChrIndex));
	if(plotMode == 1)
		hh = strcat(runDirectory, '/', targetChr);
	else
		hh = strcat(runDirectory, '/', targetChr,'_region_',int2str(regionNo), '_fromBin_', num2str(plotStart), '_toBin_', num2str(plotEnd));	
	end
	print(ff,hh,'-dpng');
	savefig(hh);
	close All;
end




