 function obj = CNVsTrackPlot(obj, saveResult, varargin)  
%plotting the CNAtra result of a chromosome, iso-copy numeric block, or genomic region.  


%% ------------- Plot Mode -------------- %%
nargin = length(varargin);
switch nargin,
    case 0,
        error('Error: No enough arguments.');    
    case 1,
		% Mode:1, plot whole chromosome. 
		plotMode = 1;
		chrNo = varargin{1};
    case 2,
                % Mode:2, plot a segment.
		plotMode = 2;
		chrNo = varargin{1};
		segmentNo = varargin{2};
	case 3,
		% Mode:3, plot an area of a chromosome (in Kb).
		plotMode = 3;
		chrNo = varargin{1};
        areaStart = floor(varargin{2}/obj.binSize);
		areaEnd   = ceil(varargin{3}/obj.binSize);
    otherwise
        error('Error: Too many arguments.')
end


%%%%%%
format long;
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
	mkdir(dir);
end
runDirectory = strcat(dir,'/','figures');
if(exist(runDirectory,'dir') ~= 7)
	mkdir(runDirectory);
end
%%%
if(strcmp(saveResult,'save')== 1)
	saveResult = 1;
else
	saveResult = 0;
end


%% ---------------- Data --------------- %%
chrNames = obj.chrNames;
CNReference  = obj.CNReference;
chrCNVsTracks = obj.chrCNVsTracks;
chrCentroTeloBoundaries  = obj.chrCentroTeloBoundaries;
%
binSize = obj.binSize;
%
if chrNo == 0
	% Plot all chromosomes
    chromosomes = cell2mat(keys(chrNames));
else
	% Plot one chromosome
    chromosomes = chrNo;
end
noChrs = length(chromosomes);



%% ----------- Processing --------------%%
for i=1:noChrs
	targetChrIndex = chromosomes(i);
	targetChr = cell2mat(chrNames(targetChrIndex));

    %%%%%%%%%%%%%%%%%%%%% Data Reading %%%%%%%%%%%%%%%%%%%
	chrFNDictionary = obj.chrFNDictionary(targetChrIndex);
	chrFIndex = obj.chrFIndex(targetChrIndex);
	chrLength = length(obj.chrDictionary(targetChrIndex));
	%
	orgData = zeros(chrLength,1);
	orgData(chrFIndex) = chrFNDictionary;
	%
	copyNumber = CNReference;
	normData = orgData*2/copyNumber;


	%%-------- Segment Data  ----------%%
	segmentsCNVector = chrCNVsTracks(targetChrIndex);


	%%---------- Clipped Data ----------%%
	%Clip the Read-Depth data to a maximum copyNumber of 10.
	maxCN = max(10, max(segmentsCNVector)+ 2);

	%
	clippedData = min(normData, maxCN);
	clippedDataLen = length(clippedData);
	

	%%------------- Filtered-Regions ---------%%
	regionsBlackGapData = nan(clippedDataLen,1);
	chrBinsNumbers = 1:clippedDataLen;
	blackBinsNumbers = setdiff(chrBinsNumbers, chrFIndex);
	regionsBlackGapData(blackBinsNumbers) = clippedData(blackBinsNumbers);


	%%------ Telomeres/Centromeres ------%%
	centroTeloLocations = chrCentroTeloBoundaries(targetChrIndex);


	%%%--------- Area Selecting for plot ---------%%%
	if(plotMode == 1)
	plotStart = 1;
	plotEnd   = clippedDataLen;
	elseif(plotMode == 2)
	plotStart = segmentsStart(segmentNo);
	plotEnd   = segmentsEnd(segmentNo);
	elseif(plotMode == 3)
	plotStart = max(1, areaStart);
	plotEnd   = min(areaEnd, clippedDataLen);
	end
	plotIndices = plotStart:plotEnd;
	


    %%%------ Ploting Copy-Number regions -------%%%
    if(saveResult == 1)
		figure(targetChrIndex);
		set(gcf, 'Position', [1 1 1200 600]);
		set(gcf,'Renderer','painters');
    else
    	figure;
		set(gcf, 'Position', get(0, 'Screensize'));
    end
    
	
	hold on;
    %Data
	plot(plotIndices, clippedData(plotIndices),'LineStyle','none','Marker','.','MarkerSize',4, 'Color', [0.4 0.4 0.4]);
	%Gap data
	plot(plotIndices, regionsBlackGapData(plotIndices),'k.','MarkerSize',4);
	% Segments CN vector
	plot(plotIndices, segmentsCNVector(plotIndices),'k','lineWidth',3);
	% Centromeres & Telomeres
	centroTelo = centroTeloLocations(:);
    
	legend('Read-Depth', 'Filtered-bins', 'Segments Copy-Number');   
	for j =1:length(centroTelo)
	linePt = centroTelo(j); 
	   if (linePt >= plotStart && linePt <= plotEnd);
		  plot([linePt,linePt],[0,5],'Color','g','LineWidth',3);   
		  legend('Read-Depth', 'Filtered-bins', 'Segments Copy-Number', 'Centromeres-Telomeres');
	   end
	end    
	
	plot(plotIndices, segmentsCNVector(plotIndices),'r','lineWidth',1.5);
	ylim([-0.5,maxCN]);
	ylabel('Copy number');
	xlabel('bin Number');
	title(targetChr);
	hold off;


	%%%------ saving output-figure --------%%%
	targetChr = cell2mat(chrNames(targetChrIndex));
	%
    if(saveResult == 1)
        ff = strcat('-f',num2str(targetChrIndex));
        if(plotMode == 1)
            hh = strcat(runDirectory, '/', targetChr, '_CNVsTrack');
        else
            hh = strcat(runDirectory, '/', targetChr, '_CNVsTrack', '_fromBin_', num2str(plotStart), '_ToBin_', num2str(plotEnd));	
        end
        print(ff,hh,'-dpng');
        savefig(hh);
        close All;
    end


end
