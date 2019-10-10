function obj = normContactMapPlot(obj, saveResult, varargin)  
% plotting the raw intra-chromosome/inter-chromosome interaction matrices.  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------- Plot Mode ----------------------- %%
nargin = length(varargin);
switch nargin,
    case 0,
        error('Error: No enough arguments.');    
    case 1,
		% Mode:1, plot intra-chromosome contact map. 
		plotMode = 1;
		chr1No = varargin{1};
		chr2No = varargin{1};
    case 2,
        % Mode:2, plot inter-chromosome contact map
		plotMode = 2;
		chr1No = varargin{1};
		chr2No = varargin{2};
	case 3,
		% Mode:3, plot an area of intra-chromosome contact map.
		plotMode = 3;
		chr1No = varargin{1};
		chr2No = varargin{1};
        area1Start = floor(varargin{2}/obj.contactMapBinSize);
		area1End   = ceil(varargin{3}/obj.contactMapBinSize);
        area2Start = floor(varargin{2}/obj.contactMapBinSize);
		area2End   = ceil(varargin{3}/obj.contactMapBinSize);
	case 6,
		% Mode:4, plot an area of intra-chromosome contact map.
		plotMode = 4;
		chr1No = varargin{1};
		chr2No = varargin{2};
        area1Start = floor(varargin{3}/obj.contactMapBinSize);
		area1End   = ceil(varargin{4}/obj.contactMapBinSize);
        area2Start = floor(varargin{5}/obj.contactMapBinSize);
		area2End   = ceil(varargin{6}/obj.contactMapBinSize);		
    otherwise
        error('Error: Number of arguments is not-valid.');
end
%%%

if(chr1No ~= 23)
	chr1Name = strcat('chr',int2str(chr1No));
else
	chr1Name = 'chrX';
end
if(chr2No ~= 23)
	chr2Name = strcat('chr',int2str(chr2No));
else
	chr2Name = 'chrX';
end
%%%

chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
binSize  = obj.contactMapBinSize;
%
chr1LengthBps   = chrLengths(chr1No);
chr1Length      = ceil(chr1LengthBps/binSize);
%
chr2LengthBps   = chrLengths(chr2No);
chr2Length      = ceil(chr2LengthBps/binSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(saveResult,'save')== 1)
	saveResult = 1;
	targetChrIndex = chr1No + 23*chr2No;
	figure(targetChrIndex);
	set(gcf, 'Position', [1 1 round(1000*chr2Length/chr1Length)+100 1000]);
	set(gcf,'Renderer','painters');
else
	saveResult = 0;
	figure;
	set(gcf, 'Position', [1 1 round(1000*chr2Length/chr1Length)+100 1000]);
	set(gcf,'Renderer','painters');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------- Data ------------------------ %%



%%----------- Input-matrices ------------%%
if(chr2No == chr1No)
	[inputFilePath, outputFilePath] = findInOutPaths(obj, chr1No, chr2No);
	contactMapFile = outputFilePath;
	%
	%%% Upper-Half
	sparsedData = dlmread(contactMapFile, '\t');
	row = sparsedData(:,1);
	col = sparsedData(:,2);
	v   = sparsedData(:,3);
	chrHalfData = zeros(chr1Length, chr2Length);
	lin_idcs = sub2ind(size(chrHalfData), row, col);
	chrHalfData(lin_idcs) = v;
	%
	%%% Lower-Half
	[n,m] = size(chrHalfData);
	chrData = chrHalfData' + chrHalfData;
	chrData(1:n+1:end) = diag(chrHalfData);
	%%%
elseif(chr2No > chr1No)
	[inputFilePath, outputFilePath] = findInOutPaths(obj, chr1No, chr2No);
	contactMapFile = outputFilePath;
	%
	sparsedData = dlmread(contactMapFile, '\t');
	row = sparsedData(:,1);
	col = sparsedData(:,2);
	v   = sparsedData(:,3);
	%
	chrData = zeros(chr1Length, chr2Length);
	lin_idcs = sub2ind(size(chrData), row, col);
	chrData(lin_idcs) = v;
	%%%
else
	[inputFilePath, outputFilePath] = findInOutPaths(obj, chr2No, chr1No);
	contactMapFile = outputFilePath;
	% Rotate
	sparsedData = dlmread(contactMapFile, '\t');
	col = sparsedData(:,1);
	row = sparsedData(:,2);
	v   = sparsedData(:,3);
	%
	chrData = zeros(chr1Length, chr2Length);
	lin_idcs = sub2ind(size(chrData), row, col);
	chrData(lin_idcs) = v;
	%%%
end


%%------------ Selected-Data ------------%%
if(plotMode == 1 || plotMode == 2)
	plot1Start = 1;
	plot1End   = chr1Length;
	plot2Start = 1;
	plot2End   = chr2Length;
elseif(plotMode == 3 || plotMode == 6)
	plot1Start = max(1, area1Start);
	plot1End   = min(area1End, chr1Length);
	plot2Start = max(1, area2Start);
	plot2End   = min(area2End, chr1Length);
end
%
plot1Indices = plot1Start:plot1End;
plot2Indices = plot2Start:plot2End;
selectedData = chrData(plot1Indices, plot2Indices);
%

nonzeroData = v(:);
%
if(obj.contactMapBinSize <= 100000)
	plotTh1 = 95;
else
	plotTh1 = 98;	
end
selectedDataThreshold1 = prctile(nonzeroData, plotTh1);
%%%

selectedIndices1 = (selectedData > selectedDataThreshold1);
selectedData(selectedIndices1) = selectedDataThreshold1;
M = selectedData;
selectedData = M - diag(diag(M));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------- Plotting ---------------------- %%
x = plot1Indices;
y = plot2Indices;
c = selectedData;
%
customColorMap = bluewhitered(256);
colormap(customColorMap)
imagesc([plot1Start:plot1End], [plot2Start:plot2End], c);
colorbar;
ylabel(chr1Name,'FontSize',12,'FontWeight','bold');
xlabel(chr2Name,'FontSize',12,'FontWeight','bold');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- Save-Figure -----------------------%%
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
%

if(saveResult == 1)
	ff = strcat('-f',num2str(targetChrIndex));
	if(plotMode == 1)
		hh = strcat(runDirectory, '/norm_', chr1Name);
	elseif(plotMode == 2)
		hh = strcat(runDirectory, '/norm_', chr1Name, '_', chr2Name);	
	elseif(plotMode == 3)
		hh = strcat(runDirectory, '/norm_', chr1Name, '_from_', num2str(varargin{2}), '_To_', num2str(varargin{3}));	
	elseif(plotMode == 4)
		hh = strcat(runDirectory, '/norm_', chr1Name, '_from_', num2str(varargin{3}), '_To_', num2str(varargin{4}), '_', chr2Name, '_from_', num2str(varargin{5}), '_To_', num2str(varargin{6}));	
	end
	print(ff,hh,'-dpng');
	savefig(hh);
	close All;
end

%%%
%%%
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inputFilePath, outputFilePath] = findInOutPaths(obj, chr1Index, chr2Index)

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% Saving Interaction-matrices %%
		format long;
		dir = obj.outputDirectory;
		if(exist(dir,'dir') ~= 7)
			mkdir(dir);
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%---- Chromosome Names ----%%%%
		if(chr1Index ~= 23)
			chr1Name = strcat('chr',int2str(chr1Index));
		else
			chr1Name = 'chrX';
		end
		%
		if(chr2Index ~= 23)
			chr2Name = strcat('chr',int2str(chr2Index));
		else
			chr2Name = 'chrX';
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%------ Input Files -------%%%%
		runDirectoryIn = strcat(dir,'/','rawContactMap_binSize',int2str(obj.contactMapBinSize));
		if(exist(runDirectoryIn,'dir') ~= 7)
			mkdir(runDirectoryIn);
		end
		%
		cisDirectory = strcat(runDirectoryIn,'/','cis_matrics');
		transDirectory = strcat(runDirectoryIn,'/','trans_matrics');
		if(exist(cisDirectory,'dir') ~= 7)
			mkdir(cisDirectory);
		end
		if(exist(transDirectory,'dir') ~= 7)
			mkdir(transDirectory);
		end
		%%%
		if(chr1Index == chr2Index)
			inputFilePath = strcat(cisDirectory, '/', chr1Name, '.txt');
		else
			inputFilePath = strcat(transDirectory, '/', chr1Name, '_', chr2Name, '.txt');

		end
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%------ Output Files ------%%%%
		runDirectoryOut = strcat(dir,'/','correctContactMap_binSize',int2str(obj.contactMapBinSize));
		if(exist(runDirectoryOut,'dir') ~= 7)
			mkdir(runDirectoryOut);
		end
		%
		cisDirectory = strcat(runDirectoryOut,'/','cis_matrics');
		transDirectory = strcat(runDirectoryOut,'/','trans_matrics');
		if(exist(cisDirectory,'dir') ~= 7)
			mkdir(cisDirectory);
		end
		if(exist(transDirectory,'dir') ~= 7)
			mkdir(transDirectory);
		end
		%%%
		if(chr1Index == chr2Index)
			outputFilePath = strcat(cisDirectory, '/', chr1Name, '.txt');
		else
			outputFilePath = strcat(transDirectory, '/', chr1Name, '_', chr2Name, '.txt');

		end		
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmap = bluewhitered(m)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
bottom = [0 0 0.5];
botmiddle = [0 0.5 1];
middle = [1 1 1];
topmiddle = [1 0 0];
top = [0.5 0 0];
% Find middle
lims = get(gca, 'CLim');
% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end
% 
% m = 64;
% new = [bottom; botmiddle; middle; topmiddle; top];
% % x = 1:m;
% 
% oldsteps = linspace(0, 1, 5);
% newsteps = linspace(0, 1, m);
% newmap = zeros(m, 3);
% 
% for i=1:3
%     % Interpolate over RGB spaces of colormap
%     newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
% end
% 
% % set(gcf, 'colormap', newmap), colorbar


end