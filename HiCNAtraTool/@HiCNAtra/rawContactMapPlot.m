function obj = rawContactMapPlot(obj, saveResult, varargin)  
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
	targetChrIndex = chr1Index + 23*chr2Index;
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
	contactMapFile = inputFilePath;
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
	contactMapFile = inputFilePath;
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
	contactMapFile = inputFilePath;
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
%%%
plot1Indices = plot1Start:plot1End;
plot2Indices = plot2Start:plot2End;
selectedData = chrData(plot1Indices, plot2Indices);
%%%

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
%
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
		hh = strcat(runDirectory, '/raw_', chr1Name);
	elseif(plotMode == 2)
		hh = strcat(runDirectory, '/raw_', chr1Name, '_', chr2Name);	
	elseif(plotMode == 3)
		hh = strcat(runDirectory, '/raw_', chr1Name, '_from_', num2str(varargin{2}), '_To_', num2str(varargin{3}));	
	elseif(plotMode == 4)
		hh = strcat(runDirectory, '/raw_', chr1Name, '_from_', num2str(varargin{3}), '_To_', num2str(varargin{4}), '_', chr2Name, '_from_', num2str(varargin{5}), '_To_', num2str(varargin{6}));	
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% CBREWER - This function produces a colorbrewer table (rgb data) for a 
% given type, name and number of colors of the colorbrewer tables. 
% For more information on 'colorbrewer', please visit
% http://colorbrewer2.org/
% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% 
% [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% INPUT:
%   - ctype: type of color table 'seq' (sequential), 'div' (diverging), 'qual' (qualitative)
%   - cname: name of colortable. It changes depending on ctype.
%   - ncol:  number of color in the table. It changes according to ctype and
%            cname
%   - interp_method: interpolation method (see interp1.m). Default is "cubic" )
% 
% A note on the number of colors: Based on the original data, there is
% only a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce 
% the "extended" colormaps.
%
% Example:  To produce a colortable CT of ncol X 3 entries (RGB) of 
%           sequential type and named 'Blues' with 8 colors:
%                   CT=cbrewer('seq', 'Blues', 8);
%           To use this colortable as colormap, simply call:
%                   colormap(CT)
% 
%           To see the various colormaps available according to their types and
%           names, simply call: cbrewer()
%
%  This product includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011
% ------------------------------
% 18.09.2015  Minor fixes, fixed a bug where the 'spectral' color table did not appear in the preview
% load colorbrewer data
load('colorbrewer.mat')
% initialise the colormap is there are any problems
colormap=[];
if (~exist('interp_method', 'var'))
    interp_method='cubic';
end
% If no arguments
if (~exist('ctype', 'var') | ~exist('cname', 'var') | ~exist('ncol', 'var'))
    disp(' ')
    disp('[colormap] = cbrewer(ctype, cname, ncol [, interp_method])')
    disp(' ')
    disp('INPUT:')
    disp('  - ctype: type of color table *seq* (sequential), *div* (divergent), *qual* (qualitative)')
    disp('  - cname: name of colortable. It changes depending on ctype.')
    disp('  - ncol:  number of color in the table. It changes according to ctype and cname')
    disp('  - interp_method:  interpolation method  (see interp1.m). Default is "cubic" )')
    
    disp(' ')
    disp('Sequential tables:')
    z={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
             'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Spectral'};
    disp(z')     
         
    disp('Divergent tables:')
    z={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
    disp(z')
    
    disp(' ')
    disp('Qualitative tables:')
    %getfield(colorbrewer, 'qual')
    z={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
    disp(z')
    plot_brewer_cmap
    return
end
% Verify that the input is appropriate
ctype_names={'div', 'seq', 'qual'};
if (~ismember(ctype,ctype_names))
    disp('ctype must be either: *div*, *seq* or *qual*')
    colormap=[];
    return
end
if (~isfield(colorbrewer.(ctype),cname))
    disp(['The name of the colortable of type *' ctype '* must be one of the following:'])
    getfield(colorbrewer, ctype)
    colormap=[];
    return
end
if (ncol>length(colorbrewer.(ctype).(cname)))
%     disp(' ')
%     disp('----------------------------------------------------------------------')
%     disp(['The maximum number of colors for table *' cname '* is ' num2str(length(colorbrewer.(ctype).(cname)))])
%     disp(['The new colormap will be extrapolated from these ' num2str(length(colorbrewer.(ctype).(cname))) ' values'])
%     disp('----------------------------------------------------------------------')
%     disp(' ')
    cbrew_init=colorbrewer.(ctype).(cname){length(colorbrewer.(ctype).(cname))};
    colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
    colormap=colormap./255;
    return
end
if (isempty(colorbrewer.(ctype).(cname){ncol}))
    
    while(isempty(colorbrewer.(ctype).(cname){ncol}))
        ncol=ncol+1;
    end        
    disp(' ')
    disp('----------------------------------------------------------------------')
    disp(['The minimum number of colors for table *' cname '* is ' num2str(ncol)])
    disp('This minimum value shall be defined as ncol instead')
    disp('----------------------------------------------------------------------')
    disp(' ')
end
colormap=(colorbrewer.(ctype).(cname){ncol})./255;
end