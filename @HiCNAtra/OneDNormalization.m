function obj = OneDNormalization(obj)


outputDirectory = obj.outputDirectory;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------- chromosomes for Analysis -------------%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------- Filtering-indices ----------------%%
chrNewFIndex = binFiltering(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------- Input-data -------------------%%
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
readLength = obj.readLength;
binSize = obj.contactMapBinSize;
%%%

genomeF = [];
genomeM = [];
genomeG = [];
genomeE = [];
genomeC = [];
%
genomeF2 = [];
genomeM2 = [];
genomeG2 = [];
genomeE2 = [];
genomeC2 = [];
%
genomeRows = [];
genomeCols = [];
genomeEles = [];
genomeEles2 = [];
genomeFIndex = [];
%%%

cDict = containers.Map({1},{[]});
remove(cDict,1);
eDict = containers.Map({1},{[]});
remove(eDict,1);
gDict = containers.Map({1},{[]});
remove(gDict,1);
mDict = containers.Map({1},{[]});
remove(mDict,1);
preNormArray  = [];
postNormArray = [];
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
for i  = 1:1:noChrs 
	%
	chr1Index = chromosomes(i);
	chr1LengthBps   = chrLengths(chr1Index);
	chr1LengthBins  = ceil(chr1LengthBps/binSize);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		j = i;
		%
		chr2Index = chromosomes(j);
		chr2LengthBps  = chrLengths(chr1Index);
		chr2LengthBins  = ceil(chr2LengthBps/binSize);
				
				
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%----------- Data ----------%%
		[rawIntFreq, outputFilePath] = rawMatrixRead(obj, chr1Index, chr2Index, chr1LengthBins, chr2LengthBins);
		%
		effLen1Track = obj.contactMapEffectiveLengthTracks(chr1Index);
		effLen2Track = obj.contactMapEffectiveLengthTracks(chr2Index);
		%
		GC1Track = obj.contactMapGCTracks(chr1Index);
		GC2Track = obj.contactMapGCTracks(chr2Index);
		%
		mapp1Track = obj.contactMapMappabilityTracks(chr1Index);
		mapp2Track = obj.contactMapMappabilityTracks(chr2Index);
		%
		CNVs1Tracks = obj.contactMapCNVsTracks(chr1Index);
		CNVs2Tracks = obj.contactMapCNVsTracks(chr2Index);
	

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%-------- Filtering ----------%%
		chr1FIndices  = chrNewFIndex(chr1Index);
		chr2FIndices  = chrNewFIndex(chr2Index);
		%
		filteredRawIntMatrix = rawIntFreq(chr1FIndices, chr2FIndices);
		%
		effLen1Track = effLen1Track(chr1FIndices);
		effLen2Track = effLen2Track(chr2FIndices);
		%
		GC1Track = GC1Track(chr1FIndices);
		GC2Track = GC2Track(chr2FIndices);
		%
		mapp1Track = mapp1Track(chr1FIndices);
		mapp2Track = mapp2Track(chr2FIndices);
		%
		CNVs1Tracks = CNVs1Tracks(chr1FIndices);
		CNVs2Tracks = CNVs2Tracks(chr2FIndices);
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%-------- Normalization --------%%		
		[f, m, g, e, c, f2, m2, g2, e2, c2, nRow, nCol, nElement, nElement2, subFilteredIndices] = normalizationFeatures(obj, filteredRawIntMatrix, mapp1Track, GC1Track, effLen1Track, CNVs1Tracks);
		%
		genomeF  = [genomeF; f];
		genomeM  = [genomeM; m];
		genomeG  = [genomeG; g];
		genomeE  = [genomeE; e];
		genomeC  = [genomeC; c];
		%
		genomeF2  = [genomeF2; f2];
		genomeM2  = [genomeM2; m2];
		genomeG2  = [genomeG2; g2];
		genomeE2  = [genomeE2; e2];
		genomeC2  = [genomeC2; c2];
		%
		genomeRows = [genomeRows; nRow];
		genomeCols = [genomeCols; nCol];
		genomeEles = [genomeEles; nElement];
		genomeEles2 = [genomeEles2; nElement2];
		genomeFIndex = [genomeFIndex; subFilteredIndices];

		%%%%%%%%%%%%%
		cDict(i) = c;
		eDict(i) = e;
		gDict(i) = g;
		mDict(i) = m;
		[cCorr,~]  = corr(f, cDict(i),'Type','Spearman');
		[eCorr,~]  = corr(f, eDict(i),'Type','Spearman');
		[gCorr,~]  = corr(f, gDict(i),'Type','Spearman');		
		[mCorr,~]  = corr(f, mDict(i),'Type','Spearman');
		preSpearmanCorrelation = [cCorr, eCorr, gCorr, mCorr]
		preNormArray = [preNormArray; preSpearmanCorrelation];

		

end
%%%



%%%%%%%%%%%%%%%%%%%%%%%% Correction %%%%%%%%%%%%%%%%%%

[cCorr,~]  = corr(genomeF, genomeC,'Type','Spearman');
[eCorr,~]  = corr(genomeF, genomeE,'Type','Spearman');
[gCorr,~]  = corr(genomeF, genomeG,'Type','Spearman');		
[mCorr,~]  = corr(genomeF, genomeM,'Type','Spearman');
preSpearmanCorrelation = [cCorr, eCorr, gCorr, mCorr];
preNormArray = [preNormArray; preSpearmanCorrelation];
%
[genomeFNorm] = genomeGLMDataNormalization(genomeF, genomeF2, genomeM2, genomeG2, genomeE2, genomeC, genomeEles2);
%
[cCorr,~]  = corr(genomeFNorm, genomeC,'Type','Spearman');
[eCorr,~]  = corr(genomeFNorm, genomeE,'Type','Spearman');
[gCorr,~]  = corr(genomeFNorm, genomeG,'Type','Spearman');		
[mCorr,~]  = corr(genomeFNorm, genomeM,'Type','Spearman');
%
postSpearmanCorrelation  = [cCorr, eCorr, gCorr, mCorr];
postNormArray = [postNormArray; postSpearmanCorrelation];







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
startIndex = 1;
startIndex2 = 1;
for i  = 1:1:noChrs 
	%
	chr1Index = chromosomes(i);
	chr1LengthBps   = chrLengths(chr1Index);
	chr1LengthBins  = ceil(chr1LengthBps/binSize);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		j = i;
		%
		chr2Index = chromosomes(j);
		chr2LengthBps  = chrLengths(chr1Index);
		chr2LengthBins  = ceil(chr2LengthBps/binSize);
		%
		chr1FIndices  = chrNewFIndex(chr1Index);
		chr2FIndices  = chrNewFIndex(chr2Index);
		
		%%%%%%%%%%%%%%%%
		matrixIndex = k;
		k = k + 1;
		%
		nRow = genomeRows(matrixIndex);
		nCol = genomeCols(matrixIndex);
		nElement = genomeEles(matrixIndex);
		nElement2 = genomeEles2(matrixIndex);
		%
		if(nElement >= 1)
			%
			fNorm = genomeFNorm(startIndex: startIndex + nElement -1);
			filteredIndices = genomeFIndex(startIndex2: startIndex2 + nElement2 -1);
			startIndex  = startIndex + nElement;
			startIndex2 = startIndex2 + nElement2;

			%%%%%
			[cCorr,~]  = corr(fNorm, cDict(i),'Type','Spearman');
			[eCorr,~]  = corr(fNorm, eDict(i),'Type','Spearman');
			[gCorr,~]  = corr(fNorm, gDict(i),'Type','Spearman');		
			[mCorr,~]  = corr(fNorm, mDict(i),'Type','Spearman');
			postSpearmanCorrelation = [cCorr, eCorr, gCorr, mCorr]
			postNormArray = [postNormArray; postSpearmanCorrelation];
				
			%%%%%
			smallMatrix = reshape(fNorm,[length(filteredIndices),length(filteredIndices)]);
			%%
			outMatrix = zeros(nRow,nCol);
			outMatrix(filteredIndices,filteredIndices) = smallMatrix;
		else
			outMatrix = [];
		end

		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		normIntFreq = zeros(chr1LengthBins, chr2LengthBins);
		normIntFreq(chr1FIndices,chr2FIndices) = outMatrix;
		normIntFreq = triu(normIntFreq);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% Saving Interaction-matrices %%
		[inputFilePath, outputFilePath] = findInOutPaths(obj, chr1Index, chr2Index);
		[row col v] = find(normIntFreq);
		dlmwrite(outputFilePath, [row col v], 'delimiter','\t');

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write correlation-results
postNormArray = [postNormArray(2:end,:); postNormArray(1,:)];
%
dir = obj.outputDirectory;
%
preNormFilePath = strcat(dir, '/', 'contactMapPreNorm_OneD',int2str(obj.contactMapBinSize),'.xls');
postNormFilePath = strcat(dir, '/', 'contactMapPostNorm_OneD',int2str(obj.contactMapBinSize),'.xls');
%


%
A = preNormArray;
fid = fopen(preNormFilePath,'wt');
for ii = 1:size(A,1)
	fprintf(fid,'%g\t',A(ii,:));
	fprintf(fid,'\n');
end
fclose(fid);
%
A = postNormArray;
fid = fopen(postNormFilePath,'wt');
for ii = 1:size(A,1)
	fprintf(fid,'%g\t',A(ii,:));
	fprintf(fid,'\n');
end
fclose(fid);
%




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------------- End "contactMap Normalization" ------------------%%%			














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis Sub-routines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------------------------------------------------------------------%%
%% Sub-routines for normalizing interaction-frequencies with GLM poisson-regression model.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chrData, outputFilePath] = rawMatrixRead(obj, chr1No, chr2No, chr1Length, chr2Length)
		%%%
		[inputFilePath, outputFilePath] = findInOutPaths(obj, chr1No, chr2No);
		contactMapFile = inputFilePath;
		%
		sparsedData = dlmread(contactMapFile, '\t');
		row = sparsedData(:,1);
		col = sparsedData(:,2);
		v   = sparsedData(:,3);
		chrHalfData = zeros(chr1Length, chr2Length);
		lin_idcs = sub2ind(size(chrHalfData), row, col);
		chrHalfData(lin_idcs) = v;
		%
		[n,m] = size(chrHalfData);
		chrData = chrHalfData' + chrHalfData;
		chrData(1:n+1:end) = diag(chrHalfData);		

end




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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chrNewFIndex = binFiltering(obj)
%%%%%%%%

%Filtering the input-chromosome to remove low-mappability, black-listed bins.
chrNewFIndex = containers.Map({1},{[]});
remove(chrNewFIndex,1);

%%--------------- Input ---------------%%
chromosomes = [1:1:length(obj.chrNames)];
noChrs = length(chromosomes);
%
chrNames = obj.chrNames;
chrLengths = obj.chrLengths;
% Bin-sizes
RdBinSize = obj.binSize;
contactMapBinSize = obj.contactMapBinSize;
%
contactMapEffectiveLengthTracks = obj.contactMapEffectiveLengthTracks;%number
contactMapGCTracks = obj.contactMapGCTracks;%fraction
contactMapMappabilityTracks = obj.contactMapMappabilityTracks;%fraction
contactMapCNVsTracks = obj.contactMapCNVsTracks;%copy number
%
minEffectiveLength = obj.minEffectiveLength;
minMappabilityThreshold = obj.minMappabilityThreshold;
minGCThreshold = obj.minGCThreshold;
maximumFalseBinsAllowed = obj.maximumFalseBinsAllowed;
%
chrBlackListedBoundaries = obj.chrBlackListedBoundaries;% dictionary of black-listed regions per chromosomes.
chrCentroTeloBoundaries = obj.chrCentroTeloBoundaries;% dictionary of centromeres/telomeres per chromosomes.
chrGapBoundaries = obj.chrGapBoundaries;% dictionary of gap regions per chromosome.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering-out homogenous deletions.
for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(chrIndex));
    chrLengthBps = chrLengths(chrIndex);
	targetChr1KbLength = ceil(chrLengthBps/1000);
	targetChrInitialLength = ceil(chrLengthBps/RdBinSize);
	targetChrFinalLength   = ceil(chrLengthBps/contactMapBinSize);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% filtering deletions/Centromeres/telomeres/black-listed %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	chrFlags = ones(targetChrInitialLength,1);
	%
	%
	%%%%---- Filtering deletion-regions ----%%%%
	regionsArray = obj.regionsInfoDic(chrIndex);
	[noRegions,~] = size(regionsArray);	
	%
	if(noRegions > 0)	
		regionsStart  = regionsArray(:,2);
		regionsEnd    = regionsArray(:,3);
		regionsCN     = regionsArray(:,5);
		regionsCategory = regionsArray(:,6);
		%
		for j=1:noRegions
			if(regionsCN(j) < 0.5)
				chrFlags(regionsStart(j):regionsEnd(j)) = 0;
			end
		end
	end
	%
	%
	%%%%----------- Centromeres/telomeres -----------%%%%
	centroTeloBoundaries = chrCentroTeloBoundaries(chrIndex);
	chrFlags(centroTeloBoundaries(1,1):centroTeloBoundaries(1,2))= 0;
	chrFlags(centroTeloBoundaries(2,1):centroTeloBoundaries(2,2))= 0;
	chrFlags(centroTeloBoundaries(3,1):centroTeloBoundaries(3,2))= 0;
	%
	%
	%%%%--------------- Gab regions -----------------%%%%
	gapBoundaries = chrGapBoundaries(chrIndex);
	[noGapRegions, ~]  = size(gapBoundaries);
	for j=1:noGapRegions
		chrFlags(gapBoundaries(j,1):gapBoundaries(j,2))= 0;    	
	end	
	%
	%
	%%%%------------- Blacklisted regions -----------%%%%
    blackListedBoundaries = chrBlackListedBoundaries(chrIndex);
    [noBlackRegions,~] = size(blackListedBoundaries);
    for j=1:noBlackRegions
		chrFlags(blackListedBoundaries(j,1):blackListedBoundaries(j,2))= 0;    	
	end		
	%
	%
	%%%%%---------------- 1Kb Bin-size --------------%%%%
	scalingRatio = RdBinSize/1000;
	chrFlags1Kb = repelem(chrFlags, scalingRatio);
	chrFlags1Kb = chrFlags1Kb(1:targetChr1KbLength);
	%
	%
	%%%%%---------- Contact-map Bin-size ------------%%%%
	a = chrFlags1Kb;
	n = contactMapBinSize/1000;
	correctThreshold = arrayfun(@(k) mean(a(k:min(k+n-1,length(a)))),1:n:length(a))';% the averaged vector
	selectedBins = (correctThreshold > maximumFalseBinsAllowed);
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%- filtering low-mappability,low-GC,low effective-length %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	mappabilityTracks = contactMapMappabilityTracks(chrIndex);
	mappabilityBins   = (mappabilityTracks > minMappabilityThreshold & ~isnan(mappabilityTracks));
	selectedBins = bitand(selectedBins, mappabilityBins);
	%
	GCTracks = contactMapGCTracks(chrIndex);
	GCBins   = (GCTracks > minGCThreshold & ~isnan(GCTracks));
	selectedBins = bitand(selectedBins, GCBins);
	%
	effTracks = contactMapEffectiveLengthTracks(chrIndex);
	effBins   = (effTracks > (minEffectiveLength * contactMapBinSize));
	selectedBins = bitand(selectedBins, effBins);
	%
	%
	%
	chrNewFIndex(chrIndex) = find(selectedBins == 1);
end    
%%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, m, g, e, c, fSum, mSum, gSum, eSum, cSum, nRow, nCol, n, nSum, filteredIndices] = normalizationFeatures(obj, inMatrix, m1, g1, e1, c1)
%%%%%%%%
	%%%%  Interaction-matrix  %%%%
	[nRow, nCol] = size(inMatrix);
	%
	intFreqBins = sum(inMatrix)';
	filteredIndices = find(intFreqBins >= obj.minInteractionFreq);

	%
	if (length(filteredIndices) > 1)

		%----- Filtering -----%
		f = inMatrix(filteredIndices,filteredIndices);
		f = f(:);
		c = c1(filteredIndices)*c1(filteredIndices)';
		c = c(:);	
		m = m1(filteredIndices)*m1(filteredIndices)';
		m = m(:);
		e = e1(filteredIndices)*e1(filteredIndices)';
		e = e(:);
		g = g1(filteredIndices)*g1(filteredIndices)';
		g = g(:);
		%
		fSum = intFreqBins(filteredIndices);
		mSum = m1(filteredIndices);
		gSum = g1(filteredIndices);
		eSum = e1(filteredIndices);
		cSum = c1(filteredIndices);		
		%
		n = length(f);
		nSum = length(filteredIndices);
	else
		%----- Filtering -----%
		f = [];
		m = [];
		g = [];
		e = [];
		c = [];
		%
		fSum = [];
		mSum = [];
		gSum = [];
		eSum = [];
		cSum = [];
		%		
		n = 0;
		nSum = 0;
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fNorm] = genomeGLMDataNormalization(fMat, f, m, g, e, cMat, nArray)
%%%%%%%%
	if (length(f) > 1)

		%%------- Standardization ---------%
		mS = (m-mean(m))/(std(m));
		gS = (g-mean(g))/(std(g));
		eS = (e-mean(e))/(std(e));		
		%	

		mS = mS - min(mS) + 1;
		gS = gS - min(gS) + 1;
		eS = eS - min(eS) + 1;		

		%%---------- Smoothing ------------%
		%warning('off','all');
		%
		%mS = splineSmooth(mS);
		%gS = splineSmooth(gS);
		%eS = splineSmooth(eS);		


		%%------- Normalization ----------%%
		%% 1) genomic biases

		features = [mS, gS, eS];
		link = 'log';
		X  = log(features);
		X2 = [ones(length(X),1), X];
		y  = f;

		%%
		stats = nbreg(X2,y);
		b = stats.b;
		yhat  = glmval(b,X,link);
		yMean = sum(yhat)/length(yhat);
		yhat  = yhat/yMean;
 
		%%%%%%%%%%%%%
		yGenome = [];
		startIndex = 1;
		for i = 1:23
			nElement = nArray(i);
			yChr = yhat(startIndex:startIndex+nElement-1);
			startIndex = startIndex + nElement;
			%%%
			yChrMatrix = yChr * yChr';
			yChrColumn = yChrMatrix(:);
			yGenome    = [yGenome; yChrColumn];
		end

		%% 2) correction
		fNorm = fMat ./sqrt(yGenome);
	else
		fNorm = [];
	end
	%%%
end
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = splineSmooth(in)
%%%%%%%%

y = in';
%
n = length(y);
x = rand(2,n);
%
st = tpaps(x,y);
ySmoothed = fnval(st,x);
%
out = ySmoothed';

end
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = nbreg( x, y, varargin )
% nbreg - Negative binomial regression 
% 
% Author: Surojit Biswas
% Email: sbiswas@live.unc.edu
% Institution: The University of North Carolina at Chapel Hill
%
% References: 
% [1] Hardin, J.W. and Hilbe, J.M. Generalized Linear Models and
% Extensions. 3rd Ed. p. 251-254. 
% 
% INPUT:
% x     =   [n x p] Design matrix. Rows are observations, Columns are
%           variables. (required)*
% y     =   [n x 1] Response matrix. Rows are observations. (required)
% 
% The following inputs are optional. For optional inputs the function call
% should include the input name as a string and then the input itself.
% For example, nbreg(x, y, 'b', myStartingRegCoeffs)
%
% alpha =   Scalar value. Starting estimate for dispersion parameter.
%           (optional, DEFAULT = 0.01)
% b     =   [p x 1] starting vector of regression coefficients. (optional,
%           DEFAULT = estimated with poisson regression)
% offset =  [n x 1] offset vector. (optional, default = 0)**
% trace =   logical value indicating if the algorithm's progress should be
%           printed to the screen. (optional, default = false)
% regularization = scalar specifying the amount of regularization (default
%                  = 1e-4). 
% distr     =   Distribution to fit ['poisson' | 'negbin' (defualt)];
% estAlpha  =   Estimate alpha or use the alpha provided? [true (default) |
%               false]
%
% OUTPUT:
% stats.b       =   [p x 1] vector of estimated regression coefficients.
% stats.alpha   =   Estimated dispersion parameter.
% stats.logL    =   Final model log-likelihood.
% stats.delta   =   Objective change at convergence. 
%
% NOTES:
% *  Supply your own column of ones for an intercept term.
% ** Make sure you have taken the log of the offset vector if necessary.
% If one wants to include an offset to control for exposure time/fit a rate
% model, then typically the offset vector should log-ed before it's
% supplied as an input. The regression routine assumes that E[y|x] ~ X*b +
% log(offset)
n = size(y,1);
p = size(x,2);
assert(n == size(x,1), 'Dimension mismatch between x and y');
df = n - p;
ALPHATHRESH = 1e-8;
%% Get variable arguments
% Offset.
offset = getVararg(varargin, 'offset');
if isempty(offset)
    offset = 0;
end
% Regression coefficients.
b = getVararg(varargin, 'b');
if isempty(b)
    mu = (y + mean(y))/2;
    eta = log(mu);
    b = Inf(size(x,2),1);
else
    eta = x*b + offset;
    mu = exp(eta);
end
% Alpha
alpha = getVararg(varargin, 'alpha');
if isempty(alpha)
    alpha = 0.01;
elseif alpha < ALPHATHRESH
    alpha = 0.01;
end
% Trace
trace = getVararg(varargin, 'trace');
if isempty(trace)
    trace = false;
end
% Regularization
reg = getVararg(varargin, 'regularization');
if isempty(reg)
    reg = 1e-4;
end
% Distribution
distr = getVararg(varargin, 'distr');
if isempty(distr)
    distr = 'negbin';
end
% Estimate alpha?
estAlpha = getVararg(varargin, 'estAlpha');
if isempty(estAlpha)
    estAlpha = true;
end
if trace
    fprintf('Delta\t');
    fprintf('alpha\t');
    for i = 1 : p
        fprintf('b_%0.0f\t', i)
    end
    fprintf('\n');
end
%% Optimize
MAXIT = 50;
eps = 1e-4;
REG = reg*eye(p);
LSOPTS.SYM = true;
LSOPTS.POSDEF = true;
del = Inf;
itr = 1;
logL = -Inf;
objChange = Inf;
while itr < MAXIT && objChange > eps
    innerItr = 0;
    while any(del > eps) && innerItr < 1000
        if strcmpi(distr, 'negbin')
            w =  repmat(mu./(1 + alpha*mu), 1, p)';     % Weight matrix
        elseif strcmpi(distr, 'poisson')
            w = repmat(mu, 1, p)';
        else
            error('Unknown distribution');
        end
        z = (eta + (y - mu)./mu) - offset;          % Working response
        xtw = x'.*w;                                % X'*W
        bold = b; 
        try
            b = linsolve(xtw*x + REG, xtw*z, LSOPTS);   % Update
        catch 
            b = linsolve(xtw*x + REG, xtw*z);
        end  
        eta = x*b + offset;                         % Linear predictor
        mu = exp(eta);                              % Mean
        del = abs(bold - b);
        innerItr = innerItr + 1;
    end
    
    if strcmpi(distr, 'poisson')
        alpha = 0;
        break
    end
    
    if ~estAlpha
        objChange = nan;
        amu = alpha*mu;
        amup1 = 1 + amu;
        logL = sum( y.*log( amu./amup1 ) - (1/alpha)*log(amup1) + gammaln(y + 1/alpha) - gammaln(y + 1) - gammaln(1/alpha));
        break
    end
    
    alpha = alphaLineSearch(y, mu, alpha);
    
    % Likelihood evaluation.
    amu = alpha*mu;
    amup1 = 1 + amu;
    plogL = logL;
    logL = sum( y.*log( amu./amup1 ) - (1/alpha)*log(amup1) + gammaln(y + 1/alpha) - gammaln(y + 1) - gammaln(1/alpha));
    objChange = (logL - plogL); if isnan(objChange); objChange = Inf; end
    
    
    % Display an update if trace switch is on.
    if trace
        fprintf('%0.6f\t', [objChange, alpha, b']);
        fprintf('\n');
    end
    if alpha < ALPHATHRESH
        break
    end    
    
    itr = itr + 1;
    del = Inf;
end
stats.alpha = alpha;
stats.b = b;
stats.logL = logL;
stats.delta = objChange;
    
    function [ alpha ] = alphaLineSearch( y, mu, alpha, varargin )
        
        gridSize = getVararg(varargin, 'gridSize');
        if isempty(gridSize)
            gridSize = 0.1;
        end
        border = 1e-10;
        interval = [Inf, -Inf];
        gflank = [Inf, -Inf];
        g = updateWindow(alpha);
        sg = g;
        % Crudely follow the gradient by stepping in that direction until you 
        % get to a point where it changes. At this point you know you have
        % bounded the domain value that maximizes the function. This is
        % assuming that the function being optimized is convex.
        while sign(sg) == sign(g)
            proposal = alpha + gridSize*sign(sg);
            % Ensure the proposal doesn't violate the constraint that alpha >
            % 0.
            if proposal <= 0
                % Check left border to see if we're looking at a corner
                % solution. If it's a corner solution (gradient at border is
                % negative) then return the border and exit. Otherwise set the
                % left edge of the interval to the border.
                [~, g] = alphaML(y,mu,border);
                if g < 0
                    alpha = border;
                    return;
                else
                    interval(1) = border;
                    gflank(1) = g;
                    break;
                end
            end
            alpha = proposal;
            sg = updateWindow(alpha);
        end
        % At this point the solution lies between interval(1) and interval(2).
        % note this means there is an interior point solution.
        while interval(2) - interval(1) > 1e-6
            %gt = sum(abs(gflank));
            %w = (gt - abs(gflank))/gt;
            alpha = mean(interval); 
            updateWindow(alpha);
        end
        alpha = mean(interval);
        
        function gr = updateWindow(a)
            [~, gr] = alphaML(y,mu,a);
            if sign(gr) > 0
                interval(1) = a;
                gflank(1) = gr;
            else
                interval(2) = a;
                gflank(2) = gr;
            end
        end
    end
    
    function [ f, g, H ] = alphaML( y, mu, a )
 
        ainv = 1/a;
        amuq = a*mu;
        amup1q = amuq + 1;
        lnamup1q = log(amup1q);
        lnamuq = log(amuq);
        % Function value:
        f = sum( y.*(lnamuq - lnamup1q) - ainv*lnamup1q + gammaln(y + ainv) - gammaln(ainv) );
        if nargout > 1
        % Gradient:
        g = (ainv^2)*sum( (lnamup1q - 1) + (a*y + 1)./(a*mu + 1) - psi(y + ainv) + psi(ainv));
        end
        if nargout > 2
        % Hessian:
        H = sum(  (ainv^3)*( (a*y + 1)./((amuq + 1).^2) - (2*a*y + 4)./(amuq + 1) - 2*lnamup1q - 3 ) ...
            + (ainv^2)*(psi(1, y + ainv) - psi(1, ainv))  );
        end
    end
    
    function a = getVararg(v, param)
        ind = find(strcmpi(v, param)) + 1;
        if isempty(ind)
            a = [];
        else
            a = v{ind};
        end
    end
        
end
