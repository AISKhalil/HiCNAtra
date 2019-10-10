function obj = transGenomeWiseNormalization(obj)


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
genomeRows = [];
genomeCols = [];
genomeEles = [];
genomeFIndex = [];

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
for i  = 1:1:noChrs 
	%
	chr1Index = chromosomes(i);
	chr1LengthBps   = chrLengths(chr1Index);
	chr1LengthBins  = ceil(chr1LengthBps/binSize);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for j = i+1:1:noChrs
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
		[f, m, g, e, c, nRow, nCol, nElement, subFilteredIndices] = normalizationFeatures(obj, filteredRawIntMatrix, mapp1Track, mapp2Track, GC1Track, GC2Track, effLen1Track, effLen2Track, CNVs1Tracks, CNVs2Tracks);
		genomeF = [genomeF; f];
		genomeM = [genomeM; m];
		genomeG = [genomeG; g];
		genomeE = [genomeE; e];
		genomeC = [genomeC; c];
		genomeRows = [genomeRows; nRow];
		genomeCols = [genomeCols; nCol];
		genomeEles = [genomeEles; nElement];
		genomeFIndex = [genomeFIndex; subFilteredIndices];
		
	end
	%%%
end
%%%


%%%%%%%%%%%%%%%%%%% Correction %%%%%%%%%%%%%%%%%%%%%%%
[cCorr,~]  = corr(genomeF, genomeC,'Type','Spearman');
[eCorr,~]  = corr(genomeF, genomeE,'Type','Spearman');
[gCorr,~]  = corr(genomeF, genomeG,'Type','Spearman');		
[mCorr,~]  = corr(genomeF, genomeM,'Type','Spearman');
preSpearmanCorrelation = [cCorr, eCorr, gCorr, mCorr];
%
[genomeFNorm] = genomeGLMDataNormalization(genomeF, genomeM, genomeG, genomeE, genomeC);
%
[cCorr,~]  = corr(genomeFNorm, genomeC,'Type','Spearman');
[eCorr,~]  = corr(genomeFNorm, genomeE,'Type','Spearman');
[gCorr,~]  = corr(genomeFNorm, genomeG,'Type','Spearman');		
[mCorr,~]  = corr(genomeFNorm, genomeM,'Type','Spearman');
postSpearmanCorrelation = [cCorr, eCorr, gCorr, mCorr];


%%%%%%%%%%%%%%%%%%%
k = 1;
startIndex = 1;
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
		%
		if(nElement >= 1)
			%
			fNorm = genomeFNorm(startIndex:startIndex + nElement -1);
			filteredIndices = genomeFIndex(startIndex:startIndex + nElement -1);
			startIndex = startIndex + nElement;			
			%%%%%
			outMatrix = zeros(nRow,nCol);
			outMatrix(filteredIndices) = fNorm;
			outMatrix = reshape(outMatrix,[nRow, nCol]);
		else
			outMatrix = [];
		end

		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		normIntFreq = zeros(chr1LengthBins, chr2LengthBins);
		normIntFreq(chr1FIndices,chr2FIndices) = outMatrix;
		

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% Saving Interaction-matrices %%
		[inputFilePath, outputFilePath] = findInOutPaths(obj, chr1Index, chr2Index);
		[row col v] = find(normIntFreq);
		dlmwrite(outputFilePath, [row col v], 'delimiter','\t');
	%%%	
	end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write correlation-results
postNormArray = [postNormArray(2:end,:); postNormArray(1,:)];
%
dir = obj.outputDirectory;
preNormFilePath = strcat(dir, '/', 'contactMapPreNorm_transGenomeWise',int2str(obj.contactMapBinSize),'.xls');
postNormFilePath = strcat(dir, '/', 'contactMapPostNorm_transGenomeWise',int2str(obj.contactMapBinSize),'.xls');
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
		chrData = zeros(chr1Length, chr2Length);
		lin_idcs = sub2ind(size(chrData), row, col);
		chrData(lin_idcs) = v;
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
%Filtering the input-chromosome to remove low-mappability, black-listed bins.
chrNewFIndex = containers.Map({1},{[]});
remove(chrNewFIndex,1);

%%-------------------- Input -------------------%%
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, m, g, e, c, nRow, nCol, nElement, filteredIndices] = normalizationFeatures(obj, inMatrix, m1, m2, g1, g2, e1, e2, c1, c2)
%%%%%%%%
	% Interaction-matrix
	[nRow, nCol] = size(inMatrix);
	%%%
	%%%
	x = [];
	y = [];
	z = [];
	for i = 1:nRow
		for j = 1:nCol
			x = [x;i];
			y = [y;j];
			z = [z; i+(j-1)*nRow];
		end
	end
	distanceIndices = z((y-x)>= 0);
	length(distanceIndices);

	intMatrix = inMatrix(:);
	filteredIndices = find(intMatrix >= obj.minInteractionFreq);
	length(filteredIndices);
	filteredIndices = intersect(filteredIndices, distanceIndices);
	length(filteredIndices);

	if (length(filteredIndices) > 1)
		%- Features -%
		m  = m1 * m2';
		m = m(:);
		g  = g1 * g2';
		g = g(:);
		e  = e1 * e2';
		e = e(:);
		c  = c1 * c2';
		c = c(:);
		
		%----- Filtering -----%
		f = intMatrix(filteredIndices);
		m = m(filteredIndices);
		g = g(filteredIndices);
		e = e(filteredIndices);
		c = c(filteredIndices);
		nElement = length(f);
	else
		%----- Filtering -----%
		f = [];
		m = [];
		g = [];
		e = [];
		c = [];
		nElement = 0;
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fNorm] = genomeGLMDataNormalization(f, m, g, e, c)
%%%%%%%%
	if (length(f) > 1)

		%%------- Standardization ---------%
		mS = (m-mean(m))/(std(m));
		gS = (g-mean(g))/(std(g));
		eS = (e-mean(e))/(std(e));		
		cS = (c-mean(c))/(std(c));
		
		mS = mS - min(mS) + 1;
		gS = gS - min(gS) + 1;
		eS = eS - min(eS) + 1;		
		cS = cS - min(cS) + 1;

		%%------- Normalization ----------%%
		if(sum(cS) == length(cS))
			features = [mS, gS, eS];
			%A = 1;
		else	
			features = [mS, gS, eS, cS];	
			%A = 2;
		end

		
		distr = 'poisson';
		link = 'log';
		X = log(features);
		y = f;

		iParam = glmfit(X,y, distr,'link', link);
		yhat   = glmval(iParam,X,link);
		fNorm = y ./ yhat;
	else
		fNorm = [];
	end
%%%
end
