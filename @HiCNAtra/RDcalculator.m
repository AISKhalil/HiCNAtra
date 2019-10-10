function obj = RDcalculator(obj)
%Computing the read depth signal from the input reads.
%And filtering the low-mappability, black-listed, 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--- Checking if the bias-features are computed ---%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(obj.chrEffectiveLengthTracks.keys) == 0)
	restrictionSites(obj);
	effectiveLengthFeatureCalculater(obj);
end
%
%
if(length(obj.chrMappabilityTracks.keys) == 0)
	if(obj.RDmethod == 1)
		mappabilityRfragFeatureCalculater(obj);
	elseif(obj.RDmethod == 2)
		mappabilityBinFeatureCalculater(obj);
	elseif(obj.RDmethod == 3 || obj.RDmethod == 4)
		mappabilityWindowsScores(obj);
		mappabilityFeatureCalculater(obj);		
	end
end
%
%
if(length(obj.chrGCTracks.keys) == 0)
	if(obj.RDmethod == 1)
		gcRfragFeatureCalculater(obj);
	elseif(obj.RDmethod == 2)
		gcBinFeatureCalculater(obj);
	elseif(obj.RDmethod == 3 || obj.RDmethod == 4)
		gcWindowsScores(obj);
		gcFeatureCalculater(obj);		
	end
end
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------- Read-Depth Calculation ------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(obj.RDmethod == 1)
	readsExtractionRestrictionFragments(obj);
elseif(obj.RDmethod == 2)
	readsExtractionDanglingFragments(obj);
elseif(obj.RDmethod == 3)
	readsExtractionExactCutPositions(obj);
elseif(obj.RDmethod == 4)	
	readsExtractionRestrictionFragmentsMid(obj);
end
%
darkRegionsExtraction(obj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------- Read-Depth Normalization ---------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filtering the centromeres, telomeres, 
% black-listed, low-mappability bins.
RDfiltering(obj);

% Normalization of GC-contents, mappability, 
% and effective length features.
RDnormalization(obj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------- Parameters Estimation ------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimating the HiCNAtra parameters.
pipelineParameters(obj);


%%%
end
%%%