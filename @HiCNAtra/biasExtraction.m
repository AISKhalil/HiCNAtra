function obj = biasExtraction(obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parsing the bias-data (GC, mappability, reference-genome).
% and compute the GC-scores & mappability-scores for each window.

% Effective-length.
restrictionSites(obj);
effectiveLengthFeatureCalculater(obj);

% Mappability 
if(obj.RDmethod == 1)
	mappabilityRfragFeatureCalculater(obj);
elseif(obj.RDmethod == 2)
	mappabilityBinFeatureCalculater(obj);
elseif(obj.RDmethod == 3 || obj.RDmethod == 4)
	mappabilityWindowsScores(obj);
	mappabilityFeatureCalculater(obj);
end

% GC-tracks 
if(obj.RDmethod == 1)
	gcRfragFeatureCalculater(obj);
elseif(obj.RDmethod == 2)
	gcBinFeatureCalculater(obj);
elseif(obj.RDmethod == 3 || obj.RDmethod == 4)
	gcWindowsScores(obj);
	gcFeatureCalculater(obj);
end


%%%
end
