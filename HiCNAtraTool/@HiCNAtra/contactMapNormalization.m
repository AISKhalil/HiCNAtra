function obj = contactMapNormalization(obj)
% Normalzing the interaction frequencies between all genomic loci.
% Interaction matrices will be save at the output directory "outputDirectory".

if(obj.toolID == 1)
%%%% HiCNAtra %%%%
	if(obj.cisTransNormMethod == 1) 
		%Genome-wise
		if(obj.cisOnly == 1)
			cisGenomeWiseNormalization(obj);
		else
			cisGenomeWiseNormalization(obj);
			transGenomeWiseNormalization(obj);
		end
	else 
		%chromosome-wise
		cisTransChromosomeWiseNormalization(obj);
	end
	
elseif(obj.toolID == 2)
%%%%%%% HiCnorm %%%%%%%
	HiCNormNormalization(obj);
	
elseif(obj.toolID == 3)
%%%%%%% oneD %%%%%%%
	oneDNormalization(obj);

elseif(obj.toolID == 4)
%%%%%%% HiCnorm %%%%%%%
	OneDplusCnNormalization(obj);
	
end
%%%




end
%%%
