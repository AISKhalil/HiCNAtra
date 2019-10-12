%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GM06690_HiC 			   %%
%% GEO-Accession : GSM455133 	   %%
%% RefernceGenome: hg19            %%
%% Restirciton-Enzyme: hindIII	   %%
%% Read-Length: 76bps              %%
%% Molecule-length: 500            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the required information for creating HiCNAtra object based on the experiment and reference genome 
restrictionEnzyme = 'hindIII';
maximumMoleculeLength = 500;
readLength = 76;
referenceGenome = 'hg19';

% Add 'HiCNAtraTool' directory to Matlab search path
HiCNAtraDirectory = 'HiCNAtraTool';
addpath(HiCNAtraDirectory);

% Define the input HDF5 file(s)
HDF5Files = {'Example/GM06990_SRR027956_Input.hdf5'};

% Create HiCNAtra object 'GM06990_HiC' with the defined parameters
GM06990_HiC = HiCNAtra(HDF5Files, HiCNAtraDirectory, readLength, restrictionEnzyme, maximumMoleculeLength, referenceGenome);

% Set more parameters (optional)
GM06990_HiC.binSize = 500000;
GM06990_HiC.outputDirectory = './GM06690_HiC';

% run 'RD calculator' module (Pipeline stage 1)
GM06990_HiC.RDcalculator;

% run 'CNV caller' module (Pipeline stage 2)
GM06990_HiC.CNVcaller;

% run 'contact map corrector' (Pipeline stage 3) module that compute and correct the contact map
GM06990_HiC.contactMapCorrector;

% save the HiCNAtraObject, so you can load it directly next time.
save('GM06990_HiC.mat');
