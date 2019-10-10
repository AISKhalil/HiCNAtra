function obj = restrictionSites(obj)
%Computing the places of the restriction-sites.



%%--------- updating names of chromosomes ---------%%
genomeName = obj.genome;

if (strcmp(genomeName,'hg19') || strcmp(genomeName,'hg18') || strcmp(genomeName,'hg38')) 
    chrNames = {'chr1';'chr2';'chr3';'chr4';'chr5';'chr6';'chr7';'chr8';'chr9';'chr10';'chr11';'chr12';'chr13';'chr14';'chr15';'chr16';'chr17';'chr18';'chr19';'chr20';'chr21';'chr22';'chrX'};
elseif(strcmp(genomeName,'mm9') || strcmp(genomeName,'mm10'))
    chrNames = {'chr1';'chr2';'chr3';'chr4';'chr5';'chr6';'chr7';'chr8';'chr9';'chr10';'chr11';'chr12';'chr13';'chr14';'chr15';'chr16';'chr17';'chr18';'chrX'};
else
    chrNames = {'chr1';'chr2';'chr3';'chr4';'chr5';'chr6';'chr7';'chr8';'chr9';'chr10';'chr11';'chr12';'chr13';'chr14';'chr15';'chr16';'chr17';'chr18';'chr19';'chr20';'chr21';'chr22';'chrX'};
end
obj.chrNames = chrNames;



%%----------- chromosomes for Analysis -----------%%
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



%%----------- Effective length feature -----------%%
restrictionEnzyme = obj.restrictionEnzyme;
referenceGenomeFolder = obj.referenceGenomeFolder;


for i  = 1:1:noChrs 
    chrIndex = chromosomes(i);
    chromosome = cell2mat(chrNames(i));
    FASTAfileName = strcat(referenceGenomeFolder, chromosome, '.fa')

    fidIn = fopen(FASTAfileName,'r');
    rawData = fread(fidIn,'*char')';
    fclose(fidIn);

    firstBlock = rawData(1:100);
    NL = regexp(firstBlock, '[\n]');
    firstCharIndex = NL(1)+1;% to remove first line
    filteredData = rawData(firstCharIndex:end);
    chrSequence  = strrep(filteredData,sprintf('\n'),'');
    chrLengthBps = length(chrSequence);
    clear rawData firstBlock filteredData

    % [fragments, cut_sites, fragmentLengths] = restrict(chrSequence, restrictionEnzyme);%cut_sites = hiclib_cut_rsites - 2, cut_sites = Bio.Restriction_cut_sites - 1
    % So, we use cut-sites2 = cut-sites + 2 to have the same cut_sites as hiclib (mid of the restriction-site).
    % First cut-site is 0 (not real cut). So, we also remove it. 
    % cut_sites2 are the sart bps of the restriction-fragments.
    [~, cut_sites, ~] = restrict(chrSequence,restrictionEnzyme);
    cut_sites2 = cut_sites(2:end) + 2;

    %%%
    chrRfragStart = [1; cut_sites2];
    chrRfragEnd   = [cut_sites2-1; chrLengthBps];
    chrRfragWidth = chrRfragEnd - chrRfragStart + 1;
    obj.chrLargeFragments(chrIndex) = find(chrRfragWidth >= obj.largeRfragmentWidth);

    %%%
    obj.rsites(chrIndex) = cut_sites2;
    obj.chrLengths(chrIndex) = chrLengthBps;
end    
%%%
end
