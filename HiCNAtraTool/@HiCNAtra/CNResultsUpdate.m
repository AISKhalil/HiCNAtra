function [segmentsData, regionsData] = CNResultsUpdate(obj, chrIndex, segmentsInfo, impRegionsPerSegment)
%%Writing CNAtra results to output files including CNV regions and IBs.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------- Output File ------------ %%
format long;
dir = obj.outputDirectory;
if(exist(dir,'dir') ~= 7)
	mkdir(dir);
end

runDirectory = strcat(dir,'/','CNVs_binSize',int2str(obj.binSize));
if(exist(runDirectory,'dir') ~= 7)
	mkdir(runDirectory);
end

filePath = strcat(runDirectory, '/', 'chr', int2str(chrIndex),'.txt');
fileID = fopen(filePath,'w');
segmentsData = [];
regionsData  = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------- File Updating -----------%%
binSize = obj.binSize;
noSegments = length(cell2mat(keys(segmentsInfo)));
%%%Centromeres/Telomeres

centroTeloBoundaries = obj.chrCentroTeloBoundaries(chrIndex);
centroStart  = centroTeloBoundaries(1,1);
centroEnd    = centroTeloBoundaries(1,2);
telo1Start   = centroTeloBoundaries(2,1);
telo1End     = centroTeloBoundaries(2,2);
telo2Start   = centroTeloBoundaries(3,1);
telo2End     = centroTeloBoundaries(3,2);

regionBaseID = 1;
for i = 1:noSegments
	%%--------- Segment Info ----------%%
	segmentInfo = segmentsInfo(i);
	segmentStart = segmentInfo(1);
	segmentEnd   = segmentInfo(2);
	segmentData = [i, segmentInfo(1:4)];
	%Segment Output
	segmentsData = [segmentsData; segmentData];
	
	if(i==1)
		fprintf(fileID,'%-20s %-20s %-20s %-20s %-20s \r\n','Segment-Number','Start-Bin(Kb)','Stop-Bin(Kb)', 'Segment-Width(Kb)', 'Copy-Number');
	end
	fprintf(fileID,'%-20u %-20u %-20u %-20u %-14.6f \r\n',segmentData');
	fprintf(fileID,'--------------------------------------------------------------------------------------------\r\n');
	
	%%--------- Regions/Segment Info ----------%% 
	regionsInfo = impRegionsPerSegment(i);
	[noRegions,~] = size(regionsInfo);
    
    
	if(noRegions > 0)
        regionStart    = regionsInfo(:,1);
        regionEnd      = regionsInfo(:,2);
        regionWidth    = regionsInfo(:,3);
        regionCN       = regionsInfo(:,4);
        regionCategory = regionsInfo(:,5);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		regionIDs = [regionBaseID:regionBaseID + noRegions-1]';
		regionBaseID = regionBaseID + noRegions;		
		% Distance to segment edge
		distanceFromRegionStartToSegmentStart = regionStart - segmentStart;
		distanceFromRegionEndToSegmentEnd     = segmentEnd  - regionEnd;
		distanceToSegmentEdge = abs(min(distanceFromRegionStartToSegmentStart, distanceFromRegionEndToSegmentEnd));
		% Distance to Centromere/Telomeres
		distanceFromRegionStartToTelo1End  = regionStart - telo1End;
		distanceFromRegionEndToCentroStart = abs(centroStart - regionEnd);
		distanceFromRegionStartToCentroEnd = abs(regionStart - centroEnd);
		distanceFromRegionEndToTelo2Start  = abs(telo2Start - regionEnd);	
		distanceToCentroTelo = abs(min(min(distanceFromRegionStartToTelo1End,distanceFromRegionEndToCentroStart),min(distanceFromRegionStartToCentroEnd,distanceFromRegionEndToTelo2Start)));
		%Region Output
		if(i==1)
			fprintf(fileID,'	%-5s %-12s %-12s %-12s %-10s %-10s %-20s %-20s \r\n','ID','Start(bp)','Stop(bp)', 'Width(bp)', 'Copy-Number', 'Type', 'SegmentEdge Distance(Kb)', 'Centro/Telo Distance(Kb)');
		end
		regionData = [regionIDs, regionStart, regionEnd, regionWidth, regionCN, regionCategory, distanceToSegmentEdge, distanceToCentroTelo];
		regionsData = [regionsData; regionData];

		bpRegionStart = (regionStart-1)*binSize;
		bpRegionEnd   = regionEnd*binSize - 1;
		bpRegionWidth = bpRegionEnd - bpRegionStart + 1;

		bpRegionData = [regionIDs, bpRegionStart, bpRegionEnd, bpRegionWidth, regionCN, regionCategory, distanceToSegmentEdge*binSize, distanceToCentroTelo*binSize];
		for j=1:noRegions

            switch bpRegionData(j,6)
                case 5
                    regionClass = 'Amplification';
                case -5
                    regionClass = 'Deletion';
                otherwise
                    regionClass = 'Neutral';    			
            end			
			fprintf(fileID,'	%-5u %-12u %-12u %-12u %-10.5f %s \t \t %-20d %-20d \r\n',bpRegionData(j,1),bpRegionData(j,2),bpRegionData(j,3),bpRegionData(j,4),bpRegionData(j,5),regionClass, bpRegionData(j,7),bpRegionData(j,8));
		end


		
	end		
	fprintf(fileID,'\r\n');
end

fclose(fileID);
end