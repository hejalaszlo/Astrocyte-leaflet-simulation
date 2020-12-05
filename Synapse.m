classdef Synapse < handle
    properties (SetObservable)
		% Segment data
		SegmentImage uint16
		SegmentId uint16
		SegmentVolume % in voxels
		SegmentData struct
		SegmentName string
		
		% Voxel dimensions
		VoxelSizeX
		VoxelSizeY
		VoxelSizeZ
		
		% Pre-defined synapse properties (mmc2.xls)
		PredefinedProperties
		
		% Simulation data
		Simulation
    end
    properties (Transient, SetObservable)
        % 
    end
    properties (SetAccess = private, SetObservable)
		% Astrocytic coverage of pre- and postsynaptic sites
		AstrocyteCoveragePresynaptic
		AstrocyteCoveragePostsynaptic
    end
    properties (Transient, SetAccess = private, SetObservable)
		% Segmented image properties after correction
		SegmentImageCorrected uint16
		SegmentImageCorrectedSurface uint16
		SegmentVolumeCorrected
    end
    properties (Dependent, SetObservable)
		% Volumes and surface of different compartments (dm3, nm2)
		Surface struct % nm2
		Volume struct % dm3
		SurfaceVolumeRatio struct % 1/um
		
		% PSD Centroid (in pixel)
		PSDCentroidX
		PSDCentroidY
		PSDCentroidZ
		
		% Astrocyte segment IDs
		AstrocyteId
		
		% Ratio of astrocytes
		AstrocyteRatio
		AstrocyteRatioCorrected
		AstrocyteCoverage
		
		% Ratio of extracellular space
		ECSRatio
		ECSRatioCorrected
    end
    
    methods (Static)
        % 
    end
    
    methods
        % Set methods
     
        % Get methods
% 		function value = get.SegmentImageCorrected(this)
% 			if isempty(this.SegmentImageCorrected)
% 				correctECS(this);
% 			end
% 		end
                
		function value = get.Surface(this)
			sur.Voxel = this.VoxelSizeX * this.VoxelSizeZ;
			sur.Axon = sum(sum(sum(this.SegmentImageCorrectedSurface == this.PredefinedProperties.AxonId))) * sur.Voxel;
			sur.Spine = sum(sum(sum(this.SegmentImageCorrectedSurface == this.PredefinedProperties.SpineId))) * sur.Voxel;
			sur.Astrocyte = sum(sum(sum(ismember(this.SegmentImageCorrectedSurface, this.AstrocyteId) > 0))) * sur.Voxel;
			
			value = sur;
		end
                
		function value = get.Volume(this)
			vol.Voxel = this.VoxelSizeX * this.VoxelSizeY * this.VoxelSizeZ / 1e24;
			vol.Extracell = sum(sum(sum(this.SegmentImageCorrected == 0))) * vol.Voxel;
			vol.Axon = sum(sum(sum(this.SegmentImageCorrected == this.PredefinedProperties.AxonId))) * vol.Voxel;
			vol.Spine = sum(sum(sum(this.SegmentImageCorrected == this.PredefinedProperties.SpineId))) * vol.Voxel;
			vol.Astrocyte = sum(sum(sum(ismember(this.SegmentImageCorrected, this.AstrocyteId) > 0))) * vol.Voxel;
			
			value = vol;
		end
                
		function value = get.SurfaceVolumeRatio(this)
			% Surface is in nm2, volume is in dm3, ratio is calculated as 1/um
			survol.Axon = (this.Surface.Axon / 1e6) / (this.Volume.Axon * 1e15);
			survol.Spine = (this.Surface.Spine / 1e6) / (this.Volume.Spine * 1e15);
			survol.Astrocyte = (this.Surface.Astrocyte / 1e6) / (this.Volume.Astrocyte * 1e15);
			
			value = survol;
		end
                
		function value = get.PSDCentroidX(this)
			value = str2double(this.PredefinedProperties.PSDCentroidPixel{:}(1:4));
		end
                
		function value = get.PSDCentroidY(this)
			value = str2double(this.PredefinedProperties.PSDCentroidPixel{:}(7:10));
		end
                
		function value = get.PSDCentroidZ(this)
			value = str2double(this.PredefinedProperties.PSDCentroidPixel{:}(13:16));
		end
                
		function value = get.AstrocyteId(this)
			astrocytes = find_str_cell(this.SegmentName, 'Glia') > 0;
			value = this.SegmentId(astrocytes);
		end
                
		function value = get.AstrocyteRatio(this)
			if ~isempty(this.AstrocyteId)
				value = sum(sum(sum(ismember(this.SegmentImage, this.AstrocyteId) > 0))) / numel(this.SegmentImage) * 100;
			else
				value = 0;
			end
		end
                
		function value = get.AstrocyteRatioCorrected(this)
			if ~isempty(this.AstrocyteId)
				value = sum(sum(sum(ismember(this.SegmentImageCorrected, this.AstrocyteId) > 0))) / numel(this.SegmentImageCorrected) * 100;
			else
				value = 0;
			end
		end
                
		function value = get.AstrocyteCoverage(this)
			value = this.AstrocyteCoveragePresynaptic + this.AstrocyteCoveragePostsynaptic;
		end
                
		function value = get.ECSRatio(this)
			value = sum(sum(sum(this.SegmentImage == 0))) / numel(this.SegmentImage) * 100;
		end
                
		function value = get.ECSRatioCorrected(this)
			value = sum(sum(sum(this.SegmentImageCorrected == 0))) / numel(this.SegmentImageCorrected) * 100;
		end
                
		% Simple methods
		function calculateAstrocyteCoverage(this)
			axonSurfacePoints = 0;
			axonSurfacePointsCovered = 0;
			spineSurfacePoints = 0;
			spineSurfacePointsCovered = 0;
			for iZ = 1:size(this.SegmentImageCorrected, 3)
				% Astrocyte
				a = imdilate(ismember(this.SegmentImageCorrected(:,:,iZ), this.AstrocyteId), ones(9));
				
				% Axon surface
				s = (this.SegmentImageCorrected(:,:,iZ) == this.PredefinedProperties.AxonId) - imerode(this.SegmentImageCorrected(:,:,iZ) == this.PredefinedProperties.AxonId, ones(3));
				axonSurfacePoints = axonSurfacePoints + sum(s(:));
				as = a&s;
				axonSurfacePointsCovered = axonSurfacePointsCovered + sum(as(:));
				
				% Spine surface
				s = (this.SegmentImageCorrected(:,:,iZ) == this.PredefinedProperties.SpineId) - imerode(this.SegmentImageCorrected(:,:,iZ) == this.PredefinedProperties.SpineId, ones(3));
				spineSurfacePoints = spineSurfacePoints + sum(s(:));
				as = a&s;
				spineSurfacePointsCovered = spineSurfacePointsCovered + sum(as(:));
			end
			
			this.AstrocyteCoveragePresynaptic = axonSurfacePointsCovered / axonSurfacePoints * 100;
			this.AstrocyteCoveragePostsynaptic = spineSurfacePointsCovered / spineSurfacePoints * 100;
		end
		
		function correctECS(this)
			% First, replace bouton segment Ids with their parent axon's segment Ids
			segImage = this.SegmentImage;
			for iSeg = 1:length(this.SegmentName)
				if contains(this.SegmentName(iSeg), " TB ")
					parentId = this.SegmentData(iSeg).hierarchy(1);
					segImage(segImage == this.SegmentId(iSeg)) = parentId;
				end
			end
			
			% Extracellular space is not preserved in aldehyde-fixed tissues.
			% In Kasthuri et al. 2015, ECS accounts for 6 % of the total volume, while in frozen tissues it is 15-25 % (Pallotto et al. 2015; van Harreveld and Khattab 1968; Harreveld and Fifkova 1975; Korogod et al. 2015).
			% Also, many segments touch each others, so there is no ECS available for diffusion.
			% To address these issues, we replace the borders of each segments by ECS.
			ics = false(size(segImage));
			for iSeg = 1:size(this.SegmentId, 1)
				ics = ics | imerode(segImage == this.SegmentId(iSeg) | segImage == this.SegmentData(iSeg).hierarchy(1), ones(3));
			end
			this.SegmentImageCorrected = segImage;
			this.SegmentImageCorrected(~ics) = 0;
			
			for iZ = 1:size(this.SegmentImageCorrected, 3)
				this.SegmentImageCorrectedSurface(:,:,iZ) = this.SegmentImageCorrected(:,:,iZ) - imerode(this.SegmentImageCorrected(:,:,iZ), ones(3));
			end
			
			this.SegmentVolumeCorrected = this.SegmentVolume;
			for iSeg = 1:size(this.SegmentId, 1)
				this.SegmentVolumeCorrected(iSeg) = sum(sum(sum(this.SegmentImageCorrected == this.SegmentId(iSeg))));
			end
		end
		
		function simulationExport(this, varargin)
			s = this.Simulation;
			if nargin > 1
				save(varargin{1}, '-struct', 's');
			else
				save(sprintf("%s.mat", inputname(1)), '-struct', 's');
			end
		end
		
		function simulationImport(this, varargin)
			if nargin > 1
				this.Simulation = load(varargin{1});
			else
				this.Simulation = load(sprintf("%s.mat", inputname(1)));
			end
		end
		
		function simulationClear(this)
			this.Simulation = [];
		end
		
		function implay(this)
			im = zeros(size(this.SegmentImageCorrected));
			im(this.SegmentImageCorrected > 0) = 1;
			im(this.SegmentImageCorrected == this.PredefinedProperties.AxonId) = 2;
			im(this.SegmentImageCorrected == this.PredefinedProperties.SpineId) = 3;
			if ~isempty(this.AstrocyteId)
				im(ismember(this.SegmentImageCorrected, this.AstrocyteId)) = 4;
			end
			cmap = [0.4	0.4	0.4;
					0	0.8	0;
					0.7	0	0;
					0.8	0.8	0];
					
			implay(permute(label2rgb3d(im, cmap), [1 2 4 3]));
		end
		
		function plot(this)
			figure;
			subplot(2,1,1);
			plot(this.Simulation.cCaAstrocyte, 'r');
			subplot(2,1,2);
			plot(this.Simulation.cGluEC, 'k');
		end

		function plot3(this)
			figure;
			% Axon
			[x, y, z] = ind2sub(size(this.SegmentImageCorrected), find(this.SegmentImageCorrected == this.PredefinedProperties.AxonId));
			plot3(x, y, z, 'g.');
			hold on;
			% Spine
			[x, y, z] = ind2sub(size(this.SegmentImageCorrected), find(this.SegmentImageCorrected == this.PredefinedProperties.SpineId));
			plot3(x, y, z, 'r.');
			% Astrocyte
			[x, y, z] = ind2sub(size(this.SegmentImageCorrected), find(ismember(this.SegmentImageCorrected, this.AstrocyteId)));
			plot3(x, y, z, 'y.');
		end

		function plotEAATstate(this)
			figure;
			imagesc((1:size(this.Simulation.EAATstate, 1)) / 1000, (1:size(this.Simulation.EAATstate, 2)),  this.Simulation.EAATstate');
			cmap = zeros(14,3);
			% Outward open states - red
			cmap(1:6,1) = linspace(0.3,1,6);
			cmap(14,1) = 0.2;
			% Inward open states - green
			cmap(7:13,2) = linspace(0.2,1,7);
			colormap(cmap);
			xlabel('Time (ms)');
		end

		function plotGATstate(this)
			figure;
			imagesc((1:size(this.Simulation.GATstate, 1)) / 1000, (1:size(this.Simulation.GATstate, 2)), this.Simulation.GATstate');
			cmap = zeros(8,3);
			% Outward open states - red
			cmap(1:4,1) = linspace(0.4,1,4);
			cmap(8,1) = 0.2;
			% Inward open states - green
			cmap(5:7,2) = linspace(0.2,1,3);
			colormap(cmap);
			xlabel('Time (ms)');
		end

		function plotNCXstate(this)
			figure;
			imagesc((1:size(this.Simulation.NCXstate, 1)) / 1000, (1:size(this.Simulation.NCXstate, 2)), this.Simulation.NCXstate');
			cmap = zeros(6,3);
			% Outward open states - red
			cmap(1:2,1) = linspace(0.6,1,2);
			cmap(6,1) = 0.2;
			% Inward open states - green
			cmap(3:5,2) = linspace(0.2,1,3);
			colormap(cmap);
			xlabel('Time (ms)');
		end

		% Complex methods
		simulationRun(this, varargin)
		simulationRunNoTransporterStates(this)
		simulationRunOnlyNCX(this)
		exportObjects(this, filename)
    end
end