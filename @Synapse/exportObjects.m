function exportObjects(this, filename)
	tic;

	fid = fopen(filename,'wt');

% 	fprintf(fid, '%s\n', 'mtllib Astrocyte.mtl');
% 	fprintf(fid, '%s\n', 'usemtl Astrocyte');

	% ASTROCYTE	
	fprintf(fid, 'o %s\n', 'Astrocyte');
	fv = isosurface(ismember(this.SegmentImageCorrected, this.AstrocyteId), 0.5);
	fv = smoothpatch(fv, 0, 5);
	f = fv.faces;
	v = fv.vertices;
	v(:,1) = v(:,1) * this.VoxelSizeX / 100;
	v(:,2) = v(:,2) * this.VoxelSizeY / 100;
	v(:,3) = v(:,3) * this.VoxelSizeZ / 100;
	
	for i = 1:10000:size(v,1)
		str = [];
		for j = i:min([i+9999 size(v,1)])
			str = [str sprintf('v %f %f %f\n',v(j,1),v(j,2),v(j,3))];
		end
		if min(size(str))>0
			fprintf(fid, str);
		end
	end

	for i = 1:10000:size(f,1)
		str = [];
		for j = i:min([i+9999 size(f,1)])
			str = [str sprintf('f %d %d %d\n',f(j,1),f(j,2),f(j,3))];
		end
		if min(size(str))>0
			fprintf(fid, str);
		end
	end
	
	vcount = size(v, 1);

	% BOUTON	
	fprintf(fid, 'o %s\n', 'Bouton');
	fv = isosurface(this.SegmentImageCorrected == this.PredefinedProperties.SpineId, 0.5);
	fv = smoothpatch(fv, 0, 5);
	f = fv.faces;
	f = f + vcount;
	v = fv.vertices;
	v(:,1) = v(:,1) * this.VoxelSizeX / 100;
	v(:,2) = v(:,2) * this.VoxelSizeY / 100;
	v(:,3) = v(:,3) * this.VoxelSizeZ / 100;
	
	for i = 1:10000:size(v,1)
		str = [];
		for j = i:min([i+9999 size(v,1)])
			str = [str sprintf('v %f %f %f\n',v(j,1),v(j,2),v(j,3))];
		end
		if min(size(str))>0
			fprintf(fid, str);
		end
	end

	for i = 1:10000:size(f,1)
		str = [];
		for j = i:min([i+9999 size(f,1)])
			str = [str sprintf('f %d %d %d\n',f(j,1),f(j,2),f(j,3))];
		end
		if min(size(str))>0
			fprintf(fid, str);
		end
	end

	vcount = vcount + size(v, 1);

	% SPINE	
	fprintf(fid, 'o %s\n', 'Spine');
	fv = isosurface(this.SegmentImageCorrected == this.PredefinedProperties.AxonId, 0.5);
	fv = smoothpatch(fv, 0, 5);
	f = fv.faces;
	f = f + vcount;
	v = fv.vertices;
	v(:,1) = v(:,1) * this.VoxelSizeX / 100;
	v(:,2) = v(:,2) * this.VoxelSizeY / 100;
	v(:,3) = v(:,3) * this.VoxelSizeZ / 100;
	
	for i = 1:10000:size(v,1)
		str = [];
		for j = i:min([i+9999 size(v,1)])
			str = [str sprintf('v %f %f %f\n',v(j,1),v(j,2),v(j,3))];
		end
		if min(size(str))>0
			fprintf(fid, str);
		end
	end

	for i = 1:10000:size(f,1)
		str = [];
		for j = i:min([i+9999 size(f,1)])
			str = [str sprintf('f %d %d %d\n',f(j,1),f(j,2),f(j,3))];
		end
		if min(size(str))>0
			fprintf(fid, str);
		end
	end

	fclose(fid);
	
	toc;
end
