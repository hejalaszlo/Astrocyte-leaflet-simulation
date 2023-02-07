function simulationRunOnlyNCX(this)
	tSim = 30;																% Simulation time (ms)
	tStep = 0.001;																% Time step (ms)

	%%%%%%%%%% PARAMETERS %%%%%%%%%%
	
	this.Simulation.parameters.diffusion.DNa = 1.15e6;                          % Diffusion coefficient of Na (nm2/ms)									Goodman et al., 2005
	this.Simulation.parameters.diffusion.DCa = 9.9e5;			                % Diffusion coefficient of Ca (nm2/ms)									Syková and Nicholson, 2008

	this.Simulation.parameters.conc.NCX = 5e-4;									% Surface density of NCX (molecule/nm2)									Chu et al., 2016
	this.Simulation.parameters.conc.Nae = 140000;                               % Extracellular Na+ concentration (uM)
	this.Simulation.parameters.conc.NaiBaseline = 15000;						% Baseline intracellular Na+ concentration in astrocytes (uM)
	this.Simulation.parameters.conc.Cae = 2000;									% Extracellular Ca2+ concentration (uM)
	this.Simulation.parameters.conc.CaiBaseline = 0.1;							% Baseline intracellular Ca2+ concentration in astrocytes (uM)

	% Calculated parameters
	dxNa = sqrt(2 * this.Simulation.parameters.diffusion.DNa * tStep);          % Mean squared displacement of a Na+ ion in 1 step
	dxCa = sqrt(2 * this.Simulation.parameters.diffusion.DCa * tStep);          % Mean squared displacement of a Ca2+ ion in 1 step

	% Kinetic Markov model for EAATs                                                                                                                    Bergles et al., 2002
	Q10 = 3;

	% Kinetic Markov model for NCX																																Chu et al., 2017
	% Each rate is multiplied by Q10 (conversion from 22-25 C to 37 C), the time step, s to ms conversion and a voltage dependence factor.
	% Concentrations are converted to M (mol/dm3)
	this.Simulation.parameters.NCX.k1f = Q10 * tStep * 1e-3 * 4.325e10 * (this.Simulation.parameters.conc.Nae / 1e6) ^ 3;										% Dependent on Na(out)
	this.Simulation.parameters.NCX.k1b = Q10 * tStep * 1e-3 * 1.0524e4;
	this.Simulation.parameters.NCX.k2f = Q10 * tStep * 1e-3 * 1.0781e4;																							% Translocation
	this.Simulation.parameters.NCX.k2b = Q10 * tStep * 1e-3 * 3.3541e3;
	this.Simulation.parameters.NCX.k3f = Q10 * tStep * 1e-3 * 1.408e4;
	this.Simulation.parameters.NCX.k3b = Q10 * tStep * 1e-3 * 2.8692e9;																							% Should be multiplied by Na(in)^3 !!!
	this.Simulation.parameters.NCX.k4f = Q10 * tStep * 1e-3 * 4.024e8;																							% Should be multiplied by Ca(in) !!!
	this.Simulation.parameters.NCX.k4b = Q10 * tStep * 1e-3 * 7.3916e3;
	this.Simulation.parameters.NCX.k5f = Q10 * tStep * 1e-3 * 7.1389e3;									
	this.Simulation.parameters.NCX.k5b = Q10 * tStep * 1e-3 * 9.4959e3;																							% Translocation
	this.Simulation.parameters.NCX.k6f = Q10 * tStep * 1e-3 * 4.5134e3;
	this.Simulation.parameters.NCX.k6b = Q10 * tStep * 1e-3 * 1.0691e10 * this.Simulation.parameters.conc.Cae / 1e6;											% Dependent on Ca(out)

	%%%%%%%%%% IDENTIFY EXTRACELLULAR SUBCOMPARTMENTS %%%%%%%%%%
	space = int16(this.SegmentImageCorrected);
	for iZ = 1:size(this.SegmentImageCorrected, 3)
		spaceZ = space(:,:,iZ);
		pre = imdilate(spaceZ == this.PredefinedProperties.AxonId, ones(9));
		post = imdilate(spaceZ == this.PredefinedProperties.SpineId, ones(9));
		% Synaptic cleft
		spaceZ(spaceZ == 0 & pre & post) = -1;
		% Perisynaptic, presynaptic
		spaceZ(spaceZ == 0 & pre & ~post) = -2;
		% Perisynaptic, postsynaptic
		spaceZ(spaceZ == 0 & ~pre & post) = -3;
		
		space(:,:,iZ) = spaceZ;
	end
	VExtracell = sum(space(:) <= 0) * this.Volume.Voxel;
	VAstrocyte = sum(ismember(space(:), int16(this.AstrocyteId))) * this.Volume.Voxel;

	%%%%%%%%%% ADD MOLECULES %%%%%%%%%%
	posNCX = [];
	indNCX = find(ismember(this.SegmentImageCorrectedSurface, this.AstrocyteId) > 0); % Indices of all astrocyte surface points in space
	NCXnum = round(this.Simulation.parameters.conc.NCX * length(indNCX) * this.VoxelSizeX * this.VoxelSizeZ); % Number of NCX molecules needed to achieve the baseline concentration
	indNCX = indNCX(randperm(numel(indNCX), NCXnum)); % Random indices where EAAT will be inserted
	[posNCX(1,:), posNCX(2,:), posNCX(3,:)] = ind2sub(size(this.SegmentImageCorrectedSurface), indNCX);
	clear indNCX;
	NCXstate = ones(1, NCXnum);
	% NCXstate(1,:) = ceil(rand(1, NCXnum) * 6);
	NCXstate = arrayfun(@(x) find(x < [0,0.04444,1,1,1,1,1,1.01], 1), rand(1, NCXnum));
	parameters.NCX.k1fConst = repmat(this.Simulation.parameters.NCX.k1f, 1, NCXnum); % These kinetic parameters are constant
	parameters.NCX.k2fConst = repmat(this.Simulation.parameters.NCX.k2f, 1, NCXnum);
	parameters.NCX.k3fConst = repmat(this.Simulation.parameters.NCX.k3f, 1, NCXnum);
	parameters.NCX.k5fConst = repmat(this.Simulation.parameters.NCX.k5f, 1, NCXnum);
	parameters.NCX.k6fConst = repmat(this.Simulation.parameters.NCX.k6f, 1, NCXnum);
	parameters.NCX.k1bConst = repmat(this.Simulation.parameters.NCX.k1b, 1, NCXnum);
	parameters.NCX.k2bConst = repmat(this.Simulation.parameters.NCX.k2b, 1, NCXnum);
	parameters.NCX.k4bConst = repmat(this.Simulation.parameters.NCX.k4b, 1, NCXnum);
	parameters.NCX.k5bConst = repmat(this.Simulation.parameters.NCX.k5b, 1, NCXnum);
	parameters.NCX.k6bConst = repmat(this.Simulation.parameters.NCX.k6b, 1, NCXnum);

	Nanum2 = round(this.Simulation.parameters.conc.NaiBaseline * VAstrocyte * 6 * 1e17); % Number of Na+ ions needed to achieve the baseline concentration

	posCa = [];
	indCa = find(ismember(this.SegmentImageCorrected, this.AstrocyteId) > 0); % Indices of all astrocyte points in space
	Canum = round(this.Simulation.parameters.conc.CaiBaseline * this.Volume.Astrocyte * 6 * 1e17); % Number of Ca2+ ions needed to achieve the baseline concentration
	indCa = indCa(randperm(numel(indCa), Canum)); % Random indices where NCX will be inserted
	[posCa(1,:), posCa(2,:), posCa(3,:)] = ind2sub(size(this.SegmentImageCorrected), indCa);
	clear indCa;
	Canum2 = round(this.Simulation.parameters.conc.CaiBaseline * VAstrocyte * 6 * 1e17);

	% Starting concentrations
	cNaAstrocyte = zeros(ceil(tSim / tStep), 1);								% Na+ concentration in the "whole" astrocyte (uM)
	cCaAstrocyte = zeros(ceil(tSim / tStep), 1);								% Ca2+ concentration in the "whole" astrocyte (uM)

	%%%%%%%%%% SIMULATION %%%%%%%%%%
	astrocyteId = int16(this.AstrocyteId);
	rng(1, 'simdTwister');
	tic;
	for i = 1:ceil(tSim / tStep)
		cNaAstrocyte(i) = this.Simulation.parameters.conc.NaiBaseline;%Nanum2 / VAstrocyte / 6 / 1e17;

		% ----------- Ca2+ -----------------
		% New positions for Ca2+
		posCaNew = posCa + dxCa * randn(3, size(posCa, 2)) ./ [this.VoxelSizeX; this.VoxelSizeY; this.VoxelSizeZ];

		% Handle particles that left the space
		[~, col] = find(sum(round(posCaNew) <= 0) > 0 | round(posCaNew(1,:)) > size(space, 1) | round(posCaNew(2,:)) > size(space, 2) | round(posCaNew(3,:)) > size(space, 3));
		if i > 1 && cCaAstrocyte(i-1) > this.Simulation.parameters.conc.CaiBaseline % Remove them if [Ca]i does not drop below baseline
			posCa(:,col) = [];
			posCaNew(:,col) = [];
		else % Move particles back if [Ca]i drops below baseline
			posCaNew(:,col) = posCa(:,col);
		end

		% Avoid particles moving to another compartment
		indNew = sub2ind(size(space), round(posCaNew(1,:)), round(posCaNew(2,:)), round(posCaNew(3,:)));
		misplacedCa = ismember(space(indNew), astrocyteId) == 0;
		posCaNew(:,misplacedCa) = posCa(:,misplacedCa);

		posCa = posCaNew;
		
		% Concentrations
		ind = sub2ind(size(space), round(posCa(1,:)), round(posCa(2,:)), round(posCa(3,:)));
		cCaAstrocyte(i) = sum(ismember(space(ind), astrocyteId)) / VAstrocyte / 6 / 1e17;
		
		% ----------- NCX -----------------
		if i > 1
			% Dynamic rate constants
			parameters.NCX.k3bCurrent = repmat(this.Simulation.parameters.NCX.k3b, 1, NCXnum) * (cNaAstrocyte(i - 1) / 1e6) ^ 3; % Concentrations are converted to M (mol/dm3)
			parameters.NCX.k4fCurrent = repmat(this.Simulation.parameters.NCX.k4f, 1, NCXnum) * cCaAstrocyte(i - 1) / 1e6; % Concentrations are converted to M (mol/dm3)
			NCXkForward = [parameters.NCX.k1fConst; parameters.NCX.k2fConst; parameters.NCX.k3fConst; parameters.NCX.k4fCurrent; parameters.NCX.k5fConst; parameters.NCX.k6fConst];
			NCXkBackward = [parameters.NCX.k6bConst; parameters.NCX.k1bConst; parameters.NCX.k2bConst; parameters.NCX.k3bCurrent; parameters.NCX.k4bConst; parameters.NCX.k5bConst];

			% Probabilities to determine the transport direction
			pForward = rand(6, NCXnum);
			pBackward = rand(6, NCXnum);

			nextState = pForward <= NCXkForward & (pBackward > NCXkBackward | pForward ./ NCXkForward < pBackward ./ NCXkBackward); % Whether each NCX cycle can move to the next state
			nextState = NCXstate + nextState(sub2ind(size(nextState), NCXstate, 1:NCXnum)); % The actual next state of each NCX
			nextState(nextState == 7) = 1;

			previousState = pBackward <= NCXkBackward & (pForward > NCXkForward | pForward ./ NCXkForward > pBackward ./ NCXkBackward); % Whether each NCX cycle can move to the previous state
			previousState = nextState - previousState(sub2ind(size(previousState), NCXstate, 1:NCXnum)); % The actual previous state of each NCX
			previousState(previousState == 0) = 6;
			newState = previousState;

			% A Ca2+ was bound from the intracellular space
			if sum(newState == 4 & newState < NCXstate) > 0
				indBoundNCX = find(newState == 4 & newState < NCXstate);
				for iNCX = 1:length(indBoundNCX)
					posCa(:,end+1) = posNCX(:,indBoundNCX(iNCX)); % Place a new Ca2+ on NCX
				end
			end

			% A Ca2+ was released to the intracellular space
			if sum(newState == 5 & newState > NCXstate) > 0
				indBoundNCX = find(newState == 5 & newState > NCXstate);
				for iNCX = 1:length(indBoundNCX)
					% TODO: this should be made better
					if size(posCa, 2) > 0
						posCa(:,end) = [];
					end
				end
			end
			
			Nanum2 = Nanum2 + sum(newState == 4 & newState > NCXstate) * 3 - sum(newState == 3 & newState < NCXstate) * 3;
			Canum2 = Canum2 + sum(newState == 4 & newState < NCXstate) - sum(newState == 5 & newState > NCXstate);

			NCXstate = newState;
		end

		% ----------- Progress -----------------
		if rem(i * tStep, 1) == 0
			fprintf('%d ms...\n', i * tStep);
			toc;
			tic;
		end
	end
	
	this.Simulation.cGluCleft = [];
	this.Simulation.cGluPerisynPre = [];
	this.Simulation.cGluPerisynPost = [];
	this.Simulation.cGluEC = [];
	this.Simulation.cNaAstrocyte = cNaAstrocyte;
	this.Simulation.cCaAstrocyte = cCaAstrocyte;
	
	this.Simulation.EAATstate = [];
	this.Simulation.GATstate = [];
	this.Simulation.NCXstate = NCXstate;
end