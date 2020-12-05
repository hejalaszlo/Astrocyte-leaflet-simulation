function simulationRun(this, varargin)
	tSim = 22;																	% Simulation time (ms)
	if nargin > 1
		tSim = varargin{1};
	end
	tStep = 0.001;																% Time step (ms)

	%%%%%%%%%% PARAMETERS %%%%%%%%%%
	
	this.Simulation.parameters.diffusion.DGlu = 3.3e5;                          % Diffusion coefficient of Glu (nm2/ms)                                 Gavrilov et al., 2018
	% DGlu = 4e5;                                                               % Diffusion coefficient of Glu (nm2/ms)                                 Zheng et al., 2008
	this.Simulation.parameters.diffusion.DGlu = 2e5;							% Diffusion coefficient of Glu (nm2/ms)                                 
	this.Simulation.parameters.diffusion.DNa = 1.15e6;                          % Diffusion coefficient of Na (nm2/ms)									Goodman et al., 2005
	this.Simulation.parameters.diffusion.DCa = 9.9e5;			                % Diffusion coefficient of Ca (nm2/ms)									Syková and Nicholson, 2008
	this.Simulation.parameters.NGluReleased = 5000;                             % Number of Glu molecules released from a single vesicle

	this.Simulation.parameters.conc.GlueBaseline = 0.3;	                        % Baseline extracellular Glu concentration (uM)
	this.Simulation.parameters.conc.Glui = 3000;								% Intracellular Glu concentration (uM)
	this.Simulation.parameters.conc.EAAT = 1.08e-2;                             % Surface density of Glu transporters (EAAT1 + EAAT2) (molecule/nm2)	Lehre and Danbolt, 1998
	this.Simulation.parameters.conc.GABAe = 0.2;								% Extracellular GABA concentration (uM)									Savtchenko et al., 2015
	this.Simulation.parameters.conc.GABAi = 2000;								% Intracellular GABA concentration (uM)
	this.Simulation.parameters.conc.GAT = 4.875e-5;								% Surface density of GAT-3 transporters (molecule/nm2)					Melone et al. 2014
	this.Simulation.parameters.conc.NCX = 5e-4;									% Surface density of NCX (molecule/nm2)									Chu et al., 2016
	this.Simulation.parameters.conc.Nae = 140000;                               % Extracellular Na+ concentration (uM)
	this.Simulation.parameters.conc.NaiBaseline = 15000;						% Baseline intracellular Na+ concentration in astrocytes (uM)
	this.Simulation.parameters.conc.Cae = 2000;									% Extracellular Ca2+ concentration (uM)
	this.Simulation.parameters.conc.CaiBaseline = 0.1;							% Baseline intracellular Ca2+ concentration in astrocytes (uM)
	this.Simulation.parameters.conc.Ke = 3000;                                  % Extracellular K+ concentration (uM)
	this.Simulation.parameters.conc.Ki = 130000;                                % Intracellular K+ concentration in astrocytes (uM)
	this.Simulation.parameters.conc.Cle = 150000;								% Extracellular Cl- concentration (uM)
	this.Simulation.parameters.conc.Cli = 7000;									% Intracellular Cl- concentration in astrocytes (uM)
	this.Simulation.parameters.conc.He = 3.98107e-02;                           % Extracellular H+ concentration (uM) - corresponding to pH 7.4
	this.Simulation.parameters.conc.Hi = 3.98107e-02;                           % Intracellular H+ concentration (uM) - corresponding to pH 7.4

	this.Simulation.parameters.EAAT.turnover = 20;                              % Turnover rate of EAAT (1/s)                                           16 1/s for EAAT1 (Wadiche and Kavanaugh, 1998); 10–40 1/s for EAAT2 (Wadiche et al., 1995; Bergles and Jahr, 1997; Danbolt, 2001)

	% Calculated parameters
	dxGlu = sqrt(2 * this.Simulation.parameters.diffusion.DGlu * tStep);		% Mean squared displacement of a Glu molecule in 1 step
	dxNa = sqrt(2 * this.Simulation.parameters.diffusion.DNa * tStep);          % Mean squared displacement of a Na+ ion in 1 step
	dxCa = sqrt(2 * this.Simulation.parameters.diffusion.DCa * tStep);          % Mean squared displacement of a Ca2+ ion in 1 step

	% Kinetic Markov model for EAATs                                                                                                                    Bergles et al., 2002
	Q10 = 3;
	VrAstrocyte = -70;                                                          % Resting membrane potential of astrocytes (mV)
	F = 96485.332;																% Faraday's constant
	R = 8.314;																	% Gas constant
	T = 273 + 37;																% Temperature (37 C)
	% Voltage dependency factors of different kinetic rates
	this.Simulation.parameters.EAAT.Vfactor1 = exp(-0.46 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.EAAT.Vfactor7 = exp(-0.55 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.EAAT.Vfactor9 = exp(-0.4 * 0.1 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.EAAT.Vfactor15f = exp(-0.59 * 0.9 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.EAAT.Vfactor15b = exp(-0.59 * (1-0.9) * F * VrAstrocyte / 1000 / (2 * R * T));
	% Each rate is multiplied by Q10 (conversion from 22-25 C to 37 C), the time step, s to ms conversion and a voltage dependence factor.
	% Concentrations are converted to M (mol/dm3)
	this.Simulation.parameters.EAAT.k1f = Q10 * tStep * 1e-3 * 1e4 / this.Simulation.parameters.EAAT.Vfactor1 * this.Simulation.parameters.conc.Nae / 1e6;		% Dependent on Na(out)
	this.Simulation.parameters.EAAT.k1b = Q10 * tStep * 1e-3 * 1e2 * this.Simulation.parameters.EAAT.Vfactor1;
	this.Simulation.parameters.EAAT.k2f = Q10 * tStep * 1e-3 * 1e4 * this.Simulation.parameters.conc.Nae / 1e6;													% Dependent on Na(out)
	this.Simulation.parameters.EAAT.k2b = Q10 * tStep * 1e-3 * 5e2;
	this.Simulation.parameters.EAAT.k4f = Q10 * tStep * 1e-3 * 6e11 * this.Simulation.parameters.conc.He / 1e6;
	this.Simulation.parameters.EAAT.k4b = Q10 * tStep * 1e-3 * 7e2;
	this.Simulation.parameters.EAAT.k6f = Q10 * tStep * 1e-3 * 6e6;																								% Should be multiplied by Glu(out)
	this.Simulation.parameters.EAAT.k6b = Q10 * tStep * 1e-3 * 5e2;
	this.Simulation.parameters.EAAT.k7f = Q10 * tStep * 1e-3 * 1e4 / this.Simulation.parameters.EAAT.Vfactor7 * this.Simulation.parameters.conc.Nae / 1e6;		% Dependent on Na(out)              
	this.Simulation.parameters.EAAT.k7b = Q10 * tStep * 1e-3 * 1e3 * this.Simulation.parameters.EAAT.Vfactor7;
	this.Simulation.parameters.EAAT.k8f = Q10 * tStep * 1e-3 * 2e3;																								% Translocation
	this.Simulation.parameters.EAAT.k8b = Q10 * tStep * 1e-3 * 1.9e3;
	this.Simulation.parameters.EAAT.k9f = Q10 * tStep * 1e-3 * 1e3 / this.Simulation.parameters.EAAT.Vfactor9;													% Na unbinding
	this.Simulation.parameters.EAAT.k9b = Q10 * tStep * 1e-3 * 4e4 * this.Simulation.parameters.EAAT.Vfactor9;													% Should be multiplied by Na(in)
	this.Simulation.parameters.EAAT.k10f = Q10 * tStep * 1e-3 * 3e3;
	this.Simulation.parameters.EAAT.k10b = Q10 * tStep * 1e-3 * 9e10 * this.Simulation.parameters.conc.Hi / 1e6;
	this.Simulation.parameters.EAAT.k11f = Q10 * tStep * 1e-3 * 3e3;
	this.Simulation.parameters.EAAT.k11b = Q10 * tStep * 1e-3 * 1e5 * this.Simulation.parameters.conc.Glui / 1e6;												% Dependent on Glu(in)
	this.Simulation.parameters.EAAT.k12f = Q10 * tStep * 1e-3 * 1e5;																							% Na unbinding
	this.Simulation.parameters.EAAT.k12b = Q10 * tStep * 1e-3 * 2e7;																							% Should be multiplied by Na(in)
	this.Simulation.parameters.EAAT.k13f = Q10 * tStep * 1e-3 * 1e5;																							% Na unbinding
	this.Simulation.parameters.EAAT.k13b = Q10 * tStep * 1e-3 * 1e8;																							% Should be multiplied by Na(in)
	this.Simulation.parameters.EAAT.k14f = Q10 * tStep * 1e-3 * 1e6 * this.Simulation.parameters.conc.Ki / 1e6;													% Dependent on K(in)
	this.Simulation.parameters.EAAT.k14b = Q10 * tStep * 1e-3 * 1e3;
	this.Simulation.parameters.EAAT.k15f = Q10 * tStep * 1e-3 * 40 / this.Simulation.parameters.EAAT.Vfactor15f;												% Reverse transport
	this.Simulation.parameters.EAAT.k15b = Q10 * tStep * 1e-3 * 10 * this.Simulation.parameters.EAAT.Vfactor15b;
	this.Simulation.parameters.EAAT.k16f = Q10 * tStep * 1e-3 * 2e4;
	this.Simulation.parameters.EAAT.k16b = Q10 * tStep * 1e-3 * 1e6 * this.Simulation.parameters.conc.Ke / 1e6;													% Dependent on K(out)

	% Kinetic Markov model for GAT-3																															Bicho et al., 2005
	% Voltage dependency factors of different kinetic rates
	this.Simulation.parameters.GAT.Vfactor1 = exp(-0.9 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.GAT.Vfactor2 = exp(-0.1 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.GAT.Vfactor4 = exp(-0.2 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.GAT.Vfactor5 = exp(-0.2 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.GAT.Vfactor6 = exp(-0.01 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.GAT.Vfactor7 = exp(-0.005 * F * VrAstrocyte / 1000 / (2 * R * T));
	this.Simulation.parameters.GAT.Vfactor8 = exp(0.2 * F * VrAstrocyte / 1000 / (2 * R * T));
	% Each rate is multiplied by Q10 (conversion from 22-25 C to 37 C), the time step and a voltage dependence factor.
	% Concentrations are converted to M (mol/dm3)
	this.Simulation.parameters.GAT.k1f = Q10 * tStep * 0.5 * this.Simulation.parameters.GAT.Vfactor1 * this.Simulation.parameters.conc.Nae / 1e6;				% Dependent on Na(out)
	this.Simulation.parameters.GAT.k1b = Q10 * tStep * 0.01 / this.Simulation.parameters.GAT.Vfactor1;
	this.Simulation.parameters.GAT.k2f = Q10 * tStep * 10 * this.Simulation.parameters.GAT.Vfactor2 * this.Simulation.parameters.conc.Nae / 1e6;				% Dependent on Na(out)
	this.Simulation.parameters.GAT.k2b = Q10 * tStep * 0.1 / this.Simulation.parameters.GAT.Vfactor2;
	this.Simulation.parameters.GAT.k3f = Q10 * tStep * 1e4 * this.Simulation.parameters.conc.GABAe / 1e6;														% Dependent on GABA(out)
	this.Simulation.parameters.GAT.k3b = Q10 * tStep * 1;
	this.Simulation.parameters.GAT.k4f = Q10 * tStep * 1 * this.Simulation.parameters.GAT.Vfactor4;																% Translocation
	this.Simulation.parameters.GAT.k4b = Q10 * tStep * 1 / this.Simulation.parameters.GAT.Vfactor4;
	this.Simulation.parameters.GAT.k5f = Q10 * tStep * 0.3 * this.Simulation.parameters.GAT.Vfactor5;															% GABA and 2 Na unbinding
	this.Simulation.parameters.GAT.k5b = Q10 * tStep * 1e6 / this.Simulation.parameters.GAT.Vfactor5 * this.Simulation.parameters.conc.GABAi / 1e6;				% Dependent on GABA(in), should be multiplied by Na(in)^2 !!!
	this.Simulation.parameters.GAT.k6f = Q10 * tStep * 200 * this.Simulation.parameters.GAT.Vfactor6 * this.Simulation.parameters.conc.Cli / 1e6;				% Dependent on Cl(in)
	this.Simulation.parameters.GAT.k6b = Q10 * tStep * 10 / this.Simulation.parameters.GAT.Vfactor6;
	this.Simulation.parameters.GAT.k7f = Q10 * tStep * 0.4 * this.Simulation.parameters.GAT.Vfactor7;															% Reverse transport
	this.Simulation.parameters.GAT.k7b = Q10 * tStep * 0.02 / this.Simulation.parameters.GAT.Vfactor7;
	this.Simulation.parameters.GAT.k8f = Q10 * tStep * 50 * this.Simulation.parameters.GAT.Vfactor8;											
	this.Simulation.parameters.GAT.k8b = Q10 * tStep * 150 / this.Simulation.parameters.GAT.Vfactor8 * this.Simulation.parameters.conc.Cle / 1e6;				% Dependent on Cl(out)

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
	Vcleft = sum(space(:) == -1) * this.Volume.Voxel;
	Vperisynpre = sum(space(:) == -2) * this.Volume.Voxel;
	Vperisynpost = sum(space(:) == -3) * this.Volume.Voxel;
	VExtracell = sum(space(:) <= 0) * this.Volume.Voxel;
	VAstrocyte = sum(ismember(space(:), int16(this.AstrocyteId))) * this.Volume.Voxel;

	%%%%%%%%%% ADD MOLECULES %%%%%%%%%%
	posGlu = [];
	Glunum = round(this.Simulation.parameters.conc.GlueBaseline * VExtracell * 6 * 1e17); % Number of Glu molecules needed to achieve the baseline concentration
	ind = find(this.SegmentImageCorrected == 0); % Indices of all extracellular points in space
	ind = ind(randperm(numel(ind), Glunum)); % Random indices where Glu will be inserted
	[posGlu(1,:), posGlu(2,:), posGlu(3,:)] = ind2sub(size(space), ind);
	Glubound = zeros(1, Glunum); % Index of the EAAT molecule to which each Glu is bound (0 if not bound)
	clear indGlu;

	posEAAT = [];
	indEAAT = find(ismember(this.SegmentImageCorrectedSurface, this.AstrocyteId) > 0); % Indices of all astrocyte surface points in space
	EAATnum = round(this.Simulation.parameters.conc.EAAT * length(indEAAT) * this.VoxelSizeX * this.VoxelSizeZ); % Number of EAAT molecules needed to achieve the baseline concentration
	indEAAT = repmat(indEAAT, 3, 1); % To get enough coordinates for the EAATs, considering the >1 nm voxel size
	indEAAT = indEAAT(randperm(numel(indEAAT), EAATnum)); % Random indices where EAAT will be inserted
	[posEAAT(1,:), posEAAT(2,:), posEAAT(3,:)] = ind2sub(size(this.SegmentImageCorrectedSurface), indEAAT);
	clear indEAAT;
	EAATstate = ones(ceil(tSim / tStep), EAATnum, 'int8');
	% EAATstate(1,:) = ceil(rand(1, EAATnum) * 14);
	EAATstate(1,:) = arrayfun(@(x) find(x < [0.00223,0.01237,0.039949,0.92942,0.93277,0.93404,0.93551,0.93809,0.94027,0.96456,0.97384,0.97581,0.99949,1.01], 1), rand(1, EAATnum));
	boxPixelX = round(50 / this.VoxelSizeX);
	boxPixelY = round(50 / this.VoxelSizeY);
	boxPixelZ = round(50 / this.VoxelSizeZ);
	EAATVEC = arrayfun(@(x,y,z) sum(sum(sum(space(max(x-boxPixelX,1):min(x+boxPixelX,size(space, 1)), max(y-boxPixelY,1):min(y+boxPixelY,size(space, 2)), max(z-boxPixelZ,1):min(z+boxPixelZ,size(space, 3))) <= 0))), posEAAT(1,:), posEAAT(2,:), posEAAT(3,:)) * this.Volume.Voxel; % Volume of the extracellular space in the surrounding 50x50x50 nm3 cube (dm3)
	EAATcGlu = zeros(tSim, EAATnum, 'single'); % Glu concentration in the surrounding 50x50x50 nm3 of each EAAT
	parameters.EAAT.k1fConst = repmat(this.Simulation.parameters.EAAT.k1f, 1, EAATnum); % These kinetic parameters are constant
	parameters.EAAT.k2fConst = repmat(this.Simulation.parameters.EAAT.k2f, 1, EAATnum);
	parameters.EAAT.k4fConst = repmat(this.Simulation.parameters.EAAT.k4f, 1, EAATnum);
	parameters.EAAT.k7fConst = repmat(this.Simulation.parameters.EAAT.k7f, 1, EAATnum);
	parameters.EAAT.k8fConst = repmat(this.Simulation.parameters.EAAT.k8f, 1, EAATnum);
	parameters.EAAT.k9fConst = repmat(this.Simulation.parameters.EAAT.k9f, 1, EAATnum);
	parameters.EAAT.k10fConst = repmat(this.Simulation.parameters.EAAT.k10f, 1, EAATnum);
	parameters.EAAT.k11fConst = repmat(this.Simulation.parameters.EAAT.k11f, 1, EAATnum);
	parameters.EAAT.k12fConst = repmat(this.Simulation.parameters.EAAT.k12f, 1, EAATnum);
	parameters.EAAT.k13fConst = repmat(this.Simulation.parameters.EAAT.k13f, 1, EAATnum);
	parameters.EAAT.k14fConst = repmat(this.Simulation.parameters.EAAT.k14f, 1, EAATnum);
	parameters.EAAT.k15fConst = repmat(this.Simulation.parameters.EAAT.k15f, 1, EAATnum);
	parameters.EAAT.k16fConst = repmat(this.Simulation.parameters.EAAT.k16f, 1, EAATnum);
	parameters.EAAT.k1bConst = repmat(this.Simulation.parameters.EAAT.k1b, 1, EAATnum);
	parameters.EAAT.k2bConst = repmat(this.Simulation.parameters.EAAT.k2b, 1, EAATnum);
	parameters.EAAT.k4bConst = repmat(this.Simulation.parameters.EAAT.k4b, 1, EAATnum);
	parameters.EAAT.k6bConst = repmat(this.Simulation.parameters.EAAT.k6b, 1, EAATnum);
	parameters.EAAT.k7bConst = repmat(this.Simulation.parameters.EAAT.k7b, 1, EAATnum);
	parameters.EAAT.k8bConst = repmat(this.Simulation.parameters.EAAT.k8b, 1, EAATnum);
	parameters.EAAT.k10bConst = repmat(this.Simulation.parameters.EAAT.k10b, 1, EAATnum);
	parameters.EAAT.k11bConst = repmat(this.Simulation.parameters.EAAT.k11b, 1, EAATnum);
	parameters.EAAT.k14bConst = repmat(this.Simulation.parameters.EAAT.k14b, 1, EAATnum);
	parameters.EAAT.k15bConst = repmat(this.Simulation.parameters.EAAT.k15b, 1, EAATnum);
	parameters.EAAT.k16bConst = repmat(this.Simulation.parameters.EAAT.k16b, 1, EAATnum);

	posGAT = [];
	indGAT = find(ismember(this.SegmentImageCorrectedSurface, this.AstrocyteId) > 0); % Indices of all astrocyte surface points in space
	GATnum = round(this.Simulation.parameters.conc.GAT * length(indGAT) * this.VoxelSizeX * this.VoxelSizeZ); % Number of GAT molecules needed to achieve the baseline concentration
	indGAT = indGAT(randperm(numel(indGAT), GATnum)); % Random indices where GAT will be inserted
	[posGAT(1,:), posGAT(2,:), posGAT(3,:)] = ind2sub(size(this.SegmentImageCorrectedSurface), indGAT);
	clear indGAT;
	GATstate = ones(ceil(tSim / tStep), GATnum, 'int8');
	% GATstate(1,:) = ceil(rand(1, GATnum) * 8);
	GATstate(1,:) = arrayfun(@(x) find(x < [0,0.04444,1,1,1,1,1,1.01], 1), rand(1, GATnum));
	parameters.GAT.k1fConst = repmat(this.Simulation.parameters.GAT.k1f, 1, GATnum); % These kinetic parameters are constant
	parameters.GAT.k2fConst = repmat(this.Simulation.parameters.GAT.k2f, 1, GATnum);
	parameters.GAT.k3fConst = repmat(this.Simulation.parameters.GAT.k3f, 1, GATnum);
	parameters.GAT.k4fConst = repmat(this.Simulation.parameters.GAT.k4f, 1, GATnum);
	parameters.GAT.k5fConst = repmat(this.Simulation.parameters.GAT.k5f, 1, GATnum);
	parameters.GAT.k6fConst = repmat(this.Simulation.parameters.GAT.k6f, 1, GATnum);
	parameters.GAT.k7fConst = repmat(this.Simulation.parameters.GAT.k7f, 1, GATnum);
	parameters.GAT.k8fConst = repmat(this.Simulation.parameters.GAT.k8b, 1, GATnum);
	parameters.GAT.k1bConst = repmat(this.Simulation.parameters.GAT.k1b, 1, GATnum);
	parameters.GAT.k2bConst = repmat(this.Simulation.parameters.GAT.k2b, 1, GATnum);
	parameters.GAT.k3bConst = repmat(this.Simulation.parameters.GAT.k3b, 1, GATnum);
	parameters.GAT.k4bConst = repmat(this.Simulation.parameters.GAT.k4b, 1, GATnum);
	parameters.GAT.k6bConst = repmat(this.Simulation.parameters.GAT.k6b, 1, GATnum);
	parameters.GAT.k7bConst = repmat(this.Simulation.parameters.GAT.k7b, 1, GATnum);
	parameters.GAT.k8bConst = repmat(this.Simulation.parameters.GAT.k8b, 1, GATnum);

	posNCX = [];
	indNCX = find(ismember(this.SegmentImageCorrectedSurface, this.AstrocyteId) > 0); % Indices of all astrocyte surface points in space
	NCXnum = round(this.Simulation.parameters.conc.NCX * length(indNCX) * this.VoxelSizeX * this.VoxelSizeZ); % Number of NCX molecules needed to achieve the baseline concentration
	indNCX = indNCX(randperm(numel(indNCX), NCXnum)); % Random indices where NCX will be inserted
	[posNCX(1,:), posNCX(2,:), posNCX(3,:)] = ind2sub(size(this.SegmentImageCorrectedSurface), indNCX);
	clear indNCX;
	NCXstate = ones(ceil(tSim / tStep), NCXnum, 'int8');
	% NCXstate(1,:) = ceil(rand(1, NCXnum) * 6);
	NCXstate(1,:) = arrayfun(@(x) find(x < [0,0.04444,1,1,1,1,1,1.01], 1), rand(1, NCXnum));
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

	% Nanumreduction = 100;
	% Nanum = round(parameters.conc.NaiBaseline * Vastrocyte * 6 * 1e17 / Nanumreduction); % Number of Na+ ions needed to achieve the baseline concentration (divided by Nanumreduction to reduce computation time)
	% posNa = nan(3, Nanum + 20000); % The extra 20000 positions are added to accumulate the potential Na+ uptakes
	% dist = rand(1, Nanum) * 100 + cleftRadius + 20;
	% angle = rand(1, Nanum) * 360;
	% posNa(1,1:Nanum) = rand(1, Nanum) * 2600 + 200;
	% posNa(2,1:Nanum) = sizeX / 2 + cosd(angle) .* dist;
	% posNa(3,1:Nanum) = sizeZ / 2 + sind(angle) .* dist;
	% clear dist angle;
	Nanum2 = round(this.Simulation.parameters.conc.NaiBaseline * VAstrocyte * 6 * 1e17); % Number of Na+ ions needed to achieve the baseline concentration

	posCa = [];
	indCa = find(ismember(this.SegmentImageCorrected, this.AstrocyteId) > 0); % Indices of all astrocyte points in space
	Canum = round(this.Simulation.parameters.conc.CaiBaseline * this.Volume.Astrocyte * 6 * 1e17); % Number of Ca2+ ions needed to achieve the baseline concentration
	indCa = indCa(randperm(numel(indCa), Canum)); % Random indices where Ca will be inserted
	[posCa(1,:), posCa(2,:), posCa(3,:)] = ind2sub(size(this.SegmentImageCorrected), indCa);
	clear indCa;
% 	posCa = zeros(3, Canum);
% 	dist = rand(1, Canum) * 100 + cleftRadius + 20;
% 	angle = rand(1, Canum) * 360;
% 	posCa(1,:) = rand(1, Canum) * 2600 + 200;
% 	posCa(2,:) = sizeX / 2 + cosd(angle) .* dist;
% 	posCa(3,:) = sizeZ / 2 + sind(angle) .* dist;
% 	clear dist angle;
	Canum2 = round(this.Simulation.parameters.conc.CaiBaseline * VAstrocyte * 6 * 1e17);

	% Starting concentrations
	cGluCleft = zeros(ceil(tSim / tStep), 1);                                   % Glu concentration in the synaptic cleft (uM)
	cGluPerisynPre = zeros(ceil(tSim / tStep), 1);                              % Glu concentration in the pre-perisynaptic region (uM)
	cGluPerisynPost = zeros(ceil(tSim / tStep), 1);                             % Glu concentration in the post-perisynaptic region (uM)
	cGluEC = zeros(ceil(tSim / tStep), 1);										% Glu concentration in the extracellular space (uM)
	cNaAstrocyte = zeros(ceil(tSim / tStep), 1);								% Na+ concentration in the "whole" astrocyte (uM)
	cCaAstrocyte = zeros(ceil(tSim / tStep), 1);								% Ca2+ concentration in the "whole" astrocyte (uM)

	% Prepare environment
% 	hpGlu = plot3(posGlu(1,:) , posGlu(2,:), posGlu(3,:), 'g.');
% 	% hold on;
% 	% hpNa = plot3(posNa(1,1:10:end) , posNa(2,1:10:end), posNa(3,1:10:end), 'm.');
% 	% hold off;
% 	hold on;
% 	hpCa = plot3(posCa(1,:) , posCa(2,:), posCa(3,:), 'r.');
% 	hold off;
% % 	xlim([1 size(this.SegmentImageCorrected, 1)]);
% 	xlabel('Y');
% % 	ylim([1 size(this.SegmentImageCorrected, 2)]);
% 	ylabel('X');
% % 	zlim([1 size(this.SegmentImageCorrected, 3)]);
% 	zlabel('Z');
	
	%%%%%%%%%% SIMULATION %%%%%%%%%%
	astrocyteId = int16(this.AstrocyteId);
	rng(1, 'simdTwister');
	tic;
	for i = 1:ceil(tSim / tStep)
		% ----------- Glu -----------------
		% Release Glu after 10 ms (this was introduced to let the transporters to reach steady-state)
		if i * tStep == 10000%rem(i * tStep, 10) == 0
			posGluReleased = ones(3, this.Simulation.parameters.NGluReleased) .* round(size(space)' / 2); % Center of the space
			posGlu = horzcat(posGlu, posGluReleased);
			Glubound = horzcat(Glubound, zeros(1, this.Simulation.parameters.NGluReleased));
		end

		% New positions for Glu (max, min is used to avoid jumping molecules through the leaflet - would be better to reduce tStep to 0.1 us)
% 		posGluNew = posGlu + max(min(dxGlu * randn(3, size(posGlu, 2)) ./ [this.VoxelSizeX; this.VoxelSizeY; this.VoxelSizeZ], 45), -45);
		posGluNew = posGlu + dxGlu * randn(3, size(posGlu, 2)) ./ [this.VoxelSizeX; this.VoxelSizeY; this.VoxelSizeZ];

		% Move bound Glu molecules back
		posGluNew(:,Glubound > 0) = posGlu(:,Glubound > 0);

		% Handle particles that left the space
		[~, col] = find(sum(round(posGluNew) <= 0) > 0 | round(posGluNew(1,:)) > size(space, 1) | round(posGluNew(2,:)) > size(space, 2) | round(posGluNew(3,:)) > size(space, 3));
		if i > 1 && cGluEC(i-1) > this.Simulation.parameters.conc.GlueBaseline % Remove them if [Glu]e does not drop below baseline
			posGlu(:,col) = [];
			posGluNew(:,col) = [];
			Glubound(:,col) = [];
		else % Move particles back if [Glu]e drops below baseline
			posGluNew(:,col) = posGlu(:,col);
		end

		% Avoid particles moving to another compartment
	%     ind = sub2ind(size(space), round(posGlu(1,:)), round(posGlu(2,:)), round(posGlu(3,:)));
		indNew = sub2ind(size(space), round(posGluNew(1,:)), round(posGluNew(2,:)), round(posGluNew(3,:)));
	%     posGluNew(:,space(ind) < 50 & space(indNew) >= 50) = posGlu(:,space(ind) < 50 & space(indNew) >= 50);
	%     posGluNew(:,space(ind) >= 50 & space(indNew) < 50) = posGlu(:,space(ind) >= 50 & space(indNew) < 50);
		posGluNew(:,space(indNew) > 0) = posGlu(:,space(indNew) > 0);

		posGlu = posGluNew;

		% Concentrations
		ind = sub2ind(size(space), round(posGlu(1,:)), round(posGlu(2,:)), round(posGlu(3,:)));
		cGluCleft(i) = sum(space(ind) == -1) / Vcleft / 6 / 1e17;
		cGluPerisynPre(i) = sum(space(ind) == -2) / Vperisynpre / 6 / 1e17;
		cGluPerisynPost(i) = sum(space(ind) == -3) / Vperisynpost / 6 / 1e17;
		cGluEC(i) = sum(space(ind) <= 0) / VExtracell / 6 / 1e17;

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
		posCaNew(:,ismember(space(indNew), astrocyteId) == 0) = posCa(:,ismember(space(indNew), astrocyteId) == 0);

		posCa = posCaNew;
		
		% Concentrations
		ind = sub2ind(size(space), round(posCa(1,:)), round(posCa(2,:)), round(posCa(3,:)));
		cCaAstrocyte(i) = sum(ismember(space(ind), astrocyteId)) / VAstrocyte / 6 / 1e17;
		
		% ----------- EAAT -----------------
		% Calculate the number of Glu molecules around each EAAT
		posGluTemp = posGlu(:,Glubound == 0);
% 		EAATnumGlu = zeros(1, EAATnum);
		EAATnumGlu = (sum(abs(repmat(posGluTemp(1,:), size(posEAAT, 2), 1) - posEAAT(1,:)') <= 50 & abs(repmat(posGluTemp(2,:), size(posEAAT, 2), 1) - posEAAT(2,:)') <= 50 & abs(repmat(posGluTemp(3,:), size(posEAAT, 2), 1) - posEAAT(3,:)') <= 50, 2))';
% 		if size(posGluTemp, 2) <= 1 % It is faster to use a normal for loop of Glu molecules if there are not so many
% 			for iGlu = 1:size(posGluTemp, 2)
% 				ind = find(abs(posGluTemp(1,iGlu) - posEAAT(1,:)) <= (50 / this.VoxelSizeX) & abs(posGluTemp(2,iGlu) - posEAAT(2,:)) <= (50 / this.VoxelSizeY) & abs(posGluTemp(3,iGlu) - posEAAT(3,:)) <= (50 / this.VoxelSizeZ));
% 				EAATnumGlu(ind) = EAATnumGlu(ind) + 1;
% 			end
% 		else % If too many Glu molecules are present, let's use a parfor loop of EAATs
% 			EAATnumGlu = (sum(abs(repmat(posGluTemp(1,:), size(posEAAT, 2), 1) - posEAAT(1,:)') <= 50 & abs(repmat(posGluTemp(2,:), size(posEAAT, 2), 1) - posEAAT(2,:)') <= 50 & abs(repmat(posGluTemp(3,:), size(posEAAT, 2), 1) - posEAAT(3,:)') <= 50, 2))';
% 			posGluTemp1 = posGluTemp(1,:);
% 			posGluTemp2 = posGluTemp(2,:);
% 			posGluTemp3 = posGluTemp(3,:);
% 			parfor iEAAT = 1:size(posEAAT, 2)
% 				EAATnumGlu(iEAAT) = sum(abs(posGluTemp1 - posEAAT(1,iEAAT)) <= 50 & abs(posGluTemp2 - posEAAT(2,iEAAT)) <= 50 & abs(posGluTemp3 - posEAAT(3,iEAAT)) <= 50);
% 			end
% 		end

	%     EAATnumGlu = arrayfun(@(x,y,z) sum(abs(posGlu(1,Glubound == 0) - x) <= 50 & abs(posGlu(2,Glubound == 0) - y) <= 50 & abs(posGlu(3,Glubound == 0) - z) <= 50), EAATpos(1,:), EAATpos(2,:), EAATpos(3,:)); % The number of Glu molecules in the 100x100x100 nm3 cube around each EAAT
		EAATcGluCurrent = EAATnumGlu ./ EAATVEC / 6 / 1e17;
		if rem(i * tStep, 1) == 0
			EAATcGlu(i * tStep,:) = EAATcGluCurrent;
		end

		if i > 1
			% Dynamic rate constants
			parameters.EAAT.k6fCurrent = this.Simulation.parameters.EAAT.k6f * EAATcGluCurrent / 1e6; % Concentrations are converted to M (mol/dm3)
			parameters.EAAT.k9bCurrent = repmat(this.Simulation.parameters.EAAT.k9b, 1, EAATnum) * cNaAstrocyte(i - 1) / 1e6;
			parameters.EAAT.k12bCurrent = repmat(this.Simulation.parameters.EAAT.k12b, 1, EAATnum) * cNaAstrocyte(i - 1) / 1e6;
			parameters.EAAT.k13bCurrent = repmat(this.Simulation.parameters.EAAT.k13b, 1, EAATnum) * cNaAstrocyte(i - 1) / 1e6;
			EAATkForward = [parameters.EAAT.k1fConst; parameters.EAAT.k2fConst; parameters.EAAT.k4fConst; parameters.EAAT.k6fCurrent; parameters.EAAT.k7fConst; parameters.EAAT.k8fConst; parameters.EAAT.k9fConst; parameters.EAAT.k10fConst; parameters.EAAT.k11fConst; parameters.EAAT.k12fConst; parameters.EAAT.k13fConst; parameters.EAAT.k14fConst; parameters.EAAT.k15fConst; parameters.EAAT.k16fConst];
			EAATkBackward = [parameters.EAAT.k16bConst; parameters.EAAT.k1bConst; parameters.EAAT.k2bConst; parameters.EAAT.k4bConst; parameters.EAAT.k6bConst; parameters.EAAT.k7bConst; parameters.EAAT.k8bConst; parameters.EAAT.k9bCurrent; parameters.EAAT.k10bConst; parameters.EAAT.k11bConst; parameters.EAAT.k12bCurrent; parameters.EAAT.k13bCurrent; parameters.EAAT.k14bConst; parameters.EAAT.k15bConst];

			% Probabilities to determine the transport direction
			pForward = rand(14, EAATnum);
			pBackward = rand(14, EAATnum);

			nextState = int8(pForward <= EAATkForward & (pBackward > EAATkBackward | pForward ./ EAATkForward < pBackward ./ EAATkBackward)); % Whether each EAAT cycle can move to the next state
			nextState = EAATstate(i-1,:) + nextState(sub2ind(size(nextState), EAATstate(i-1,:), 1:EAATnum)); % The actual next state of each EAAT
			nextState(nextState == 15) = 1;
			EAATstate(i,:) = nextState;

			previousState = int8(pBackward <= EAATkBackward & (pForward > EAATkForward | pForward ./ EAATkForward > pBackward ./ EAATkBackward)); % Whether each EAAT cycle can move to the previous state
			previousState = EAATstate(i,:) - previousState(sub2ind(size(previousState), EAATstate(i-1,:), 1:EAATnum)); % The actual previous state of each EAAT
			previousState(previousState == 0) = 14;
			EAATstate(i,:) = previousState;

			% A Glu molecule was bound extracellularly
			if sum(EAATstate(i,:) == 5 & EAATstate(i,:) > EAATstate(i-1,:)) > 0
				indBoundEAAT = find(EAATstate(i,:) == 5 & EAATstate(i,:) > EAATstate(i-1,:));
				for iEAAT = 1:length(indBoundEAAT)
					posGluTemp = posGlu;
					posGluTemp(:,Glubound > 0) = NaN; % Mask the currently bound Glu molecules
					indBoundGlu = find(abs(posGluTemp(1,:) - posEAAT(1,indBoundEAAT(iEAAT))) <= 50 & abs(posGluTemp(2,:) - posEAAT(2,indBoundEAAT(iEAAT))) <= 50 & abs(posGluTemp(3,:) - posEAAT(3,indBoundEAAT(iEAAT))) <= 50, 1);
					if ~isempty(indBoundGlu)
						Glubound(indBoundGlu) = indBoundEAAT(iEAAT);
						posGlu(:,indBoundGlu) = posEAAT(:,indBoundEAAT(iEAAT)); % Place Glu on EAAT
					end
				end
			end
			% A Glu molecule was unbound extracellularly
			if sum(EAATstate(i,:) == 4 & EAATstate(i,:) < EAATstate(i-1,:)) > 0
				indBoundEAAT = find(EAATstate(i,:) == 4 & EAATstate(i,:) < EAATstate(i-1,:));
	% 			indBoundGlu = arrayfun(@(x) find(posGlu(1,:) == EAATpos(1,x) & posGlu(2,:) == EAATpos(2,x) & posGlu(3,:) == EAATpos(3,x), 1), indBoundEAAT);
				for iEAAT = 1:length(indBoundEAAT)
					Glubound(Glubound == indBoundEAAT(iEAAT)) = 0;
				end
			end
			% A Glu molecule was transported to the intracellular space
			if sum(EAATstate(i,:) == 7 & EAATstate(i,:) > EAATstate(i-1,:)) > 0
				indBoundEAAT = find(EAATstate(i,:) == 7 & EAATstate(i,:) > EAATstate(i-1,:));
				for iEAAT = 1:length(indBoundEAAT)
					posGlu(:,Glubound == indBoundEAAT(iEAAT)) = [];
					Glubound(:,Glubound == indBoundEAAT(iEAAT)) = [];
				end
			end
			% A Glu molecule was reverse transported to the extracellular space
			if sum(EAATstate(i,:) == 6 & EAATstate(i,:) < EAATstate(i-1,:)) > 0
				indBoundEAAT = find(EAATstate(i,:) == 6 & EAATstate(i,:) < EAATstate(i-1,:));
				for iEAAT = 1:length(indBoundEAAT)
					posGlu(:,end+1) = posEAAT(:,indBoundEAAT(iEAAT)); % Place a new Glu on EAAT
					Glubound(:,end+1) = indBoundEAAT(iEAAT);
				end
			end

	% 		% A Na+ was transported to the intracellular space
	% 		if sum((EAATstate(i,:) == 8 | EAATstate(i,:) == 11 | EAATstate(i,:) == 12) & EAATstate(i,:) > EAATstate(i-1,:)) > 0
	% 			indBoundEAAT = find((EAATstate(i,:) == 8 | EAATstate(i,:) == 11 | EAATstate(i,:) == 12) & EAATstate(i,:) > EAATstate(i-1,:));
	% 			posNa(:,find(isnan(posNa(1,:)), ceil(length(indBoundEAAT) / Nanumreduction))) = posEAAT(:,indBoundEAAT(1:Nanumreduction:end)); % Place a new Na+ on EAAT
	% % 			for iEAAT = 1:length(indBoundEAAT) / Nanumreduction
	% % 				posNa(:,find(isnan(posNa(1,:)), 1)) = posEAAT(:,indBoundEAAT(iEAAT)); % Place a new Na+ on EAAT
	% % 			end
	% 		end
	% 		% A Na+ was transported to the extracellular space
	% 		if sum((EAATstate(i,:) == 7 | EAATstate(i,:) == 10 | EAATstate(i,:) == 11) & EAATstate(i,:) < EAATstate(i-1,:)) > 0
	% % 			indBoundEAAT = find((EAATstate(i,:) == 7 | EAATstate(i,:) == 10 | EAATstate(i,:) == 11) & EAATstate(i,:) < EAATstate(i-1,:));
	% 			ind = find(~isnan(posNa(1,:)), sum((EAATstate(i,:) == 7 | EAATstate(i,:) == 10 | EAATstate(i,:) == 11) & EAATstate(i,:) < EAATstate(i-1,:)));
	% 			posNa(:,ind(1:Nanumreduction:end)) = NaN;
	% % 			for iEAAT = 1:length(indBoundEAAT) / Nanumreduction
	% % 				% TODO: this should be made better
	% % 				if size(posNa, 2) > 0
	% % 					posNa(:,end) = [];
	% % 				end
	% % 			end
	% 		end
			Nanum2 = Nanum2 + sum((EAATstate(i,:) == 8 | EAATstate(i,:) == 11 | EAATstate(i,:) == 12) & EAATstate(i,:) > EAATstate(i-1,:)) - sum((EAATstate(i,:) == 7 | EAATstate(i,:) == 10 | EAATstate(i,:) == 11) & EAATstate(i,:) < EAATstate(i-1,:));
		end
	
		cNaAstrocyte(i) = Nanum2 / VAstrocyte / 6 / 1e17;
		cCaAstrocyte(i) = size(posCa, 2) / VAstrocyte / 6 / 1e17;

		% ----------- GAT -----------------
		if i > 1
			% Dynamic rate constants
			parameters.GAT.k5bCurrent = repmat(this.Simulation.parameters.GAT.k5b, 1, GATnum) * (cNaAstrocyte(i - 1) / 1e6) ^ 2; % Concentrations are converted to M (mol/dm3)
			GATkForward = [parameters.GAT.k1fConst; parameters.GAT.k2fConst; parameters.GAT.k3fConst; parameters.GAT.k4fConst; parameters.GAT.k5fConst; parameters.GAT.k6fConst; parameters.GAT.k7fConst; parameters.GAT.k8fConst];
			GATkBackward = [parameters.GAT.k8bConst; parameters.GAT.k1bConst; parameters.GAT.k2bConst; parameters.GAT.k3bConst; parameters.GAT.k4bConst; parameters.GAT.k5bCurrent; parameters.GAT.k6bConst; parameters.GAT.k7bConst];

			% Probabilities to determine the transport direction
			pForward = rand(8, GATnum);
			pBackward = rand(8, GATnum);

			nextState = int8(pForward <= GATkForward & (pBackward > GATkBackward | pForward ./ GATkForward < pBackward ./ GATkBackward)); % Whether each GAT cycle can move to the next state
			nextState = GATstate(i-1,:) + nextState(sub2ind(size(nextState), GATstate(i-1,:), 1:GATnum)); % The actual next state of each GAT
			nextState(nextState == 9) = 1;
			GATstate(i,:) = nextState;

			previousState = int8(pBackward <= GATkBackward & (pForward > GATkForward | pForward ./ GATkForward > pBackward ./ GATkBackward)); % Whether each EAAT cycle can move to the previous state
			previousState = GATstate(i,:) - previousState(sub2ind(size(previousState), GATstate(i-1,:), 1:GATnum)); % The actual previous state of each GAT
			previousState(previousState == 0) = 8;
			GATstate(i,:) = previousState;
		end

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

			nextState = int8(pForward <= NCXkForward & (pBackward > NCXkBackward | pForward ./ NCXkForward < pBackward ./ NCXkBackward)); % Whether each NCX cycle can move to the next state
			nextState = NCXstate(i-1,:) + nextState(sub2ind(size(nextState), NCXstate(i-1,:), 1:NCXnum)); % The actual next state of each NCX
			nextState(nextState == 7) = 1;
			NCXstate(i,:) = nextState;

			previousState = int8(pBackward <= NCXkBackward & (pForward > NCXkForward | pForward ./ NCXkForward > pBackward ./ NCXkBackward)); % Whether each NCX cycle can move to the previous state
			previousState = NCXstate(i,:) - previousState(sub2ind(size(previousState), NCXstate(i-1,:), 1:NCXnum)); % The actual previous state of each NCX
			previousState(previousState == 0) = 6;
			NCXstate(i,:) = previousState;

			% A Ca2+ was bound from the intracellular space
			if sum(NCXstate(i,:) == 4 & NCXstate(i,:) < NCXstate(i-1,:)) > 0
				indBoundNCX = find(NCXstate(i,:) == 4 & NCXstate(i,:) < NCXstate(i-1,:));
				for iNCX = 1:length(indBoundNCX)
					posCa(:,end+1) = posNCX(:,indBoundNCX(iNCX)); % Place a new Ca2+ on NCX
				end
			end

			% A Ca2+ was released to the intracellular space
			if sum(NCXstate(i,:) == 5 & NCXstate(i,:) > NCXstate(i-1,:)) > 0
				indBoundNCX = find(NCXstate(i,:) == 5 & NCXstate(i,:) > NCXstate(i-1,:));
				for iNCX = 1:length(indBoundNCX)
					% TODO: this should be made better
					if size(posCa, 2) > 0
						posCa(:,end) = [];
					end
				end
			end

	% 		% 3 Na+ were transported to the intracellular space
	% 		if sum(NCXstate(i,:) == 4 & NCXstate(i,:) > NCXstate(i-1,:)) > 0
	% 			indBoundNCX = find(NCXstate(i,:) == 4 & NCXstate(i,:) > NCXstate(i-1,:));
	% 			posNa(:,find(isnan(posNa(1,:)), 3 * ceil(length(indBoundNCX) / Nanumreduction))) = repmat(posNCX(:,indBoundNCX(1:Nanumreduction:end)), [1 3]); % Place 3 new Na+ on EAAT
	% % 			for iNCX = 1:length(indBoundNCX) / Nanumreduction
	% % 				posNa(:,find(isnan(posNa(1,:)), 3)) = repmat(posNCX(:,indBoundNCX(iNCX * Nanumreduction)), [1 3]); % Place 3 new Na+ on NCX
	% % 			end
	% 		end
	% 		
	% 		% 3 Na+ were transported to the extracellular space
	% 		if sum(NCXstate(i,:) == 3 & NCXstate(i,:) < NCXstate(i-1,:)) > 0
	% % 			indBoundNCX = find(NCXstate(i,:) == 3 & NCXstate(i,:) < NCXstate(i-1,:));
	% 			ind = find(~isnan(posNa(1,:)), 3 * sum(NCXstate(i,:) == 3 & NCXstate(i,:) < NCXstate(i-1,:)));
	% 			posNa(:,ind(1:Nanumreduction:end)) = NaN;
	% % 			posNa = posNa(:,1:end-length(indBoundNCX)+1);
	% % 			for iNCX = 1:length(indBoundNCX) * 3 / Nanumreduction
	% % 				% TODO: this should be made better
	% % 				if size(posNa, 2) > 0
	% % 					posNa(:,end) = [];
	% % 				end
	% % 			end
	% 		end
			Nanum2 = Nanum2 + sum(NCXstate(i,:) == 4 & NCXstate(i,:) > NCXstate(i-1,:)) * 3 - sum(NCXstate(i,:) == 3 & NCXstate(i,:) < NCXstate(i-1,:)) * 3;
			Canum2 = Canum2 + sum(NCXstate(i,:) == 4 & NCXstate(i,:) < NCXstate(i-1,:)) - sum(NCXstate(i,:) == 5 & NCXstate(i,:) > NCXstate(i-1,:));
		end

		% ----------- Draw -----------------
% 		set(hpGlu, 'XData', posGlu(1,:), 'YData', posGlu(2,:), 'ZData', posGlu(3,:));
% % 	    set(hpNa, 'XData', posNa(1,1:10:end), 'YData', posNa(2,1:10:end), 'ZData', posNa(3,1:10:end));
% 		set(hpCa, 'XData', posCa(1,:), 'YData', posCa(2,:), 'ZData', posCa(3,:));
% 		xlim([1 size(this.SegmentImageCorrected, 1)]);
% 		ylim([1 size(this.SegmentImageCorrected, 2)]);
% 		zlim([1 size(this.SegmentImageCorrected, 3)]);
% 		drawnow;

		% ----------- Progress -----------------
		if rem(i * tStep, 1) == 0
			fprintf('%d ms...\n', i * tStep);
			toc;
			tic;
		end
	end
	
	this.Simulation.cGluCleft = cGluCleft;
	this.Simulation.cGluPerisynPre = cGluPerisynPre;
	this.Simulation.cGluPerisynPost = cGluPerisynPost;
	this.Simulation.cGluEC = cGluEC;
	this.Simulation.cNaAstrocyte = cNaAstrocyte;
	this.Simulation.cCaAstrocyte = cCaAstrocyte;
	
	this.Simulation.EAATstate = EAATstate;
	this.Simulation.GATstate = GATstate;
	this.Simulation.NCXstate = NCXstate;
end