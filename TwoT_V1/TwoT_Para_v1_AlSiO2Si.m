function [h,dx,nx,z0,LambdaL,LambdaE,CL,gammaE,T0,G,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp,Data,fold_name,sample_name,para] = TwoT_Para_v1_AlSiO2Si()

%% sample geometry: Al / SiO2 / Si
h  = [75.1 105 50000]*1e-9;      % layer thicknesses (m)
dx = [2 5 100]*1e-9;        % spatial increments (m)
nx = round(h./dx);

z0 = zeros(sum(nx),1);
I = 1:nx(1);
z0(I) = dx(1)/2:dx(1):(nx(1)*dx(1)-dx(1)/2);
for iL = 2:length(nx)
    ln1 = sum(nx(1:iL-1));
    rn1 = sum(nx(1:iL));
    I = (ln1+1):rn1;
    z0(I) = z0(ln1) + dx(iL-1)/2 + (dx(iL)/2:dx(iL):(nx(iL)*dx(iL)-dx(iL)/2));
end

%% material properties
T0 = 295;

LambdaL = [15, 1.32, 130];      % lattice thermal conductivity [W/m/K]
LambdaE = [140, 1e-4, 1e-4];    % electron thermal conductivity [W/m/K]
CL = [2.43e6, 1.62e6, 1.64e6];  % lattice heat capacity [J/m^3/K]
gammaE = [135, 1e-6, 1e-6];      % Ce = gammaE*Te [J/m^3/K^2]
gEL = [2.4e17, 1, 1];       % electron-lattice coupling [W/m^3/K]

% interface conductance matrix [W/m^2/K], one row per interface
% columns correspond to [G(1,1), G(1,2), G(2,1), G(2,2)] for each interface
% (1,1): lattice/phonon channel, (2,2): electron channel
G = [
    3e8, 0, 0, 1e20;  % Al -- SiO2
    1e8, 0, 0, 1e20   % SiO2 -- Si
];

%% time discretization
dt_max = 5e-12;
pulseIndex = 1500;
tdelay_model = [linspace(-5e-12,2e-12,pulseIndex), 2e-12+0.05e-12:0.05e-12:80e-12, 80e-12+dt_max:dt_max:6e-9];
delta_time = diff(tdelay_model);
nt = length(delta_time);

%% define heat absorption
heatL = zeros(sum(nx),1);
heatE = zeros(sum(nx),1);

iheat = 1; % pump absorption to Al layer only
I_heat = 1:sum(nx(1:iheat));
zAbs = z0(I_heat);

Absprofile = 2;
% 1 = import absorption profile from AbsCal_AlSiO2Si.mat
%     (columns: position[nm], normalized absorbance)
% 2 = single exponential, pump from top (air/Al)
% 3 = single exponential, pump from substrate side

switch Absprofile
    case 1
        load('AbsCal_AlSiO2Si.mat','AbsCal')
        heatE(I_heat) = interp1(AbsCal(:,1)*1e-9, AbsCal(:,2), zAbs, 'linear', 'extrap');
    case 2
        opd = 40e-9;
        heatE(I_heat) = exp(-zAbs/opd);
    case 3
        opd = 12e-9;
        heatE(I_heat) = exp((zAbs-sum(h(1:iheat)))/opd);
end

%% intensity profile of pump pulse
PulseTemp = zeros(nt,1);
importPulse = 1;
% 0 = use Gaussian pulse by FWHM
% 1 = import pulse profile from PulseWidth_AlSiO2Si.mat

if importPulse == 1
    load('PulseWidth_20250530_kbj.mat','PulseWidth')
    td = PulseWidth(:,1)*1e-12;
    TempProf = PulseWidth(:,2);
    for i = 1:nt
        tt = tdelay_model(i+1);
        if tt >= td(1) && tt <= td(end)
            PulseTemp(i) = interp1(td,TempProf,tt);
        end
    end
else
    FWHM_corr = 0.6e-12;
    sigma_corr = FWHM_corr/2/sqrt(2*log(2));
    for i = 1:nt
        tt = tdelay_model(i+1);
        PulseTemp(i) = exp(-tt^2/(2*sigma_corr^2));
    end
end

%% scale source by absorbed fluence
Pump = 18.9e-3;   % [W]
Abs = 0.12;    % total absorbance
w0 = 12.0e-6;   % [m]
LensT = 0.9;
Fluence = Pump*LensT*Abs/(80e6/2*pi*w0^2); % [J/m^2]

F2 = sum(dx(1)*heatE(I_heat)) * sum(PulseTemp.*delta_time(:));
heatE(I_heat) = heatE(I_heat)*Fluence/F2;

%% load raw data (optional comparison)
fold_name = 'new_folder';
sample_name = 'AlSiO2';

Data = zeros(4);
data_name = [sample_name, '.mat'];
curr_path = fileparts(mfilename('fullpath'));
fold_path = fullfile(curr_path, fold_name);
data_path = fullfile(fold_path, data_name);

if ~exist(fold_path, 'dir')
    mkdir(fold_path)
    warning('TwoT:DataFolderCreated', 'Data folder created at %s. Place %s and re-run for model-data comparison.', fold_path, data_name);
elseif ~exist(data_path, 'file')
    warning('TwoT:DataNotFound', 'Raw data file %s not found. Simulation runs without data comparison.', data_path);
else
    disp(['TTM calculation starts with ', sample_name, ' data in ', fold_name,'.'])
    load(data_path,'RawData')
    Data = RawData;

    % optional lock-in phase correction (same style as ThreeT)
    fmod = 10.9e6;
    ii = sqrt(-1);
    ttt = RawData(:,2);
    Vin = RawData(:,3);
    Vout = RawData(:,4);
    V = (Vin + ii*Vout).*exp(-ii*2*pi*fmod*ttt);
    Data(:,3) = abs(real(V));
    Data(:,4) = -abs(imag(V));
end

%% build parameter struct
para.h = h;
para.dx = dx;
para.nx = nx;
para.z0 = z0;
para.LambdaL = LambdaL;
para.LambdaE = LambdaE;
para.CL = CL;
para.gammaE = gammaE;
para.T0 = T0;
para.G = G;
para.gEL = gEL;
para.tdelay_model = tdelay_model;
para.delta_time = delta_time;
para.nt = nt;
para.pulseIndex = pulseIndex;
para.iheat = iheat;
para.heatL = heatL;
para.heatE = heatE;
para.Fluence = Fluence;
para.PulseTemp = PulseTemp;
para.Absprofile = Absprofile;
para.importPulse = importPulse;
para.fold_name = fold_name;
para.sample_name = sample_name;
para.data_name = data_name;
para.Data = Data;

end
