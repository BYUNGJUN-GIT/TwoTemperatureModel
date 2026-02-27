function [h,dx,nx,z0,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp,Data,fold_name,sample_name,para] = TwoT_Para_v1_AlSiO2Si()

%% sample geometry: Al / SiO2 / Si
h  = [80 300 5e5]*1e-9;      % layer thicknesses (m)
dx = [2 10 500]*1e-9;        % spatial increments (m)
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

LambdaL = [237, 1.4, 130];      % lattice thermal conductivity [W/m/K]
LambdaE = [180, 1e-4, 1e-4];    % electron thermal conductivity [W/m/K]
CL = [2.42e6, 1.55e6, 1.65e6];  % lattice heat capacity [J/m^3/K]
gammaE = [97, 1e-6, 1e-6];      % Ce = gammaE*Te [J/m^3/K^2]
gEL = [2.4e17, 1e6, 1e6];       % electron-lattice coupling [W/m^3/K]

% interface lattice conductance [W/m^2/K], interface i-(i+1)
Gll = [1.5e8; 1.2e8];

%% time discretization
dt_max = 5e-12;
pulseIndex = 1500;
tdelay_model = [linspace(-1.2e-12,2e-12,pulseIndex), 2e-12+0.05e-12:0.05e-12:80e-12, 80e-12+dt_max:dt_max:6e-9];
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
        opd = 12e-9;
        heatE(I_heat) = exp(-zAbs/opd);
    case 3
        opd = 12e-9;
        heatE(I_heat) = exp((zAbs-sum(h(1:iheat)))/opd);
end

%% intensity profile of pump pulse
PulseTemp = zeros(nt,1);
importPulse = 0;
% 0 = use Gaussian pulse by FWHM
% 1 = import pulse profile from PulseWidth_AlSiO2Si.mat

if importPulse == 1
    load('PulseWidth_AlSiO2Si.mat','PulseWidth')
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
Pump = 4e-3;   % [W]
Abs = 0.12;    % total absorbance
w0 = 5.0e-6;   % [m]
LensT = 0.8;
Fluence = Pump*LensT*Abs/(80e6/2*pi*w0^2); % [J/m^2]

F2 = sum(dx(1)*heatE(I_heat)) * sum(PulseTemp.*delta_time(:));
heatE(I_heat) = heatE(I_heat)*Fluence/F2;

%% load raw data (optional comparison)
fold_name = 'new_folder';
sample_name = 'AlSiO2Si_default';

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
    fmod = 2e6;
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
para.Gll = Gll;
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
