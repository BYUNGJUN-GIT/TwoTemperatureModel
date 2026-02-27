function [h,dx,nx,z0,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,T0,G,g12,g13,g23,...
    tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,Data,fold_name,sample_name,para] = ThreeT_Para_v4_PtCoPtSapp()

%% sample properties
h = [1.5 0.5 4 50 5000]*1e-9;  % layer thicknesses (m) 1.5 0.5 4 100 5000
dx = [0.1 0.05 0.1 2 100]*1e-9; % spatial increments (m) 0.1 0.05 0.1 2 100
nx = round(h./dx);    % vector listing number of nodes of the layers

z0 = zeros(sum(nx),1);
I = 1:nx(1);
z0(I) = dx(1)/2:dx(1):(nx(1)*dx(1)-dx(1)/2); % spatial coordinate for plotting temperatures

for iL = 2:length(nx)
    ln1 = sum(nx(1:iL-1));
    ln2 = ln1+1;
    rn1 = sum(nx(1:iL));
    I = ln2:rn1;
    z0(I) = z0(ln1) + dx(iL-1)/2 + (dx(iL)/2:dx(iL):(nx(iL)*dx(iL)-dx(iL)/2));
end

T0 = 295;  %ambient temperature (K)

Lambda1 = [7 14 7 33 33];  %thermal conductivity phonons (W/m/K) 33
Lambda2 = [50 20 50 0.01 0.01];  %thermal conductivity electrons (W/m/K)
Lambda3 = [0.01 0.01 0.01 0.01 0.01];  %thermal conductivity magnons

C1 = [2.82 3.75 2.82 3.08 3.08]*1e6;  %specific heat (total) (J/m^3/K) 2.82 3.75 2.82 3.08 3.08
gamma2 = [400 680 400 1e-3 1e-3];   %specific heat coefficient of electrons (C2 = gamma*Te) (J/m^3/K^2) 250307 400
C3 = [1 0.02e6 1 1 1];  %specific heat of magnons
C1 = C1-T0*gamma2-C3;

G = [240, 0, 0,... % Pt1.5 -- Co0.5 %% 250306
     0, 8e3, 0,...
     0, 0, 0;...  %interface thermal conductance (W/m^2/K) 1->1 1->2 2->1 2->2 first interface; 1->1 1->2 2->1 2->2 second interface; etc.
     240, 0, 0,... % Co0.5 -- Pt4
     0, 8e3, 0,...
     0, 0, 0;...   % neglect thermal resistance of interface for electrons
     110, 0, 0,... % Pt4 -- Sapp 110
     0, 0, 0,...
     0, 0, 0;...
     33e3, 0, 0,... % Sapp -- Sapp 33e3
     0, 0, 0,...
     0, 0, 0]*1e6;

g12 = [6e17 2e18 6e17 1 1];  % phonon-electron coupling constant (W/m^3/K)
g13 = ones(1,length(nx));    % phonon-magnon
g23 = [1 0.5e17 1 1 1];      % electron-magnon 0.9e17 : PHYS. REV. APPLIED 13, 024007 (2020)

%% time discretization
dt_max = 4e-12;
pulseIndex = 2000;
tdelay_model = [linspace(-1.5e-12,2e-12,pulseIndex) 2e-12+0.1e-12:dt_max/5/2:200e-12 200e-12+dt_max:dt_max:10e-9];

delta_time = diff(tdelay_model);
nt = length(delta_time);

%% define heat absorption 
heat1 = zeros(1,sum(nx));
heat2 = zeros(1,sum(nx));
heat3 = zeros(1,sum(nx));

iheat = 3;    % index of the last absorption layer
I_heat = 1:sum(nx(1:iheat));
zAbs = z0(1:sum(nx(1:iheat)));

Absprofile = 1; 
% 1 = import an absorption profile, AbsCal.mat (Column 1: position, Column 2: absorbance) 
% 2 = a single exponential curve with pump incident on the top surface
% 3 = a single exponential curve with pump incident on the substrate 

switch Absprofile
    case 1
        load('AbsCal_PtCoPt_kbj.mat','AbsCal')
        heat2(I_heat) = interp1(AbsCal(:,1)*1e-9, AbsCal(:,2), zAbs);
        
    case 2
        opd = 12e-9;  % optical penetration depth        
        heat2(I_heat) = exp(-zAbs/opd);
        
    case 3
        opd = 12e-9;  % optical penetration depth
        heat2(I_heat) = exp((zAbs-sum(h(1:iheat)))/opd);   
end

%% intensity profile of a pump pulse 
PulseTemp = zeros(1,nt);
F2 = 0;

importPulse = 1; 
% 0 = use a FWHM pulse duration
% 1 = import a temporal intensity profile of pump-probe cross-correlation, PulseWidth.mat

if importPulse == 1
    load('PulseWidth_20250530_kbj.mat','PulseWidth')
    td = PulseWidth(:,1)*1e-12;
    TempProf = PulseWidth(:,2);
    for i = 2:pulseIndex
         tt = tdelay_model(i);   
         if tt < td(end)
            PulseTemp(i) = interp1(td,TempProf,tt); 
         end
         F2 = F2 + PulseTemp(i)*delta_time(i-1);
    end
    
else
    FWHM_corr = 0.824e-12;
    sigma_corr = FWHM_corr/2/sqrt(2*log(2));  
    for i = 2:pulseIndex
         tt = tdelay_model(i);   
         PulseTemp(i) = exp(-tt^2/(2*sigma_corr^2));     
         F2 = F2 + PulseTemp(i)*delta_time(i-1);
    end
end

% Pump Fluence
Pump = 4e-3; % Pump laser power in W
Abs = 0.303; % Total absorbance
w0 = 5.74e-6; % 20X, 10X: 2.96, 5.74e-6
LensT = 0.8; % 20X, 10X: 0.7, 0.8
Fluence = Pump*LensT*Abs/(80e6/2*pi*w0^2);

F2 = sum(dx(1)*heat2(I_heat))*F2;
heat2(I_heat) = heat2(I_heat)*Fluence/F2;

%% Load raw measurement data
% 250606 bj
fold_name = '250601';
sample_name = 'Co05Sapp_F11_V0145';

Data = zeros(4);
data_name = [sample_name, '.mat'];
curr_path = fileparts(mfilename('fullpath'));
fold_path = fullfile(curr_path, fold_name);
data_path = fullfile(fold_path, data_name);

if ~exist(fold_path, 'dir')
    mkdir(fold_path)
    message = 'NEW FOLDER CREATED! Please place the data in the folder and RE-RUN the code.';
    error('FolderNotFound:DataRequired', message);
elseif ~exist(data_path, 'file')
    message = 'DATAFILE NEEDED! Please place the data in the folder and RE-RUN the code.';
    error('FolderNotFound:DataRequired', message);
else
    disp(['TTM calculation starts with ', sample_name, ' data in ', fold_name,'.'])
end

load(data_path,'RawData')
Data = RawData;

fmod = 2e6;
ii = sqrt(-1);

ttt = RawData(:,2);
Vin = RawData(:,3);
Vout = RawData(:,4);
V = (Vin + ii*Vout).*exp(-ii*2*pi*fmod*ttt); % correct for pump being advanced in experiment

Data(:,3) = abs(real(V)); 
Data(:,4) = -abs(imag(V));

%% build a parameter array
para.h = h;
para.dx = dx;
para.nx = nx;
para.z0 = z0;
para.Lambda1 = Lambda1;
para.Lambda2 = Lambda2;
para.Lambda3 = Lambda3;
para.C1 = C1;
para.gamma2 = gamma2;
para.C3 = C3;
para.T0 = T0;
para.G = G;
para.g12 = g12;
para.g23 = g23;
para.tdelay_model = tdelay_model;
para.delta_time = delta_time;
para.nt = nt;
para.pulseIndex = pulseIndex;
para.iheat = iheat;
para.heat1 = heat1;
para.heat2 = heat2;
para.heat3 = heat3;
para.Fluence = Fluence;
para.PulseTemp = PulseTemp;
para.Data = Data;
para.fold_name = fold_name;
para.data_name = data_name;

