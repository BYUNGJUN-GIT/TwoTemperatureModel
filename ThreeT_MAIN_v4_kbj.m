%% ThreeT_v3 package for calculation of three-temperature model (3TM).
% ThreeT_MAIN_v3.m will load parameters from ThreeT_Para_v3.m and
% simulates temperatures of (1) phonons, (2) electrons, and (3) magnons as
% a function of time and position, based on 3TM. It can also compare the model
% calculation with raw data and calculate sensitivities of magnon temperature 
% to parameters in the thermal model.
%
% T1 = phonons
% T2 = electrons
% T3 = magnons
%
% Important outputs of this script:
% T_fit : Temperatures of the designated position (Ind_Fit) as a function of time 
% Z : Goodness of fit, i.e., sum of squares of residuals. 
% Sensnorm, Sensnorm_z : Sensitivities of magnon temperatures (T3) to
%                        parameters in the thermal model.
%                         
% MAT files for optional input
% Data.mat : Measurement data (time delay, Vin, Vout)
% AbsCal.mat : Absorption profile of multilayers (position, absorbance)
% PulseWidth.mat : Temporal intensity profile of cross-correlation of pump
%                  and probe pulses (time, intensity)
%
% The original scripts were written by Richard Wilson (richard.wilson@ucr.edu)
% and developed by Johannes Kimling.
% The version 3 is released by Hyejin Jang (hjang32@illinois.edu)
% at the University of Illinois, Urbana-Champaign.
% Date: Aug 11, 2019.

%% V4 Updates(kbj)
% 1. modVin: Substracting the negative time delay offset for low f_mod(line 47, 72)
% 2. Save the scaled data, T_e, T_ph, and T_m(line 39, 116)

%% user inputx
[h,dx,nx,z0, Lambda1,Lambda2,Lambda3,C1,gamma2,C3,T0,G,g12,g13,g23,tdelay_model,delta_time,nt,...,
heat1,heat2,heat3,PulseTemp,Data,fold_name,sample_name,para] = ThreeT_Para_v4_PtCoPtSapp();

para_path = strcat(pwd,'\',fold_name,'\',['Para_', sample_name, '.mat']); 
T_path = strcat(pwd,'\',fold_name,'\', ['T_', sample_name, '.csv']);

% nodes of laser heated layers for convergence criterion in iterative calculation at each time increment
IndStart = 1;  % 1
IndEnd = sum(nx(1:3)); % 1:3

importData = 1; % Import raw data for comparison with 3TM, ex) Data_PtCoPt.mat
sim = 1;        % Calculate three temperatures
modVin = 0;     % 2025/06/05 kbj: For low f_mod, subtract the negative time delay offset

% Calculation of sensitivities takes a long time! 
sens_th = 0;   % sensitivities to thermal properties
sens_z = 0;    % sensitivities to thickness, absorption profile, pulse width.

IndFit = nx(1)+round(nx(2)/2); % position of the thermometer
% IndFit = sum(nx(1:3))+round(nx(4)/2); % position of the thermometer

%% solve heat diffusion problem
if sim
    [T1,T2,T3] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,...
        gamma2,C3,T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,...
        PulseTemp,IndStart,IndEnd);
    save(para_path)
    
    T_Fit = horzcat(1e12*tdelay_model', T1(:,IndFit), T2(:,IndFit), T3(:,IndFit));
end
%% compare measurement data and model 
if importData == 1
    
   tdelay_start = -1;  
   tdelay_end = 5;
   tdelay_norm = 5;  % T3 of Raw data wicll be scaled to match with T3 of 3TM at this time delay.
   
   % 2025/06/05 kbj
   if modVin == 1
       negtd_ind = Data(:,2) < -1;
       negtd_offset = mean(Data(negtd_ind, 3));
       Data(:,3) = Data(:,3)-negtd_offset;
   end

   [~,t_start] = min(abs(Data(:,2) - tdelay_start));
   [~,t_end]   = min(abs(Data(:,2) - tdelay_end));
   [~,t_data]  = min(abs(Data(:,2) - tdelay_norm));
   [~,t_model] = min(abs(tdelay_model - tdelay_norm*1e-12));
   
   Data_norm = Data(t_start:t_end,3)./Data(t_data,3)*T3(t_model,IndFit);
   Data_norm = horzcat(Data(t_start:t_end,2),Data_norm);

   % T1 & T2 & T3 plots
   T1_model = interp1(T_Fit(:,1),T_Fit(:,2),Data_norm(:,1));
   T2_model = interp1(T_Fit(:,1),T_Fit(:,3),Data_norm(:,1));
   T3_model = interp1(T_Fit(:,1),T_Fit(:,4),Data_norm(:,1));
   res = (T3_model - Data_norm(:,2)).^2;
   Z = sum(res)/length(Data_norm(:,2));
   disp(Z)
   DataFit = horzcat(Data_norm,T2_model,res);
   
   figure(1)
   plot(Data_norm(:,1),Data_norm(:,2),'ok',Data_norm(:,1),T2_model,'-r',...
       Data_norm(:,1), T1_model, '-b',...
       Data_norm(:,1), T3_model, '-g')
   xlabel('Time delay (ps)')
   ylabel('\DeltaT (K)')
   xlim([tdelay_start tdelay_end])
   ylim([0 inf])
   legend('Data', 'T2 (Electrons)', 'T1 (Phonons)', 'T3 (Magnons)');
   
   figure(2)
   loglog(Data(:,2),Data(:,3)./Data(t_data,3)*T2(t_model,IndFit),'ok',...
       tdelay_model*1e12,T2(:,IndFit),'-r', tdelay_model*1e12, T1(:, IndFit),'-b',...
       tdelay_model*1e12, T3(:, IndFit),'-g')
   xlabel('Time delay (ps)')
   ylabel('\DeltaT (K)')
   xlim([0.1 inf])
   ylim([min(T3(:, IndFit)), inf])
   legend('Data', 'T2 (Electrons)', 'T1 (Phonons)', 'T3 (Magnons)');

   % 2025/06/06 kbj
   header = {'t_d', 'Normalized Data', 't_d', 'ΔT_ph', 'ΔT_e', 'ΔT_m'};
   Data_norm = [Data_norm; zeros(size(T_Fit,1)- size(Data_norm,1), size(Data_norm,2))];
   DataNorm_TFit = [Data_norm, T_Fit];
   DataNorm_TFit_H = [header; num2cell(DataNorm_TFit)];
   writecell(DataNorm_TFit_H, T_path, 'Delimiter', '\t');

   disp(['Scaled data, Te, Tph, and Tm are saved as T_', sample_name, ' in ', fold_name,'.'])
end

%% sensplots
if sens_th
    [Sens,Sensnorm] = ThreeT_CeLeTdep_Pulse_SENS_v3(dx,nx,Lambda1,Lambda2,Lambda3,...
        C1,gamma2,C3,T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,...
        PulseTemp,IndStart,IndEnd,IndFit);
end

if sens_z
    opd_sens = 15e-9;   % assuming a single exponential decay form for absorption profile
    opd_top = 1;        % 1 for pump incident on the top surface; 0 for pump on substrate
    PulseFWHM = 1e-12;  % assuming a Gaussian for temporal profile of the laser pulse
    
    [Sens_z,Sensnorm_z] = ThreeT_CeLeTdep_Pulse_SENS_z_v3(dx,nx,Lambda1,Lambda2,Lambda3,...
        C1,gamma2,C3,T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,...
        PulseTemp,IndStart,IndEnd,IndFit,h,opd_sens,opd_top,PulseFWHM,para); 
end