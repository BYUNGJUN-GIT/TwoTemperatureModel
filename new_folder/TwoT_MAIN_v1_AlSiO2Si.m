%% Two-temperature model package for Al/SiO2/Si
% T_L = lattice (phonon) temperature
% T_E = electron temperature

[h,dx,nx,z0,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp,Data,fold_name,sample_name,para] = TwoT_Para_v1_AlSiO2Si();

para_path = fullfile(pwd, fold_name, ['Para_', sample_name, '.mat']);
T_path = fullfile(pwd, fold_name, ['T_', sample_name, '.csv']);

sim = 1;
importData = size(Data,1) > 4;
modVin = 0; % set 1 for offset correction at negative delay

% default readout point: center of Al
IndFit = max(1, round(nx(1)/2));

if sim
    [TL,TE] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
    save(para_path);

    T_Fit = [1e12*tdelay_model(:), TL(:,IndFit), TE(:,IndFit)];
end

if importData
    tdelay_start = -5;
    tdelay_end = 10;
    tdelay_norm = 5;

    if modVin == 1
        negtd_ind = Data(:,2) < -1;
        negtd_offset = mean(Data(negtd_ind, 3));
        Data(:,3) = Data(:,3)-negtd_offset;
    end

    [~,t_start] = min(abs(Data(:,2) - tdelay_start));
    [~,t_end]   = min(abs(Data(:,2) - tdelay_end));
    [~,t_data]  = min(abs(Data(:,2) - tdelay_norm));
    [~,t_model] = min(abs(tdelay_model - tdelay_norm*1e-12));

    Data_norm = Data(t_start:t_end,3)./Data(t_data,3)*TE(t_model,IndFit);
    Data_norm = [Data(t_start:t_end,2), Data_norm];

    TL_model = interp1(T_Fit(:,1),T_Fit(:,2),Data_norm(:,1));
    TE_model = interp1(T_Fit(:,1),T_Fit(:,3),Data_norm(:,1));
    res = (TE_model - Data_norm(:,2)).^2;
    Z = sum(res)/length(Data_norm(:,2));
    disp(['Goodness-of-fit Z = ', num2str(Z)])

    figure(1)
    plot(Data_norm(:,1),Data_norm(:,2),'ok',Data_norm(:,1),TE_model,'-r',Data_norm(:,1),TL_model,'-b')
    xlabel('Time delay (ps)')
    ylabel('\DeltaT (K)')
    xlim([tdelay_start tdelay_end])
    ylim([0 inf])
    legend('Data', 'T_E (electrons)', 'T_L (lattice)');

    figure(2)
    loglog(Data(:,2),Data(:,3)./Data(t_data,3)*TE(t_model,IndFit),'ok',...
        tdelay_model*1e12,TE(:,IndFit),'-r', tdelay_model*1e12, TL(:, IndFit),'-b')
    xlabel('Time delay (ps)')
    ylabel('\DeltaT (K)')
    xlim([0.1 inf])
    ylim([min(TL(:, IndFit)), inf])
    legend('Data', 'T_E (electrons)', 'T_L (lattice)');

    header = {'t_d_data_ps', 'Normalized Data', 't_d_model_ps', 'DeltaT_lattice_K', 'DeltaT_electron_K'};
    Data_norm = [Data_norm; zeros(size(T_Fit,1)-size(Data_norm,1), size(Data_norm,2))];
    DataNorm_TFit = [Data_norm, T_Fit];
    DataNorm_TFit_H = [header; num2cell(DataNorm_TFit)];
    writecell(DataNorm_TFit_H, T_path, 'Delimiter', '\t');

    disp(['Scaled data + TL/TE saved as T_', sample_name, ' in ', fold_name,'.'])
else
    header = {'t_d_model_ps','DeltaT_lattice_K','DeltaT_electron_K'};
    out = [header; num2cell(T_Fit)];
    writecell(out, T_path, 'Delimiter', '\t');

    figure(1);
    plot(tdelay_model*1e12, TE(:,IndFit), '-r', tdelay_model*1e12, TL(:,IndFit), '-b', 'LineWidth', 1.2);
    xlabel('Time delay (ps)');
    ylabel('\DeltaT (K)');
    legend('T_E (electrons)','T_L (lattice)','Location','best');
    xlim([min(tdelay_model)*1e12, 10]);
    grid on;

    figure(2);
    semilogx(tdelay_model(tdelay_model>0)*1e12, TE(tdelay_model>0,IndFit), '-r', ...
             tdelay_model(tdelay_model>0)*1e12, TL(tdelay_model>0,IndFit), '-b', 'LineWidth', 1.2);
    xlabel('Time delay (ps)');
    ylabel('\DeltaT (K)');
    legend('T_E (electrons)','T_L (lattice)','Location','best');
    grid on;

    disp(['TTM done. Outputs saved to ', T_path]);
end
