function [Sens,Sensnorm] = ThreeT_CeLeTdep_Pulse_SENS_z_v3(dx,nx,Lambda1,Lambda2,Lambda3,...
    C1,gamma2,C3,T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,...
    PulseTemp,IndStart,IndEnd,IndFit,h,opd_sens,opd_top,PulseFWHM,para)
%% SENS to optical penetration depth 
importPulse = 1;
S_opd = zeros(length(tdelay_model),1);

[~,~,~,~,heat1A,heat2A,heat3A,~] = ThreeT_Para_z_v3(para,h,opd_sens,opd_top,importPulse,PulseFWHM);
[~,~,~,~,heat1B,heat2B,heat3B,~] = ThreeT_Para_z_v3(para,h,opd_sens*1.01,opd_top,importPulse,PulseFWHM);

[~,~,T3_opd_ref] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
    T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1A,heat2A,heat3A,PulseTemp,IndStart,IndEnd);  

[~,~,T3_opd_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
    T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1B,heat2B,heat3B,PulseTemp,IndStart,IndEnd);  

Num = log(T3_opd_temp(:,IndFit)) - log(T3_opd_ref(:,IndFit));
S_opd(:,1) = Num./log(1.01);
  
%% SENS to pulse width 
importPulse = 0;
S_pulse = zeros(length(tdelay_model),1);

[~,~,~,~,~,~,~,Pulsetemp_ref] = ThreeT_Para_z_v3(para,h,opd_sens,opd_top,importPulse,PulseFWHM);
[~,~,~,~,~,~,~,Pulsetemp_temp] = ThreeT_Para_z_v3(para,h,opd_sens,opd_top,importPulse,PulseFWHM*1.01);

[~,~,T3_pulse_ref] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
    T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,Pulsetemp_ref,IndStart,IndEnd);  

[~,~,T3_pulse_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
    T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,Pulsetemp_temp,IndStart,IndEnd);  

Num = log(T3_pulse_temp(:,IndFit)) - log(T3_pulse_ref(:,IndFit));
S_pulse(:,1) = Num./log(1.01);

%% SENS to thickness 
importPulse = 1;
S_h = zeros(length(tdelay_model),length(h));

for j = 1:length(h)
    h_temp = h;
    h_temp(j) = h(j)+dx(j);
    
    [~,~,nx_ref,~,heat1A,heat2A,heat3A,~] = ThreeT_Para_z_v3(para,h,opd_sens,opd_top,importPulse,PulseFWHM);      
    [~,~,nx_temp,~,heat1B,heat2B,heat3B,~] = ThreeT_Para_z_v3(para,h_temp,opd_sens,opd_top,importPulse,PulseFWHM);
   
    [~,~,T3_h_ref] = ThreeT_CeLeTdep_Pulse_v3(dx,nx_ref,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
        T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1A,heat2A,heat3A,PulseTemp,IndStart,IndEnd);  

    [~,~,T3_h_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx_temp,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
        T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1B,heat2B,heat3B,PulseTemp,IndStart,IndEnd+1);  
    
    if j <= 2
        IndFit = IndFit + 1;
    end
    
    Num = log(T3_h_temp(:,IndFit)) - log(T3_h_ref(:,IndFit));
    Denom = log(h_temp(j)) - log(h(j));
    S_h(:,j) = Num./Denom;
end

%% Return Sensitivity
Sens.tdelay_model = 1e12*tdelay_model';
Sens.S_opd = S_opd;
Sens.S_pulse = S_pulse;
Sens.S_h = S_h;

Sensnorm.tdelay_model = 1e12*tdelay_model';
Sensnorm.S2_opd = S_opd.*T3_opd_ref(:,IndFit)./max(T3_opd_ref(:,IndFit));
Sensnorm.S2_pulse = S_pulse.*T3_pulse_ref(:,IndFit)./max(T3_pulse_ref(:,IndFit));
Sensnorm.S2_h = S_h.*T3_h_ref(:,IndFit)./max(T3_h_ref(:,IndFit));

end