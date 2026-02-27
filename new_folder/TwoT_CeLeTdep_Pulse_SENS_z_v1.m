function [Sens,Sensnorm] = TwoT_CeLeTdep_Pulse_SENS_z_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,para)
% Sensitivity to optical penetration depth, pulse width, and thickness

scale = 1.01;

if isfield(para,'opd')
    opd_ref = para.opd;
else
    opd_ref = 12e-9;
end
if isfield(para,'PulseFWHM')
    pulse_ref = para.PulseFWHM;
else
    pulse_ref = 0.6e-12;
end

opd_top = 1;

[~,~,nxA,~,heatL_A,heatE_A,PulseA,IndFit_A] = TwoT_Para_z_v1(para,para.h,opd_ref,opd_top,1,pulse_ref);
[~,TE_opd_ref] = TwoT_CeLeTdep_Pulse_v1(dx,nxA,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL_A,heatE_A,PulseA);

[~,~,nxB,~,heatL_B,heatE_B,PulseB,IndFit_B] = TwoT_Para_z_v1(para,para.h,opd_ref*scale,opd_top,1,pulse_ref);
[~,TE_opd_tmp] = TwoT_CeLeTdep_Pulse_v1(dx,nxB,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL_B,heatE_B,PulseB);

S_opd = (log(max(TE_opd_tmp(:,IndFit_B),eps))-log(max(TE_opd_ref(:,IndFit_A),eps)))./log(scale);

[~,~,nxC,~,heatL_C,heatE_C,PulseC,IndFit_C] = TwoT_Para_z_v1(para,para.h,opd_ref,opd_top,0,pulse_ref);
[~,TE_pw_ref] = TwoT_CeLeTdep_Pulse_v1(dx,nxC,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL_C,heatE_C,PulseC);

[~,~,nxD,~,heatL_D,heatE_D,PulseD,IndFit_D] = TwoT_Para_z_v1(para,para.h,opd_ref,opd_top,0,pulse_ref*scale);
[~,TE_pw_tmp] = TwoT_CeLeTdep_Pulse_v1(dx,nxD,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL_D,heatE_D,PulseD);

S_pulse = (log(max(TE_pw_tmp(:,IndFit_D),eps))-log(max(TE_pw_ref(:,IndFit_C),eps)))./log(scale);

S_h = zeros(length(tdelay_model),length(para.h));
TE_h_ref = TE_opd_ref(:,IndFit_A);
for j = 1:length(para.h)
    h_tmp = para.h;
    h_tmp(j) = h_tmp(j) + para.dx(j);

    [~,~,nx_ref,~,heatL_ref,heatE_ref,Pulse_ref,IndFit_ref] = TwoT_Para_z_v1(para,para.h,opd_ref,opd_top,1,pulse_ref);
    [~,~,nx_tmp,~,heatL_tmp,heatE_tmp,Pulse_tmp,IndFit_tmp] = TwoT_Para_z_v1(para,h_tmp,opd_ref,opd_top,1,pulse_ref);

    [~,TE_ref] = TwoT_CeLeTdep_Pulse_v1(dx,nx_ref,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL_ref,heatE_ref,Pulse_ref);
    [~,TE_tmp] = TwoT_CeLeTdep_Pulse_v1(dx,nx_tmp,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL_tmp,heatE_tmp,Pulse_tmp);

    Num = log(max(TE_tmp(:,IndFit_tmp),eps)) - log(max(TE_ref(:,IndFit_ref),eps));
    Den = log(h_tmp(j)) - log(para.h(j));
    S_h(:,j) = Num./Den;
end

Sens.tdelay_model = 1e12*tdelay_model(:);
Sens.S_opd = S_opd;
Sens.S_pulse = S_pulse;
Sens.S_h = S_h;

Sensnorm.tdelay_model = Sens.tdelay_model;
Sensnorm.S2_opd = S_opd.*TE_h_ref./max(TE_h_ref);
Sensnorm.S2_pulse = S_pulse.*TE_h_ref./max(TE_h_ref);
Sensnorm.S2_h = S_h.*(TE_h_ref./max(TE_h_ref));
end
