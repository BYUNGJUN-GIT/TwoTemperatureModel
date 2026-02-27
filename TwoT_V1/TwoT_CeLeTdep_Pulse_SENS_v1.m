function [Sens,Sensnorm] = TwoT_CeLeTdep_Pulse_SENS_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,...
    T0,G,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp,IndFit)
% Sensitivity calculator for the 2-temperature model.
% Sensitivity definition (logarithmic):
%   S_p = dlog(TE)/dlog(p) ~= [log(TE_perturbed)-log(TE_ref)]/log(1.01)

[~,TE_ref] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,G,gEL,...
    tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);

TE_fit = TE_ref(:,IndFit);
TEmax = max(TE_fit);

S_gEL = zeros(length(tdelay_model),length(gEL));
for j = 1:length(gEL)
    gEL_temp = gEL;
    if gEL(j) > 1
        gEL_temp(j) = gEL_temp(j)*1.01;
        [~,TE_temp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,G,gEL_temp,...
            tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        S_gEL(:,j) = log_ratio(TE_temp(:,IndFit), TE_fit)/log(1.01);
    end
end

S_LambdaL = zeros(length(tdelay_model),length(LambdaL));
for j = 1:length(LambdaL)
    LambdaL_temp = LambdaL;
    if LambdaL(j) > 0.01
        LambdaL_temp(j) = LambdaL_temp(j)*1.01;
        [~,TE_temp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL_temp,LambdaE,CL,gammaE,T0,G,gEL,...
            tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        S_LambdaL(:,j) = log_ratio(TE_temp(:,IndFit), TE_fit)/log(1.01);
    end
end

S_LambdaE = zeros(length(tdelay_model),length(LambdaE));
for j = 1:length(LambdaE)
    LambdaE_temp = LambdaE;
    if LambdaE(j) > 0.01
        LambdaE_temp(j) = LambdaE_temp(j)*1.01;
        [~,TE_temp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE_temp,CL,gammaE,T0,G,gEL,...
            tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        S_LambdaE(:,j) = log_ratio(TE_temp(:,IndFit), TE_fit)/log(1.01);
    end
end

S_CL = zeros(length(tdelay_model),length(CL));
for j = 1:length(CL)
    CL_temp = CL;
    if CL(j) > 1
        CL_temp(j) = CL_temp(j)*1.01;
        [~,TE_temp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL_temp,gammaE,T0,G,gEL,...
            tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        S_CL(:,j) = log_ratio(TE_temp(:,IndFit), TE_fit)/log(1.01);
    end
end

S_gammaE = zeros(length(tdelay_model),length(gammaE));
for j = 1:length(gammaE)
    gammaE_temp = gammaE;
    if gammaE(j) > 0.01
        gammaE_temp(j) = gammaE_temp(j)*1.01;
        [~,TE_temp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE_temp,T0,G,gEL,...
            tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        S_gammaE(:,j) = log_ratio(TE_temp(:,IndFit), TE_fit)/log(1.01);
    end
end

nIf = size(G,1);
S_G11 = zeros(length(tdelay_model),nIf);
S_G22 = zeros(length(tdelay_model),nIf);
for j = 1:nIf
    G_temp = G;
    if G_temp(j,1) > 0
        G_temp(j,1) = G_temp(j,1)*1.01;
        [~,TE_temp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,G_temp,gEL,...
            tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        S_G11(:,j) = log_ratio(TE_temp(:,IndFit), TE_fit)/log(1.01);
    end

    G_temp = G;
    if G_temp(j,4) > 0
        G_temp(j,4) = G_temp(j,4)*1.01;
        [~,TE_temp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,G_temp,gEL,...
            tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        S_G22(:,j) = log_ratio(TE_temp(:,IndFit), TE_fit)/log(1.01);
    end
end

Sens.tdelay_model = 1e12*tdelay_model(:);
Sens.S_gEL = S_gEL;
Sens.S_LambdaL = S_LambdaL;
Sens.S_LambdaE = S_LambdaE;
Sens.S_CL = S_CL;
Sens.S_gammaE = S_gammaE;
Sens.S_G11 = S_G11;
Sens.S_G22 = S_G22;

scale = TE_fit./max(TEmax,eps);
Sensnorm.tdelay_model = 1e12*tdelay_model(:);
Sensnorm.S2_gEL = S_gEL.*scale;
Sensnorm.S2_LambdaL = S_LambdaL.*scale;
Sensnorm.S2_LambdaE = S_LambdaE.*scale;
Sensnorm.S2_CL = S_CL.*scale;
Sensnorm.S2_gammaE = S_gammaE.*scale;
Sensnorm.S2_G11 = S_G11.*scale;
Sensnorm.S2_G22 = S_G22.*scale;

end

function out = log_ratio(num,den)
num_safe = max(num, eps);
den_safe = max(den, eps);
out = log(num_safe) - log(den_safe);
end
