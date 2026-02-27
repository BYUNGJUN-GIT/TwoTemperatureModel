function [Sens,Sensnorm] = TwoT_CeLeTdep_Pulse_SENS_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp,IndFit)
% Sensitivity of TTM response (TE at IndFit) to thermal parameters

[~,TE_ref] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
TE_trace = TE_ref(:,IndFit);
TE_max = max(TE_trace);

scale = 1.01;
denom = log(scale);

S_LambdaL = calc_vec_sens(LambdaL, 1e-9, @(x) run_te(x,LambdaE,CL,gammaE,Gll,gEL));
S_LambdaE = calc_vec_sens(LambdaE, 1e-9, @(x) run_te(LambdaL,x,CL,gammaE,Gll,gEL));
S_CL      = calc_vec_sens(CL,      1e-9, @(x) run_te(LambdaL,LambdaE,x,gammaE,Gll,gEL));
S_gammaE  = calc_vec_sens(gammaE,  1e-12,@(x) run_te(LambdaL,LambdaE,CL,x,Gll,gEL));
S_gEL     = calc_vec_sens(gEL,     1e-9, @(x) run_te(LambdaL,LambdaE,CL,gammaE,Gll,x));
S_Gll     = calc_vec_sens(Gll(:)', 1e-9, @(x) run_te(LambdaL,LambdaE,CL,gammaE,x(:),gEL));

Sens.tdelay_model = 1e12*tdelay_model(:);
Sens.S_LambdaL = S_LambdaL;
Sens.S_LambdaE = S_LambdaE;
Sens.S_CL = S_CL;
Sens.S_gammaE = S_gammaE;
Sens.S_gEL = S_gEL;
Sens.S_Gll = S_Gll;

normFactor = TE_trace./max(TE_max, eps);
Sensnorm.tdelay_model = Sens.tdelay_model;
Sensnorm.S2_LambdaL = S_LambdaL.*normFactor;
Sensnorm.S2_LambdaE = S_LambdaE.*normFactor;
Sensnorm.S2_CL = S_CL.*normFactor;
Sensnorm.S2_gammaE = S_gammaE.*normFactor;
Sensnorm.S2_gEL = S_gEL.*normFactor;
Sensnorm.S2_Gll = S_Gll.*normFactor;

    function Sout = calc_vec_sens(vec,thresh,run_fn)
        Sout = zeros(length(tdelay_model), numel(vec));
        for j = 1:numel(vec)
            if vec(j) > thresh
                vtemp = vec;
                vtemp(j) = vtemp(j)*scale;
                TE_temp = run_fn(vtemp);
                Sout(:,j) = (log(max(TE_temp,eps)) - log(max(TE_trace,eps)))./denom;
            end
        end
    end

    function TE_col = run_te(L1,L2,C1,g2,G1,g1)
        [~,TE_tmp] = TwoT_CeLeTdep_Pulse_v1(dx,nx,L1,L2,C1,g2,T0,G1,g1,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp);
        TE_col = TE_tmp(:,IndFit);
    end
end
