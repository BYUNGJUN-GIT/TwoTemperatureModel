function [h,dx,nx,z0,heatL,heatE,PulseTemp,IndFit] = TwoT_Para_z_v1(para,h,opd,opd_top,importPulse,PulseFWHM)
% Build geometry/source/pulse variants for z-related sensitivities

dx = para.dx;
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

heatL = zeros(sum(nx),1);
heatE = zeros(sum(nx),1);
iheat = para.iheat;
I_heat = 1:sum(nx(1:iheat));
zAbs = z0(I_heat);

if opd_top == 1
    heatE(I_heat) = exp(-zAbs/opd);
else
    heatE(I_heat) = exp((zAbs-sum(h(1:iheat)))/opd);
end

PulseTemp = zeros(para.nt,1);
if importPulse == 1
    PulseTemp = para.PulseTemp;
else
    sigma_corr = PulseFWHM/2/sqrt(2*log(2));
    for i = 1:para.nt
        tt = para.tdelay_model(i+1);
        PulseTemp(i) = exp(-tt^2/(2*sigma_corr^2));
    end
end

F2 = sum(dx(1)*heatE(I_heat)) * sum(PulseTemp.*para.delta_time(:));
heatE(I_heat) = heatE(I_heat)*para.Fluence/F2;

IndFit = max(1, round(nx(1)/2));
end
