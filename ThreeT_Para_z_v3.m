function [h,dx,nx,z0,heat1,heat2,heat3,PulseTemp]...
    = ThreeT_Para_z_v3(para,h_sens,opd_sens,opd_top,importPulse,PulseFWHM)
%% samle properties

h = h_sens;
dx = para.dx;
nx = round(h./dx);    %vector listing number of nodes of the layers

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

%% time discretization
tdelay_model = para.tdelay_model;
delta_time = para.delta_time;
nt = para.nt;

%% define heat absorption 
iheat = para.iheat;
I_heat = 1:sum(nx(1:iheat));
zAbs = z0(1:sum(nx(1:iheat)));

opd = opd_sens; %10e-9; 

heat1 = zeros(1,sum(nx));
heat2 = zeros(1,sum(nx));
heat3 = zeros(1,sum(nx));

if opd_top ==1
    heat2(I_heat) = exp(-zAbs/opd);
else 
    heat2(I_heat) = exp((zAbs-sum(h(1:iheat)))/opd); 
end

%% intensity profile of a pump pulse 

PulseTemp = zeros(1,nt);
F2 = 0;

pulseIndex = para.pulseIndex;

if importPulse == 1
    load('PulseWidth_20180224.mat','PulseWidth')
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
    sigma_pump = PulseFWHM/sqrt(2*log(2));
    for i = 2:pulseIndex
         tt = tdelay_model(i);   
         PulseTemp(i) = exp(-tt^2/sigma_pump^2);     
         F2 = F2 + PulseTemp(i)*delta_time(i-1);
    end
end

Fluence = para.Fluence;

F2 = sum(dx(1)*heat2(I_heat))*F2;
heat2(I_heat) = heat2(I_heat)*Fluence/F2;

end
