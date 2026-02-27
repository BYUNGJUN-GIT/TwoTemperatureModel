function [Sens,Sensnorm] = ThreeT_CeLeTdep_Pulse_SENS_v3(dx,nx,Lambda1,Lambda2,Lambda3,...
    C1,gamma2,C3,T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,...
    PulseTemp,IndStart,IndEnd,IndFit)
%% Normal simulation 
[~,~,T3_ref] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
    T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
T3max = max(T3_ref(:,IndFit));

%% Sens to g12
clear T3_temp
S_g12 = zeros(length(tdelay_model),length(g12));

for j=1:length(g12) 
    g12_temp = g12;
    if g12(j) > 1
        g12_temp(j) = g12_temp(j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
            T0,G,g12_temp,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_g12(:,j) = Num./log(1.01);
    end    
end

%% Sens to g23
clear T3_temp
S_g23 = zeros(length(tdelay_model),length(g23));

for j=1:length(g23) 
    g23_temp = g23;
    if g23(j) > 1
        g23_temp(j) = g23_temp(j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
            T0,G,g12,g13,g23_temp,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_g23(:,j) = Num./log(1.01);
    end
end

%% Sens to Lambda1
clear T3_temp
S_Lambda1 = zeros(length(tdelay_model),length(Lambda1));

for j=1:length(Lambda1) 
    Lambda1_temp = Lambda1;
    if Lambda1(j) > 0.01
        Lambda1_temp(j) = Lambda1_temp(j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1_temp,Lambda2,Lambda3,C1,gamma2,C3,...
            T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_Lambda1(:,j) = Num./log(1.01);
    end    
end

%% Sens to Lambda2
clear T3_temp
S_Lambda2 = zeros(length(tdelay_model),length(Lambda2));

for j=1:length(Lambda2) 
    Lambda2_temp = Lambda2;
    if Lambda2(j) > 0.01
        Lambda2_temp(j) = Lambda2_temp(j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2_temp,Lambda3,C1,gamma2,C3,...
            T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_Lambda2(:,j) = Num./log(1.01);
    end    
end

%% Sens to C1
clear T3_temp
S_C1 = zeros(length(tdelay_model),length(C1));

for j=1:length(C1) 
    C1_temp = C1;
    if C1(j) > 1
        C1_temp(j) = C1_temp(j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1_temp,gamma2,C3,...
            T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_C1(:,j) = Num./log(1.01);
    end    
end

%% Sens to gamma2
clear T3_temp
S_gamma2 = zeros(length(tdelay_model),length(gamma2));

for j=1:length(gamma2) 
    gamma2_temp = gamma2;
    if gamma2(j) > 0.01
        gamma2_temp(j) = gamma2_temp(j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2_temp,C3,...
            T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_gamma2(:,j) = Num./log(1.01);
    end    
end

%% Sens to C3
clear T3_temp
S_C3 = zeros(length(tdelay_model),length(C3));

for j=1:length(C3) 
    C3_temp = C3;
    if C3(j) > 1
        C3_temp(j) = C3_temp(j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3_temp,...
            T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit))-log(T3_ref(:,IndFit));
        S_C3(:,j) = Num./log(1.01);
    end    
end

%% Sens to G_12
clear T3_temp
iLL = 1;
S_G12 = zeros(length(tdelay_model),length(G(iLL,:)));

for j=1:length(G(iLL,:)) 
    G_temp = G;
    if G(iLL,j) > 0
        G_temp(iLL,j)=G_temp(iLL,j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
            T0,G_temp,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_G12(:,j) = Num./log(1.01);
    end    
end

%% Sens to G_23
clear T3_temp
iLL = 2;
S_G23 = zeros(length(tdelay_model),length(G(iLL,:)));

for j=1:length(G(iLL,:)) 
    G_temp = G;
    if G(iLL,j) > 0
        G_temp(iLL,j)=G_temp(iLL,j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
            T0,G_temp,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
         Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
         S_G23(:,j) = Num./log(1.01);
    end    
end

%% Sens to G_34
clear T3_temp
iLL = 3;
S_G34 = zeros(length(tdelay_model),length(G(iLL,:)));

for j=1:length(G(iLL,:)) 
    G_temp = G;
    if G(iLL,j) > 0
        G_temp(iLL,j)=G_temp(iLL,j)*1.01;

        [~,~,T3_temp] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,C1,gamma2,C3,...
            T0,G_temp,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,PulseTemp,IndStart,IndEnd);  
        
        Num = log(T3_temp(:,IndFit)) - log(T3_ref(:,IndFit));
        S_G34(:,j) = Num./log(1.01);
    end    
end

%% Return Sensitivity
Sens.tdelay_model = 1e12*tdelay_model';
Sens.S_g12 = S_g12;
Sens.S_g23 = S_g23;
Sens.S_Lambda1 = S_Lambda1;
Sens.S_Lambda2 = S_Lambda2;
Sens.S_C1 = S_C1;
Sens.S_gamma2 = S_gamma2;
Sens.S_C3 = S_C3;
Sens.S_G12 = S_G12;
Sens.S_G23 = S_G23;
Sens.S_G34 = S_G34;

Sensnorm.tdelay_model = 1e12*tdelay_model';
Sensnorm.S2_g12 = S_g12.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_g23 = S_g23.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_Lambda1 = S_Lambda1.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_Lambda2 = S_Lambda2.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_C1 = S_C1.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_gamma2 = S_gamma2.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_C3 = S_C3.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_G12 = S_G12.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_G23 = S_G23.*T3_ref(:,IndFit)./T3max;
Sensnorm.S2_G34 = S_G34.*T3_ref(:,IndFit)./T3max;
end