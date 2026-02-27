function [TL,TE] = TwoT_CeLeTdep_Pulse_v1(dx,nx,LambdaL,LambdaE,CL,gammaE,T0,Gll,gEL,tdelay_model,delta_time,nt,heatL,heatE,PulseTemp)
% Two-temperature model solver (lattice + electron)

N = sum(nx);
TL_vec = T0*ones(N,1);
TE_vec = T0*ones(N,1);

TL = zeros(nt+1,N);
TE = zeros(nt+1,N);
TL(1,:) = TL_vec;
TE(1,:) = TE_vec;

layer_of_node = zeros(N,1);
node = 1;
for iL = 1:length(nx)
    layer_of_node(node:node+nx(iL)-1) = iL;
    node = node + nx(iL);
end

for m = 1:nt
    dt = delta_time(m);

    Ce = zeros(N,1);
    kL = zeros(N,1);
    kE = zeros(N,1);
    CLat = zeros(N,1);
    g = zeros(N,1);
    for i = 1:N
        iL = layer_of_node(i);
        Ce(i) = max(gammaE(iL)*TE_vec(i), gammaE(iL)*T0);
        kL(i) = LambdaL(iL);
        kE(i) = LambdaE(iL);
        CLat(i) = CL(iL);
        g(i) = gEL(iL);
    end

    AL = sparse(N,N);
    AE = sparse(N,N);
    BL = zeros(N,1);
    BE = zeros(N,1);

    for i = 1:N
        % left link
        if i == 1
            GL_left = 0; GE_left = 0; % adiabatic top
        else
            GL_left = link_cond(i-1,i,kL,dx,layer_of_node,Gll,true);
            GE_left = link_cond(i-1,i,kE,dx,layer_of_node,Gll,false);
        end

        % right link
        if i == N
            GL_right = 2*kL(i)/dx(layer_of_node(i)); % Dirichlet bath at T0
            GE_right = 2*kE(i)/dx(layer_of_node(i));
            T_right = T0;
        else
            GL_right = link_cond(i,i+1,kL,dx,layer_of_node,Gll,true);
            GE_right = link_cond(i,i+1,kE,dx,layer_of_node,Gll,false);
            T_right = 0;
        end

        vol = dx(layer_of_node(i)); % unit-area control volume

        AL(i,i) = CLat(i)/dt + (GL_left + GL_right)/vol + g(i);
        AE(i,i) = Ce(i)/dt + (GE_left + GE_right)/vol + g(i);

        if i > 1
            AL(i,i-1) = -GL_left/vol;
            AE(i,i-1) = -GE_left/vol;
        end
        if i < N
            AL(i,i+1) = -GL_right/vol;
            AE(i,i+1) = -GE_right/vol;
        end

        srcL = heatL(i)*PulseTemp(m);
        srcE = heatE(i)*PulseTemp(m);

        BL(i) = CLat(i)/dt*TL_vec(i) + g(i)*TE_vec(i) + srcL + GL_right/vol*T_right;
        BE(i) = Ce(i)/dt*TE_vec(i) + g(i)*TL_vec(i) + srcE + GE_right/vol*T_right;
    end

    % staggered implicit coupling for robustness
    TL_new = AL\BL;

    BE2 = BE + g.*(TL_new - TL_vec);
    TE_new = AE\BE2;

    TL_vec = TL_new;
    TE_vec = TE_new;

    TL(m+1,:) = TL_vec;
    TE(m+1,:) = TE_vec;
end

end

function G = link_cond(i,j,k,dx,layer_of_node,Gll,isLattice)
li = layer_of_node(i);
lj = layer_of_node(j);
ri = dx(li)/2;
rj = dx(lj)/2;

if li == lj
    G = 1/(ri/max(k(i),eps) + rj/max(k(j),eps));
else
    if isLattice
        gInt = Gll(min(li,lj));
        G = 1/(ri/max(k(i),eps) + 1/max(gInt,eps) + rj/max(k(j),eps));
    else
        G = 1/(ri/max(k(i),eps) + rj/max(k(j),eps));
    end
end
end
