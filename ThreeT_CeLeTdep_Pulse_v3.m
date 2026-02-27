function [T1,T2,T3] = ThreeT_CeLeTdep_Pulse_v3(dx,nx,Lambda1,Lambda2,Lambda3,...
    C1,gamma2,C3,T0,G,g12,g13,g23,tdelay_model,delta_time,nt,heat1,heat2,heat3,...
    PulseTemp,IndStart,IndEnd)
%% initialize parameters
Tnew = zeros(3*sum(nx),1);   %placeholder for the three temperature profiles at a certain time; T(1:3:end) is Tp, T(2:3:end) is Te; T(3:3:end) is Tm
Told = Tnew;
q_old = Tnew;
q_new = Tnew;
T1 = zeros(nt+1,sum(nx));
T2 = T1;
T3 = T1;
a_matrix = zeros(3*sum(nx),3*sum(nx));
b_matrix = a_matrix;

tic

%% compute matrix elements of heat diffusion model
theta = 0.5;
currentTime = tdelay_model(1);
C2 = zeros(sum(nx),1);
Ke = zeros(sum(nx),1);
for iL = 1:length(nx)
    node_start = sum(nx(1:iL-1))+1;
    node_end = sum(nx(1:iL));
    I = node_start:node_end;
    C2(I) = gamma2(iL)*T0*ones(nx(iL),1);
    Ke(I) = Lambda2(iL);
end
Tnew2 = Told(3*IndStart-1:3:3*IndEnd);
for m = 1:nt
    Z = 1;
    while Z > 1e-5
        dt = delta_time(m);

        %------------------matrix elements for surface layer-------------------     
        iL = 1;
        node_start = sum(nx(1:iL-1))+1;
        node_end = sum(nx(1:iL));
        I = node_start:node_end;
        K1 = Lambda1(iL)*ones(length(I),1);
        % K2 = Lambda2(iL)*ones(length(I),1);
        K2 = Ke(I);
        K3 = Lambda3(iL)*ones(length(I),1);
        HC1 = C1(iL)*ones(length(I),1);
        HC2 = C2(I);
        HC3 = C3(iL)*ones(length(I),1);

        %left boundary matrix elements (adiabatic)
        [row1,row2,row3] = pop_row_adb_left_ThreeT(K1(1), K2(1), K3(1), HC1(1), HC2(1), HC3(1), g12(iL), g13(iL), g23(iL), dt ,dx(iL), theta-1);
        [ROW1,ROW2,ROW3] = pop_row_adb_left_ThreeT(K1(1), K2(1), K3(1), HC1(1), HC2(1), HC3(1), g12(iL), g13(iL), g23(iL), dt ,dx(iL), theta);
        node = node_start;
        a_matrix(node,node:node+5) = row1;
        a_matrix(node+1,node:node+5) = row2;
        a_matrix(node+2,node:node+5) = row3;       
        b_matrix(node,node:node+5) = ROW1;
        b_matrix(node+1,node:node+5) = ROW2;
        b_matrix(node+2,node:node+5) = ROW3;

        %interior nodes matrix elements   
        clear row1 row2 row3 ROW1 ROW2 ROW3
        [row1,row2,row3] = pop_row_int_vec_ThreeT(K1(2:end-1), K2(2:end-1), K3(2:end-1), HC1(2:end-1), HC2(2:end-1), HC3(2:end-1), g12(1), g13(1), g23(1), dt ,dx(1), theta-1);
        [ROW1,ROW2,ROW3] = pop_row_int_vec_ThreeT(K1(2:end-1), K2(2:end-1), K3(2:end-1), HC1(2:end-1), HC2(2:end-1), HC3(2:end-1), g12(1), g13(1), g23(1), dt ,dx(1), theta);
        index = 1;
        for node = node_start+1:node_end-1
            a_matrix((node*3-2),(node*3-5):(node*3+3)) = row1(index,:);
            a_matrix((node*3-1),(node*3-5):(node*3+3)) = row2(index,:);
            a_matrix((node*3),(node*3-5):(node*3+3)) = row3(index,:);
            b_matrix((node*3-2),(node*3-5):(node*3+3)) = ROW1(index,:);
            b_matrix((node*3-1),(node*3-5):(node*3+3)) = ROW2(index,:);
            b_matrix((node*3),(node*3-5):(node*3+3)) = ROW3(index,:);
            index = index+1;        
        end

        %right boundary matrix elements (conductance)
        clear row1 row2 row3 ROW1 ROW2 ROW3    
        [row1,row2,row3] = pop_row_cond_right_ThreeT(K1(end), K2(end), K3(end), HC1(end), HC2(end), HC3(end), g12(1), g13(1), g23(1),G(1,:), dt ,dx(1), theta-1); 
        [ROW1,ROW2,ROW3] = pop_row_cond_right_ThreeT(K1(end), K2(end), K3(end), HC1(end), HC2(end), HC3(end), g12(1), g13(1), g23(1),G(1,:), dt ,dx(1), theta);
        node = node_end;
        a_matrix((node*3-2),(node*3-5):(node*3+3)) = row1;
        a_matrix((node*3-1),(node*3-5):(node*3+3)) = row2;
        a_matrix((node*3),(node*3-5):(node*3+3)) = row3;     
        b_matrix((node*3-2),(node*3-5):(node*3+3)) = ROW1;
        b_matrix((node*3-1),(node*3-5):(node*3+3)) = ROW2;
        b_matrix((node*3),(node*3-5):(node*3+3)) = ROW3;

        %------------------matrix elements for inner layers--------------------        
        for iL = 2:(length(nx)-1)
            node_start = sum(nx(1:iL-1))+1;
            node_end = sum(nx(1:iL));
            I = node_start:node_end;
            K1 = Lambda1(iL)*ones(length(I),1);
            % K2 = Lambda2(iL)*ones(length(I),1);
            K2 = Ke(I);
            K3 = Lambda3(iL)*ones(length(I),1);
            HC1 = C1(iL)*ones(length(I),1);
            HC2 = C2(I);
            HC3 = C3(iL)*ones(length(I),1);

            %left boundary matrix elements (conductance)   
            clear row1 row2 row3 ROW1 ROW2 ROW3        
            [row1,row2,row3] = pop_row_cond_left_ThreeT(K1(1), K2(1), K3(3), HC1(1), HC2(1), HC3(1), g12(iL), g13(iL), g23(iL),G(iL-1,:), dt ,dx(iL), theta-1);
            [ROW1,ROW2,ROW3] = pop_row_cond_left_ThreeT(K1(1), K2(1), K3(3), HC1(1), HC2(1), HC3(1), g12(iL), g13(iL), g23(iL),G(iL-1,:), dt ,dx(iL), theta);        
            node = node_start;
            a_matrix((node*3-2),(node*3-5):(node*3+3)) = row1;
            a_matrix((node*3-1),(node*3-5):(node*3+3)) = row2;
            a_matrix((node*3),(node*3-5):(node*3+3)) = row3;     
            b_matrix((node*3-2),(node*3-5):(node*3+3)) = ROW1;
            b_matrix((node*3-1),(node*3-5):(node*3+3)) = ROW2;
            b_matrix((node*3),(node*3-5):(node*3+3)) = ROW3;

            %interior nodes matrix elements
            clear row1 row2 row3 ROW1 ROW2 ROW3
            [row1,row2, row3]=pop_row_int_vec_ThreeT(K1(2:end-1), K2(2:end-1), K3(2:end-1), HC1(2:end-1), HC2(2:end-1), HC3(2:end-1), g12(iL), g13(iL), g23(iL), dt ,dx(iL), theta-1);
            [ROW1,ROW2, ROW3]=pop_row_int_vec_ThreeT(K1(2:end-1), K2(2:end-1), K3(2:end-1), HC1(2:end-1), HC2(2:end-1), HC3(2:end-1),  g12(iL), g13(iL), g23(iL), dt ,dx(iL), theta);
            index = 1;
            for node = node_start+1:node_end-1
                a_matrix((node*3-2),(node*3-5):(node*3+3)) = row1(index,:);
                a_matrix((node*3-1),(node*3-5):(node*3+3)) = row2(index,:);
                a_matrix((node*3),(node*3-5):(node*3+3)) = row3(index,:);
                b_matrix((node*3-2),(node*3-5):(node*3+3)) = ROW1(index,:);
                b_matrix((node*3-1),(node*3-5):(node*3+3)) = ROW2(index,:);
                b_matrix((node*3),(node*3-5):(node*3+3)) = ROW3(index,:);
                index = index+1;        
            end

            %right boundary matrix elements (conductance)
            clear row1 row2 row3 ROW1 ROW2 ROW3
            [row1,row2,row3] = pop_row_cond_right_ThreeT(K1(end), K2(end), K3(end), HC1(end), HC2(end), HC3(end), g12(iL), g13(iL), g23(iL),G(iL,:), dt ,dx(iL), theta-1);
            [ROW1,ROW2,ROW3] = pop_row_cond_right_ThreeT(K1(end), K2(end), K3(end), HC1(end), HC2(end), HC3(end), g12(iL), g13(iL), g23(iL),G(iL,:), dt ,dx(iL), theta);
            node = node_end;
            a_matrix((node*3-2),(node*3-5):(node*3+3)) = row1;
            a_matrix((node*3-1),(node*3-5):(node*3+3)) = row2;
            a_matrix((node*3),(node*3-5):(node*3+3)) = row3;     
            b_matrix((node*3-2),(node*3-5):(node*3+3)) = ROW1;
            b_matrix((node*3-1),(node*3-5):(node*3+3)) = ROW2;
            b_matrix((node*3),(node*3-5):(node*3+3)) = ROW3;      
        end

        %-----------------matrix elements for bottom layer---------------------
        iL = length(nx);
        node_start = sum(nx(1:iL-1))+1;
        node_end = sum(nx(1:iL));
        I = node_start:node_end;
        K1 = Lambda1(iL)*ones(length(I),1);
        % K2 = Lambda2(iL)*ones(length(I),1);
        K2 = Ke(I);
        K3 = Lambda3(iL)*ones(length(I),1);
        HC1 = C1(iL)*ones(length(I),1);
        HC2 = C2(I);
        HC3 = C3(iL)*ones(length(I),1);

        %left boundary matrix elements (conductance)   
        clear row1 row2 row3 ROW1 ROW2 ROW3
        [row1,row2,row3]=pop_row_cond_left_ThreeT(K1(1), K2(1), K3(1), HC1(1), HC2(1), HC3(1), g12(end), g13(end), g23(end),G(length(nx)-1,:), dt ,dx(iL), theta-1);    
        [ROW1,ROW2,ROW3]=pop_row_cond_left_ThreeT(K1(1), K2(1), K3(1), HC1(1), HC2(1), HC3(1), g12(end), g13(end), g23(end),G(length(nx)-1,:), dt ,dx(iL), theta);    
        node = node_start;
        a_matrix((node*3-2),(node*3-5):(node*3+3)) = row1;
        a_matrix((node*3-1),(node*3-5):(node*3+3)) = row2;
        a_matrix((node*3),(node*3-5):(node*3+3)) = row3;     
        b_matrix((node*3-2),(node*3-5):(node*3+3)) = ROW1;
        b_matrix((node*3-1),(node*3-5):(node*3+3)) = ROW2;
        b_matrix((node*3),(node*3-5):(node*3+3)) = ROW3;

        %interior nodes matrix elements   
        clear row1 row2 row3 ROW1 ROW2 ROW3
        [row1,row2, row3]=pop_row_int_vec_ThreeT(K1(2:end-1), K2(2:end-1), K3(2:end-1), HC1(2:end-1), HC2(2:end-1), HC3(2:end-1), g12(iL), g13(iL), g23(iL), dt ,dx(iL), theta-1);
        [ROW1,ROW2, ROW3]=pop_row_int_vec_ThreeT(K1(2:end-1), K2(2:end-1), K3(2:end-1), HC1(2:end-1), HC2(2:end-1), HC3(2:end-1),  g12(iL), g13(iL), g23(iL), dt ,dx(iL), theta);
        index = 1;
        for node = node_start+1:node_end-1
            a_matrix((node*3-2),(node*3-5):(node*3+3)) = row1(index,:);
            a_matrix((node*3-1),(node*3-5):(node*3+3)) = row2(index,:);
            a_matrix((node*3),(node*3-5):(node*3+3)) = row3(index,:);
            b_matrix((node*3-2),(node*3-5):(node*3+3)) = ROW1(index,:);
            b_matrix((node*3-1),(node*3-5):(node*3+3)) = ROW2(index,:);
            b_matrix((node*3),(node*3-5):(node*3+3)) = ROW3(index,:);
            index = index+1;        
        end

        %right boundary matrix element (isothermal)
        clear row1 row2 row3 ROW1 ROW2 ROW3
        [row1,row2,row3] = pop_row_iso_right_ThreeT(K1(end), K2(end), K3(end), HC1(end), HC2(end), HC3(end), g12(end), g13(end), g23(end), dt ,dx(iL), theta-1);
        [ROW1,ROW2,ROW3] = pop_row_iso_right_ThreeT(K1(end), K2(end), K3(end), HC1(end), HC2(end), HC3(end), g12(end), g13(end), g23(end), dt ,dx(iL), theta);
        node = node_end;
        a_matrix((node*3-2),(node*3-5):(node*3)) = row1;
        a_matrix((node*3-1),(node*3-5):(node*3)) = row2;
        a_matrix((node*3),(node*3-5):(node*3)) = row3;     
        b_matrix((node*3-2),(node*3-5):(node*3)) = ROW1;
        b_matrix((node*3-1),(node*3-5):(node*3)) = ROW2;
        b_matrix((node*3),(node*3-5):(node*3)) = ROW3;        

        %----------------------------source term-------------------------------
        node = 1;  
        for iL = 1:length(nx)-1       
            for n = 1:nx(iL)     

                q_new(node*3-2)= dt/C1(iL)*heat1(node)*PulseTemp(m);
                q_new(node*3-1)= dt/C2(node)*heat2(node)*PulseTemp(m);
                q_new(node*3)= dt/C3(iL)*heat3(node)*PulseTemp(m);
                
                node = node+1;
            end
        end        
        Amatrix = sparse(a_matrix);
        Bmatrix = b_matrix*Told+theta*q_old+(1-theta)*q_new;
        Tnew = Amatrix\Bmatrix; %linsolve(Amatrix,Bmatrix)   

        %----------------recalculate electronic heat capacity--------------
        for iL = 1:length(nx)
            node_start = sum(nx(1:iL-1))+1;
            node_end = sum(nx(1:iL));
            I = node_start:node_end;
            C2(I) = gamma2(iL)*0.5*(T0 + Tnew(3*node_start-1:3:3*node_end) + T0 + Told(3*node_start-1:3:3*node_end));
            Ke(I) = Lambda2(iL)/T0*0.5*(T0 + Tnew(3*node_start-1:3:3*node_end) + T0 + Told(3*node_start-1:3:3*node_end));
        end   

        res = (Tnew2-Tnew(3*IndStart-1:3:3*IndEnd)).^2./((IndEnd-IndStart)*Tnew(3*IndStart-1:3:3*IndEnd).^2);
        Z = sqrt(sum(res));
        
        Tnew2 = Tnew(3*IndStart-1:3:3*IndEnd);

    end

    Told = Tnew;
    q_old = q_new;
    currentTime = currentTime+dt;

    T1(m+1,:) = Tnew(1:3:end);
    T2(m+1,:) = Tnew(2:3:end);  
    T3(m+1,:) = Tnew(3:3:end); 
end

toc
