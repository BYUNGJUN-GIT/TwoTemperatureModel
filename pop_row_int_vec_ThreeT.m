function [ROW1,ROW2,ROW3] = pop_row_int_vec_ThreeT(lambda1, lambda2, lambda3, Cp1, Cp2, Cp3, g12, g13, g23, dt, dx, theta)

diff1 = lambda1./Cp1*dt/dx^2;
ch12 = dt*g12./Cp1;
ch13 = dt*g13./Cp1;
    
diff2 = lambda2./Cp2*dt/dx^2;
ch21 = dt*g12./Cp2;
ch23 = dt*g23./Cp2;
    
diff3 = lambda3./Cp3*dt/dx^2;
ch31 = dt*g13./Cp3;
ch32 = dt*g23./Cp3;

A1 = theta*diff1; %T1,i-1 to T1,i
B1 = zeros(size(A1)); %T2,i-1 to T1,i
C1 = zeros(size(A1)); %T3,i-1 to T1,i
    
D1 = 1-theta*2*diff1 - theta*ch12 - theta*ch13; %T1i to T1i
E1 = theta*ch12; %T1,i to T2,i
F1 = theta*ch13; %T1,i to T3,i
H1 = theta*diff1; %T1,i to T1,i+1
K1 = zeros(size(A1)); %T1,i to T2,i+1
L1 = zeros(size(A1)); %T1,i to T3,i+1
    
A2 = zeros(size(A1)); %T1,i-1 to T2,i
B2 = theta*diff2; %T2,i-1 to T2,i
C2 = zeros(size(A1)); %T3,i-1 to T2,i
D2 = theta*ch21; %T2i to T1i
E2 = 1-2*theta*diff2 - theta*ch21 - theta*ch23; %T2,i to T2,i
F2 = theta*ch23; %T2,i to T3,i
H2 = zeros(size(A1)); %T2,i to T1,i+1
K2 = theta*diff2; %T2,i to T2,i+1
L2 = zeros(size(A1)); %T2,i to T3,i+1
    
A3 = zeros(size(A1)); %T1,i-1 to T3,i
B3 = zeros(size(A1)); %T2,i-1 to T3,i
C3 = theta*diff3; %T3,i-1 to T3,i
D3 = theta*ch31; %T3i to T1i
E3 = theta*ch32; %T3,i to T2,i
F3 = 1-2*theta*diff3 - theta*ch31 - theta*ch32; %T3,i to T3,i
H3 = zeros(size(A1)); %T3,i to T1,i+1
K3 = zeros(size(A1)); %T3,i to T2,i+1
L3 = theta*diff3; %T3,i to T3,i+1


ROW1 = [A1 B1 C1 D1 E1 F1 H1 K1 L1];

ROW2 = [A2 B2 C2 D2 E2 F2 H2 K2 L2];
ROW3 = [A3 B3 C3 D3 E3 F3 H3 K3 L3];