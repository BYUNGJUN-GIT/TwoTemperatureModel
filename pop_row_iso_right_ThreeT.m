function [ROW1, ROW2, ROW3] = pop_row_iso_right_ThreeT(lambda1, lambda2, lambda3, Cp1, Cp2, Cp3, g12, g13, g23, dt ,dx, theta)

diff1 = lambda1/Cp1*dt/dx^2;
ch12 = dt*g12/Cp1;
ch13 = dt*g13/Cp1;
    
diff2 = lambda2/Cp2*dt/dx^2;
ch21 = dt*g12/Cp2;
ch23 = dt*g23/Cp2;
    
diff3 = lambda3/Cp3*dt/dx^2;
ch31 = dt*g13/Cp3;
ch32 = dt*g23/Cp3;

Ar1 = theta*diff1; %T1,i-1 to T1,i
Br1 = 0; %T2,i-1 to T1,i
Cr1 = 0; %T3,i-1 to T1,i
Dr1 = 1-2*diff1*theta - theta*(ch13+ch12); %T1i to T1i
Er1 = theta*ch12; %T1,i to T2,i
Fr1 = theta*ch13; %T1,i to T3,i
    
Ar2 = 0; %T1,i-1 to T2,i
Br2 = theta*diff2; %T2,i-1 to T2,i
Cr2 = 0; %T3,i-1 to T2,i
Dr2 = theta*ch21; %T2i to T1i
Er2 = 1-2*diff2*theta - theta*(ch21+ch23); %T2,i to T2,i
Fr2 = theta*ch23; %T2,i to T3,i

Ar3 = 0; %T1,i-1 to T3,i
Br3 = 0; %T2,i-1 to T3,i
Cr3 = theta*diff3; %T3,i-1 to T3,i
Dr3 = theta*ch31; %T3i to T1i
Er3 = theta*ch32; %T3,i to T2,i
Fr3 = 1-2*diff3*theta - theta*(ch31+ch32); %T3,i to T3,i

ROW1= [Ar1 Br1 Cr1 Dr1 Er1 Fr1 ];
ROW2 = [Ar2 Br2 Cr2 Dr2 Er2 Fr2 ];
ROW3 = [Ar3 Br3 Cr3 Dr3 Er3 Fr3 ];