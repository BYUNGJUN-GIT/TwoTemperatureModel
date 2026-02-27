function [ROW1, ROW2, ROW3] = pop_row_cond_right_ThreeT(lambda1, lambda2, lambda3, Cp1, Cp2, Cp3, g12, g13, g23,G, dt ,dx, theta)

cond11 = dt/dx*1/Cp1*G(1);
cond12 = dt/dx*1/Cp1*G(2);
cond13 = dt/dx*1/Cp1*G(3);
cond21 = dt/dx*1/Cp2*G(4);
cond22 = dt/dx*1/Cp2*G(5);
cond23 = dt/dx*1/Cp2*G(6);
cond31 = dt/dx*1/Cp3*G(7);
cond32 = dt/dx*1/Cp3*G(8);
cond33 = dt/dx*1/Cp3*G(9);

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
Dr1 = 1-theta*diff1 - theta*ch12 - theta*ch13 - theta*(cond11+cond12+cond13); %T1i to T1i
Er1 = theta*ch12; %T1,i to T2,i
Fr1 = theta*ch13; %T1,i to T3,i
Hr1 = theta*cond11; %T1,i to T1,i+1
Kr1 = theta*cond12; %T1,i to T2,i+1
Lr1 = theta*cond13; %T1,i to T3,i+1
    
Ar2 = 0; %T1,i-1 to T2,i
Br2 = theta*diff2; %T2,i-1 to T2,i
Cr2 = 0; %T3,i-1 to T2,i
Dr2 = theta*ch21; %T2i to T1i
Er2 = 1-theta*diff2 - theta*ch21 - theta*ch23 - theta*(cond21+cond22+cond23); %T2,i to T2,i
Fr2 = theta*ch23; %T2,i to T3,i
Hr2 = theta*cond21; %T2,i to T1,i+1
Kr2 = theta*cond22; %T2,i to T2,i+1
Lr2 = theta*cond23; %T2,i to T3,i+1
 
Ar3 = 0; %T1,i-1 to T3,i
Br3 = 0; %T2,i-1 to T3,i
Cr3 = theta*diff3; %T3,i-1 to T3,i
Dr3 = theta*ch31; %T3i to T1i
Er3 = theta*ch32; %T3,i to T2,i
Fr3 = 1-theta*diff3 - theta*ch31 - theta*ch32 - theta*(cond31+cond32+cond33); %T3,i to T3,i
Hr3 = theta*cond31; %T3,i to T1,i+1
Kr3 = theta*cond32; %T3,i to T2,i+1
Lr3 = theta*cond33; %T3,i to T3,i+1

ROW1= [Ar1 Br1 Cr1 Dr1 Er1 Fr1 Hr1 Kr1 Lr1];
ROW2 = [Ar2 Br2 Cr2 Dr2 Er2 Fr2 Hr2 Kr2 Lr2];
ROW3 = [Ar3 Br3 Cr3 Dr3 Er3 Fr3 Hr3 Kr3 Lr3];