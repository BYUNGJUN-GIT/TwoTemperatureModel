function [ROW1, ROW2, ROW3] =  pop_row_left_ThreeT(lambda1, lambda2, lambda3, Cp1, Cp2, Cp3, g12, g13, g23, dt ,dx, theta)

diff1 = lambda1/Cp1*dt/dx^2;
ch12 = dt*g12/Cp1;
ch13 = dt*g13/Cp1;
    
diff2 = lambda2/Cp2*dt/dx^2;
ch21 = dt*g12/Cp2;
ch23 = dt*g23/Cp2;
    
diff3 = lambda3/Cp3*dt/dx^2;
ch31 = dt*g13/Cp3;
ch32 = dt*g23/Cp3;

Dl1 = 1-theta*diff1 - theta*ch12 - theta*ch13; %T1i to T1i
El1 = theta*ch12; %T1,i to T2,i
Fl1 = theta*ch13; %T1,i to T3,i
Hl1 = diff1*(theta); %T1,i to T1,i+1
Kl1 = 0; %T1,i to T2,i+1
Ll1 = 0; %T1,i to T3,i+1
   
ROW1 = [Dl1 El1 Fl1 Hl1 Kl1 Ll1];

Dl2 = ch21*theta; %T2i to T1i
El2 = 1-theta*diff2 - theta*ch21 - theta*ch23; %T2,i to T2,i
Fl2 = theta*ch23; %T2,i to T3,i
Hl2 = 0; %T2,i to T1,i+1
Kl2 = theta*diff2; %T2,i to T2,i+1
Ll2 = 0; %T2,i to T3,i+1

ROW2 = [Dl2 El2 Fl2 Hl2 Kl2 Ll2];  

Dl3 = theta*ch31; %T3i to T1i
El3 = theta*ch32; %T3,i to T2,i
Fl3 = 1-theta*diff3 - theta*ch31 - theta*ch32; %T3,i to T3,i
Hl3 = 0; %T3,i to T1,i+1
Kl3 = 0; %T3,i to T2,i+1
Ll3 = theta*diff3; %T3,i to T3,i+1

ROW3 = [Dl3 El3 Fl3 Hl3 Kl3 Ll3];