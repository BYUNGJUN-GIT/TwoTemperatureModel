function [ROW1, ROW2, ROW3] = pop_row_cond_left_ThreeT(lambda1, lambda2, lambda3, Cp1, Cp2, Cp3, g12, g13, g23,G, dt ,dx, theta)

diff1 = lambda1(end)/Cp1(end)*dt/dx^2;
ch12 = dt*g12(end)/Cp1(end);
ch13 = dt*g13(end)/Cp1(end);
  
diff2 = lambda2(end)/Cp2(end)*dt/dx^2;
ch21 = dt*g12(end)/Cp2(end);
ch23 = dt*g23(end)/Cp2(end);
    
diff3 = lambda3(end)/Cp3(end)*dt/dx^2;
ch31 = dt*g13(end)/Cp3(end);
ch32 = dt*g23(end)/Cp3(end);
    
cond11 = dt/dx*1/Cp1(end)*G(1);
cond12 = dt/dx*1/Cp1(end)*G(2);
cond13 = dt/dx*1/Cp1(end)*G(3);
cond21 = dt/dx*1/Cp2(end)*G(4);
cond22 = dt/dx*1/Cp2(end)*G(5);
cond23 = dt/dx*1/Cp2(end)*G(6);
cond31 = dt/dx*1/Cp3(end)*G(7);
cond32 = dt/dx*1/Cp3(end)*G(8);
cond33 = dt/dx*1/Cp3(end)*G(9);

Al1 = theta*cond11;
Bl1 = theta*cond12;
Cl1 = theta*cond13;
Dl1 = 1-theta*diff1- theta*ch12 - theta*ch13 - theta*(cond11+cond12+cond13);
El1 = theta*ch12;
Fl1 = theta*ch13;
Hl1 = theta*diff1;
Kl1 = 0;
Ll1 = 0;
      
Al2 = theta*cond21;
Bl2 = theta*cond22;
Cl2 = theta*cond23;
Dl2 = theta*ch21;
El2 = 1-theta*diff2 - theta*ch21-theta*ch23 - theta*(cond21+cond22+cond23);
Fl2 = theta*ch23;
Hl2 = 0;
Kl2 = theta*diff2;
Ll2 = 0;
    
Al3 = theta*cond31;
Bl3 = theta*cond32;
Cl3 = theta*cond33;
Dl3 = theta*ch31;
El3 =  theta*ch32;
Fl3 = 1-theta*diff3 - theta*ch31 - theta*ch32 - theta*(cond31+cond32+cond33);
Hl3 = 0;
Kl3 = 0;
Ll3 = theta*diff3;

ROW1= [Al1 Bl1 Cl1 Dl1 El1 Fl1 Hl1 Kl1 Ll1];
ROW2 = [Al2 Bl2 Cl2 Dl2 El2 Fl2 Hl2 Kl2 Ll2];
ROW3 = [Al3 Bl3 Cl3 Dl3 El3 Fl3 Hl3 Kl3 Ll3];