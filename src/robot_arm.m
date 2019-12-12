%% Linearization
clear;
% Parameters
L = 0.0036;
R = 0.5;
d11 = 1;
d12 = -0.13;
d22 = 0.24;
b1 = 3.6;
b2 = 0.7;
% Variables
syms i1 i2 di1 di2 q1 q2 dq1 dq2 ddq1 ddq2 u1 u2
% Calculate df/dx and df/du
x = [i1,q1,dq1,i2,q2,dq2]';
u = [u1,u2]';
s21 = sin(q2-q1);
c21 = cos(q2-q1);
ddq1 = (d12*s21*(b2*dq2 + d12*c21*dq1^2 - i2) - d22*(b1*dq1 + d12*c21*dq2^2 - i1))/(d11*d22-d12^2*s21^2);
ddq2 = (d12*s21*(b1*dq1 + d12*c21*dq2^2 - i1) - d22*(b1*dq2 + d12*c21*dq1^2 - i2))/(d11*d22-d12^2*s21^2);
di1 = (u1 - R*i1 - dq1)/L;
di2 = (u2 - R*i2 - dq2)/L;
f = [di1,dq1,ddq1,di2,dq2,ddq2]';
dfdx = jacobian(f,x');
dfdu = jacobian(f,u');
% For xs = 0 and us = 0
i1 = 0;
i2 = 0;
q1 = 0;
dq1 = 0;
q2 = 0;
dq2 = 0;
u1 = 0;
u2 = 0;
% Calculatec coefficient A and B
A = eval(dfdx);
B = eval(dfdu);
f = A*x + B*u;

%% LQR

%% Simulation
