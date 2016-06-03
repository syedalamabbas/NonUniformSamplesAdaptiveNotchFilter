clear;
clc;

% Given  ANF parameters and initial conditions
% x_1_init0 = 0;
% x_2_init0 = 0;
% theta_init0 = 2*pi*60;
% gamma = .001;
% xi  = .05;

% The system of differential Equations
% Dx_1 = x_2;
% Dx_2 = -2*xi*theta*x_2 - theta^2*x_1 + theta^2*y;
% Dtheta = gamma*(2*xi*theta*x_2 - theta^2*y)*x_1;

% Using truncated solution of upto 5 terms for each variable
syms t x_1_0 x_1_1 x_1_2 x_1_3 x_1_4 x_2_0 x_2_1 x_2_2 x_2_3 x_2_4  theta_0 theta_1 theta_2 theta_3 theta_4

syms A B C D E      % Functional Adomain Polynomials

A = x_2_0 + x_2_1*t + x_2_2* t^2 + x_2_3 *t^3 + x_2_4* t^4

B = collect( expand ( (theta_0 + theta_1*t + theta_2* t^2+ theta_3*t^3+ theta_4*t^4)*(x_2_0 + x_2_1*t + x_2_2* t^2 + x_2_3 *t^3 + x_2_4* t^4) ),t) 

C = collect( expand ( (theta_0 + theta_1*t + theta_2* t^2+ theta_3*t^3+ theta_4*t^4)^2*(x_1_0 + x_1_1*t + x_1_2* t^2 + x_1_3 *t^3 + x_1_4* t^4) ),t)

D = collect( expand ( (theta_0 + theta_1*t + theta_2* t^2+ theta_3*t^3+ theta_4*t^4)^2),t)

E = collect( expand ( (theta_0 + theta_1*t + theta_2* t^2+ theta_3*t^3+ theta_4*t^4)*(x_1_0 + x_1_1*t + x_1_2* t^2 + x_1_3 *t^3 + x_1_4* t^4) * (x_2_0 + x_2_1*t + x_2_2* t^2 + x_2_3 *t^3 + x_2_4* t^4) ),t)

% Therefore the Adomain Polynomial Coefficients are
A_0 = x_2_0;
A_1 = x_2_1;
A_2 = x_2_2; 
A_3 = x_2_3;
A_4 = x_2_4;

B_0 =  theta_0*x_2_0;
B_1 = (theta_0*x_2_1 + theta_1*x_2_0);
B_2 = (theta_0*x_2_2 + theta_1*x_2_1 + theta_2*x_2_0);
B_3 = (theta_0*x_2_3 + theta_1*x_2_2 + theta_2*x_2_1 + theta_3*x_2_0);
B_4 = (theta_0*x_2_4 + theta_1*x_2_3 + theta_2*x_2_2 + theta_3*x_2_1 + theta_4*x_2_0);

C_0 =  theta_0^2*x_1_0;
C_1 = (x_1_1*theta_0^2 + 2*theta_1*x_1_0*theta_0);
C_2 = (x_1_2*theta_0^2 + 2*x_1_1*theta_0*theta_1 + 2*theta_2*x_1_0*theta_0 + x_1_0*theta_1^2);
C_3 = (theta_1^2*x_1_1 + theta_0^2*x_1_3 + 2*theta_0*theta_1*x_1_2 + 2*theta_0*theta_2*x_1_1 + 2*theta_0*theta_3*x_1_0 + 2*theta_1*theta_2*x_1_0);
C_4 = (theta_2^2*x_1_0 + theta_1^2*x_1_2 + theta_0^2*x_1_4 + 2*theta_0*theta_1*x_1_3 + 2*theta_0*theta_2*x_1_2 + 2*theta_0*theta_3*x_1_1 + 2*theta_0*theta_4*x_1_0 + 2*theta_1*theta_2*x_1_1 + 2*theta_1*theta_3*x_1_0);

D_0 = theta_0^2;
D_1 = (2*theta_0*theta_1);
D_2 = (theta_1^2 + 2*theta_0*theta_2);
D_3 = (2*theta_0*theta_3 + 2*theta_1*theta_2);
D_4 = (theta_2^2 + 2*theta_0*theta_4 + 2*theta_1*theta_3);

E_0 = theta_0*x_1_0*x_2_0;
E_1 = (theta_0*x_1_0*x_2_1 + theta_0*x_1_1*x_2_0 + theta_1*x_1_0*x_2_0);
E_2 = (theta_0*x_1_0*x_2_2 + theta_0*x_1_1*x_2_1 + theta_0*x_1_2*x_2_0 + theta_1*x_1_0*x_2_1 + theta_1*x_1_1*x_2_0 + theta_2*x_1_0*x_2_0);
E_3 = (theta_0*x_1_0*x_2_3 + theta_0*x_1_1*x_2_2 + theta_0*x_1_2*x_2_1 + theta_0*x_1_3*x_2_0 + theta_1*x_1_0*x_2_2 + theta_1*x_1_1*x_2_1 + theta_1*x_1_2*x_2_0 + theta_2*x_1_0*x_2_1 + theta_2*x_1_1*x_2_0 + theta_3*x_1_0*x_2_0);
E_4 = (theta_0*x_1_0*x_2_4 + theta_0*x_1_1*x_2_3 + theta_0*x_1_2*x_2_2 + theta_0*x_1_3*x_2_1 + theta_0*x_1_4*x_2_0 + theta_1*x_1_0*x_2_3 + theta_1*x_1_1*x_2_2 + theta_1*x_1_2*x_2_1 + theta_1*x_1_3*x_2_0 + theta_2*x_1_0*x_2_2 + theta_2*x_1_1*x_2_1 + theta_2*x_1_2*x_2_0 + theta_3*x_1_0*x_2_1 + theta_3*x_1_1*x_2_0 + theta_4*x_1_0*x_2_0);












