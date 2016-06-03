clear 
clc
% syms x_1(t) x_2(t) theta(t) xi gamma y(t)
% 
% 
% A = x_2;
% B = -2*xi*theta*x_2 - theta^2*x_1 + theta^2*y;
% C = -gamma*(theta^2*y - 2*xi*theta*x_2)*x_1;
% 
% 
% %% First order derivatives
% fprintf ( '\n -------------First order derivatives--------- \n ')
% disp (['Dx_1=',char(A)])
% disp (['Dx_2=',char(B)])
% disp (['D\theta=',char(C)])
% 
% order = 1;
% %% Second order derivatives
% fprintf ( '\n -----------Second order derivatives---------- \n ')
% disp (['D2x_1=',char(diff(A,order))])
% disp (['D2x_2=',char(diff(B,order))])
% disp (['D2\theta=',char(diff(C,order))])
% 
% order = order + 1;
% %% Third order derivatives
% fprintf ( '\n --------------Third order derivatives----------- \n ')
% disp (['D3x_1=',char(diff(A,order))])
% disp (['D3x_2=',char(diff(B,order))])
% disp (['D3\theta=',char(diff(C,order))])
% 
% order = order + 1;
% %% Fourth order derivatives
% fprintf ( '\n ---------------Fourth order derivatives------------- \n ')
% disp (['D4x_1=',char(diff(A,order))])
% disp (['D4x_2=',char(diff(B,order))])
% disp (['D4\theta=',char(diff(C,order))])
% 
% % order = order + 1;
% % %% Fifth order derivatives
% % disp (['D5x_1=',char(diff(A,4))])
% % disp (['D5x_2=',char(diff(B,4))])
% % disp (['D5\theta=',char(diff(C,4))])


syms x_1(t) x_2(t) x_3(t) x_4(t) theta(t) theta_c xi_c xi gamma y(t) Dx_4(t)

A = x_2; 
B = -2*xi_c*theta_c*x_2 - theta_c^2*x_1 + 2*xi_c*theta_c*diff(y(t), t);
C = x_4;
D = -2*xi*theta*x_4 - theta^2*x_3 + theta^2*x_1;
E = -gamma*(Dx_4 + theta^2*x_3)*x_3;
 
%% First order derivatives
fprintf ( '\n -------------MATLAB::First order derivatives--------- \n ')
disp (['Dx_1=',char(A)])
disp (['Dx_2=',char(B)])
disp (['Dx_3=',char(C)])
disp (['Dx_4=',char(D)])
disp (['D\theta=',char(E)])

order = 1;
%% Second order derivatives
fprintf ( '\n -----------MATLAB::Second order derivatives---------- \n ')
disp (['D^2x_1=',char(diff(A,order))])
disp (['D^2x_2=',char(diff(B,order))])
disp (['D^2x_3=',char(diff(C,order))])
disp (['D^2x_4=',char(diff(D,order))])
disp (['D^2\theta=',char(diff(E,order))])

order = order + 1;
%% Third order derivatives
fprintf ( '\n --------------MATLAB::Third order derivatives----------- \n ')
disp (['D^3x_1=',char(diff(A,order))])
disp (['D^3x_2=',char(diff(B,order))])
disp (['D^3x_3=',char(diff(C,order))])
disp (['D^3x_4=',char(diff(D,order))])
disp (['D^3\theta=',char(diff(E,order))])

order = order + 1;
%% Fourth order derivatives
fprintf ( '\n ---------------MATLAB::Fourth order derivatives------------- \n ')
disp (['D^4x_1=',char(diff(A,order))])
disp (['D^4x_2=',char(diff(B,order))])
disp (['D^4x_3=',char(diff(C,order))])
disp (['D^4x_4=',char(diff(D,order))])
disp (['D^4\theta=',char(diff(E,order))])

