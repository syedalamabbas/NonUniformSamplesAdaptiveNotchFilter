function [ xnew ] = FunctionXiaAdaptive( t, x )
%FUNCTIONANF Computes the value of the state-space system defined by diff equations
% x = state of system
% t = time value 
% intial values and excitation input are external supplied

b = .3;
k_1 = .5;
k_2 = .5;
gamma = .6;

xnew = zeros (4,1);
freq  = 60;
A = 1;
y = A*sin(2*pi*freq*t + pi/3);     % Pure sinusoid y(t)

% SNRdB = 5;
% y = y - A + awgn(zeros(size(t))+ A ,SNRdB,'measured')  ;
 
%--------------------------Xia 2002---------------------------%
xnew(1) = -b*x(1)-y;
xnew(2) = b*x(1)*x(4)+k_1*(y-x(3));
xnew(3) = x(2)+x(1)*x(4)+k_2*(y-x(3));
xnew(4) = gamma*x(1)*(y-x(3));
%-------------------------------------------------------------%
end