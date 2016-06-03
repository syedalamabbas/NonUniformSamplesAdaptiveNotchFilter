function [ xnew ] = FunctionANF( t, x )
%  FUNCTIONANF Computes the value of the state-space system defined by diff
%  equations of Hsu ANF
%  x = state of system
%  t = time value 
%  intial values and excitation input are external supplied

% gamma = .9;
% xi = 1.4;

gamma = .01;
xi = .15;
xnew = zeros (3,1);

changingfreq  = 72;
if (t <= .33)
    changingfreq  = 72;
end
if (t > .33 && t <= .66)
    changingfreq  = 60 + 0*10 +0*(t-2);
end
if (t > .66)
     changingfreq  = 80 - 0* 15+ 0*(t-4);
end

A = 1; 
excitation = A*sin(2*pi*changingfreq*t + pi/3);     % Pure sinusoid y(t)

SNRdB = 20;
excitation = excitation - A + awgn(zeros(size(t))+A,SNRdB,'measured');

%----------------------Direct Hsu type (1999)------------------------%
xnew(1) = x(2);
xnew(2) = -x(3)^2*x(1)-2*xi*x(3)*x(2) + x(3)^2*excitation;
xnew(3) = -gamma*x(1)*x(3)^2*excitation + gamma*2*xi*x(1)*x(2)*x(3);
%-------------------------------------------------------------%
end