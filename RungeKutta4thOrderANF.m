function SolutionX = RungeKutta4thOrderANF (t,X,freq,SNRdB ,xi, gamma)
% This function computes the classical 4th order numerical approximation
% method popularly known as Runge-Kutta Method
N = length(t);
SolutionX = zeros(N,3);
SolutionX(1,:) = X;   % Intial conditions given
for i = 2:N
    h = double(t(i)- t(i-1));
    k1 = h*f(t(i),X, freq, SNRdB, xi, gamma);
    k2 = h*f(t(i)+h/2, X+k1/2,freq, SNRdB, xi, gamma);
    k3 = h*f(t(i)+h/2, X+k2/2,freq,SNRdB, xi, gamma);
    k4 = h*f(t(i)+h, X+k3,freq,SNRdB, xi, gamma);
    X = X + (k1+2*k2+2*k3+k4)/6;
    SolutionX(i,:) = X;
end
end

% Local function, describing the dynamic system of Hsu Adaptive Notch Filter

function dX = f (t,X,freq,SNRdB, xi, gamma)
% Observe here the Hard coded Excitation to compute the inbetween samples
A = 1;
% f= 60;        % No hard coded frequency , it is passed as a parameter to the function , so is the SNR
phi = pi/2;
excitation = A*sin(2*pi*freq*t + phi)  - A + awgn(zeros(size(t)) + A ,SNRdB,'measured');

dX = zeros(size(X));
dX(1) = X(2);
dX(2) = -2*xi*X(3)*X(2)-X(3)^2*X(1)+ X(3)^2*excitation;
dX(3) = -gamma*(X(3)^2*excitation - 2*xi*X(3)*X(2))*X(1);
end
