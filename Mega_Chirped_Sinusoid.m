clear;
clc; 
close all;

%------------------------------ANF parameters --------------------------%
gamma = .006;
xi = .4;
%----------------------------------------------------------------------%

%-------------------------Few initializations -----------------------%
T = 0.001;                  % Sampling rate
a = 0;                      % t minimum
b = 1;                      % t maximum
t= [a:T:b];                 % total time interval
N = length(t);              % Number of samples
%----------------------------------------------------------------------%


%----------------------Preparing the signal------------------------------%
Y = zeros(1,N);

A = 1; 
phase = pi/3;
changingfreq = 60+ 30*t.^2;
PureSignal = A * sin (2*pi*changingfreq.*t+ phase);        % Persistant Excitation Signal

SNRdB = 20;
noise = - A + awgn(zeros(size(t))+A,SNRdB,'measured');
Y =PureSignal + noise;
figure, plot(t, Y, t, Y - PureSignal, 'LineWidth',.5)
noise = strcat ('Noise,SNR =',num2str(SNRdB), 'dB')
Ax = legend ('Quad chirped signal', noise)
leg = findobj(Ax,'type','text')
set(leg,'fontsize', 18)
xlabel ('Time (seconds)')
ylabel ('Amplitude (units)')
grid on 
axis tight
% print -depsc2 MegaChirpedSignal
%-----------------------------------------------------------------------%

%------------------------------ANF----------------------------------------%
initialFreq = 2*pi*60;
[X1, X2,X3]= NonUniform4thOrderANFFixedBlock(Y, t, initialFreq, gamma, xi);
%-------------------------------------------------------------------------%


%------------------------------Plotting the solutions ---------------------------%
figure, plot3(X1,X2,X3/(2*pi), 'LineWidth',.5)
% axis([-.25 .47 -100 230 30 80])
% title('3D state space plot from Agent based ODE solver Hsu ANF')
grid on
xlabel('x1')
ylabel('x2')
zlabel('theta (Hz)')
axis tight
% print -depsc2 MegaChirpedSignalANFEvolution

Indexes = find(X2 < 10^-3);     % locating point lying on the plane x_2 = 0 

figure, scatter (X1(Indexes ),X3 (Indexes)/ (2*pi) , '.')
% title('Poincare Section')
xlabel ('x_1')
ylabel ('\theta (Hz)')
grid on


% figure, plot(t, changingfreq,'m', t, x1(:,3)'/(2*pi), ':k' , t,X3/(2*pi), '-.r', t,MX3/(2*pi), '-.b')
figure, plot(t, changingfreq,':k','LineWidth',3)
title('(b)')
hold on
plot(t,X3/(2*pi),'-.g','LineWidth',2.5)
hold off
grid on
axis tight
Ax = legend('True value','Proposed discrete method','Location', 'northwest' )
leg = findobj(Ax,'type','text')
set(leg,'fontsize', 18)
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
% print -depsc2 MegaChirpedSignalANFTracking

%------------------------------------------------------------------------%