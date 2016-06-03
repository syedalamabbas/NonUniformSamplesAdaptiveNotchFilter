clear;
clc;
close all;


%---------------------------ANF parameters--------------------------------%
gamma = .0011;                  % Adaptation speed
xi = .1;                      % Adaptation bandwidth
%-------------------------------------------------------------------------%

%--------------------------- Few initializations--------------------------%
T = 0.001;                          % sampling rate
a = 0;                              % t -minimum    
b = 1;                              % t- maximum
t= [a+T:T:b];                       % total time interval
N = length(t);
%-------------------------------------------------------------------------%


%----------------------------Signal Parameters --------------------------%
Y = zeros(1,N);

A = 1;
phase = pi/3;

startFreq = 60;
ramp_rate = 30;

changingfreq = startFreq  + ramp_rate *t;                  % Linear chirp

PureSignal = A * sin (2*pi*changingfreq.*t+ phase);        % Persistant Excitation Signal

SNRdB = 20;
noise = - A + awgn(zeros(size(t))+A,SNRdB,'measured');
Y = PureSignal + noise;
%-----------------------------------------------------------------------%


%--------------------------Plotting the ssignal-----------------------%
figure, plot(t, Y, t, Y - PureSignal, 'LineWidth',.5)
noise = strcat ('Noise,SNR =',num2str(SNRdB), 'dB')
Ax = legend ('Linear chirped signal', noise)
leg = findobj(Ax,'type','text')
set(leg,'fontsize', 18)
xlabel ('Time (seconds)')
ylabel ('Amplitude (units)')
grid on 
axis tight
% print -depsc2 ChirpedSignal
%---------------------------------------------------------------------%

%----------------------------ANF solution---------------------------------%
initialFreq = 2*pi*60;
[X1, X2,X3]= NonUniform4thOrderANFFixedBlock(Y, t, initialFreq, gamma, xi);
%------------------------------------------------------------------------%


%------------------------------Plotting the solutions--------------------%
figure, plot3(X1,X2,X3/(2*pi), 'LineWidth',.5)
% axis([-.25 .47 -100 230 30 80])
% title('3D state space plot from Agent based ODE solver Hsu ANF')
grid on
xlabel('x_1')
ylabel('x_2')
zlabel('\theta (Hz)')
axis tight
% print -depsc2 ChirpedSignalANFEvolution


Indexes = find(X2 < 10^-3);     % locating point lying on the plane x_2 = 0 

figure, scatter (X1(Indexes ),X3 (Indexes)/ (2*pi) , '.')
% title('Poincare Section')
xlabel ('x_1')
ylabel ('\theta (Hz)')
grid on


figure, plot(t, changingfreq,':k','LineWidth',3)
title('(a)')
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
% print -depsc2 ChirpedSignalANFTracking
%-------------------------------------------------------------------------%