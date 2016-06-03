

clear all;
close all;
clc


N=1200;

fact = 8;
Fs=8000/fact;                % Comes out to be 1kHz
tm=[0:1/Fs:(N-1)/Fs]';
t=[0:400:N]+1;               % Time Partition 
snr=10^(-18/20);
F=[1000 1075 975]/fact;      % Frequency Partition 
x=[];
true=[];
for i=1:3
    T=2*pi*F(i).*tm(t(i):t(i+1)-1);
%     sig=sin(T)+0.5*cos(T*2)+0.25*cos(T*3)+randn(N/3,1).*snr;
    sig=sin(T)+randn(N/3,1).*snr;
    x=[x;sig];
    true=[true;ones(N/3,1)*F(i)];
end


%% Plotting the signal
 modifiedt = (0:1:length(x)-1)*1/Fs;
figure, plot( modifiedt , x)
legend('Noisy Signal')
title('The signal with changing frequencies')
xlabel('Time(sec)')
ylabel('Amplitude')
grid on

%% My solution initializations
gamma = 0.005;
xi  = .14; 
initialFreq = 2*pi*127;



%% Results 
[theta,theta_curve,b,a]=harmonic_est(x,1,Fs);                                          % This is the 2009 paper in IEEE Signal Processing Letters
[H1, H2,H3]= NonUniform4thOrderANFFixedBlock(x, modifiedt, initialFreq, gamma, xi);    % This is my solution, it can be a non-uniformly sampled system

our_Estimate = H3'/(2*pi);
% subplot(211)
figure, plot(tm,theta_curve)
hold on
plot(tm, our_Estimate, 'g')
plot(tm,true,'r--','LineWidth',3)
grid on
xlabel('Time (sec)')
ylabel('Fundamental Frequency Estimate (Hz)')
legend('Tracking IIR (Jiang)','Proposed Method','True')
% subplot(212)
% [H,F]=freqz(b,a,N,Fs);
% plot(F,log10(abs(H)))
% title('Final Comb Filter')
% xlabel('Frequency')
% ylabel('Magnitude')