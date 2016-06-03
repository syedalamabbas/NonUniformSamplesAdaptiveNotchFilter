clear all;
close all;
clc;

%----------------------------Few Initializations--------------------------%
T = 0.001;
a = 0;
b = 3;
t= [a:T:b];
linewidth = 2.5;
%-------------------------------------------------------------------------%


%-----------------------------Tuned ANF parameters------------------------%
gamma = 0.00001;
xi  = .01;
%-------------------------------------------------------------------------%


%--------------------------Preparing the signal with 2 Sines--------------%
Y = zeros(size(t));
FirstFrequencySignal = zeros(size(t));
SecondFrequencySignal = zeros(size(t));
freqs = [61 65];
A1 = 1; 
phi1 = pi/3;
A2 = 1/2;
phi2 = pi/7;

for i =1:length(t)
    FirstFrequencySignal(i) =  A1 * sin (2*pi*freqs(1,1)*t(i)+ phi1);
    SecondFrequencySignal(i) =  A2 * sin (2*pi*freqs(1,2)*t(i)+ phi2);
    Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i); % Composite Excitation Signal
end

SNRdB = 0;
PureSignal = Y;
Y =awgn(Y,SNRdB,'measured');   % infecting the signal with noise
figure,plot(t, Y, t,FirstFrequencySignal, 'g', t,SecondFrequencySignal, 'r',t , Y - PureSignal, '-.', 'LineWidth', linewidth)
axis([0.5-.03 0.5+.03 -3  3])
ylabel ('Amplitude (units)')
xlabel('Time(seconds)')
noise = strcat('Noise,SNR=', num2str(SNRdB), ' dB')
legend('Composite','1^{st} Signal, 61 Hz','2^{nd} Signal, 65 Hz',noise)
grid on

initialFreq1 =  2*pi*63;


%----------------------Prefiltered Configuration -------------------------%
r = 0.94;
omega0 = initialFreq1*T;                  % Computing Notch at specific Freqs
b1 = [1 -2*cos(omega0) 1];
a1 = [1 -2*r*cos(omega0) r^2];
NewY= filter(b1-a1,a1,Y);                           % Prefiltering
[X1, X2,X3]= FourthOrderANFFixedBlock(NewY, T, initialFreq1, gamma, xi);


% h = fvtool (b1-a1,a1, 'Fs', 1/T);
% set(h,'DesignMask','off'); % Turn off design mask
% hchildren = get(h,'children');
% haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
% hline = get(haxes,'children');
% set(hline,'linewidth',linewidth)
% str1 = strcat('Bandpass, f_1 = ', num2str(initialFreq1/(2*pi)), 'Hz');
% legend (str1)
% print -depsc2 TwoBandPassFilters


figure, plot(t, freqs(1,1)*ones(size(t)), 'g',t, freqs(1,2)*ones(size(t)),'r' , t,X3/(2*pi), '-.k', 'LineWidth', linewidth)
axis([a b 60 66])
% title ('Prefiltered Configuration')
legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate' )
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on
% Case 1 : A1 = 1; A2 = 1/3;
% Case 2 : A1 = 1/3; A2 = 1;
% Case 3 : A1 = 1; A2 = 1;
%-------------------------------------------------------------------------%