clear all;
close all;
clc;

%----------------------------Few Initializations-----------------------%
T = 0.001;
a = 0;
b = 1;
t= [a:T:b];
linewidth = 2.5;
%----------------------------------------------------------------------%

%---------------------Preparing the signal with 2 sinusoids---------------%
Y = zeros(size(t));
FirstFrequencySignal = zeros(size(t));
SecondFrequencySignal = zeros(size(t));
changingfreqs = zeros(2,length(t));
A1 = 1;(4/pi);
phi1 = pi/3;
A2 = 1/2;(4/pi);
phi2 = pi/7;

for i =1:length(t)
    if (t(i) <= .33)
%         changingfreqs(1,i)  = 72;
        changingfreqs(1,i)  = 60;
        changingfreqs(2,i)  = 120;
    end
    if (t(i) > .33 && t(i) <= .66)
%         changingfreqs(1,i)  = 57;
        changingfreqs(1,i)  = 60;
        changingfreqs(2,i)  = 120;
    end
    if (t(i) > .66)
%         changingfreqs(1,i)  = 80;
        changingfreqs(1,i)  = 60;
        changingfreqs(2,i)  = 120;
    end
    FirstFrequencySignal(i) =  A1 * sin (2*pi*changingfreqs(1,i)*t(i)+ phi1);
    SecondFrequencySignal(i) =  A2 * sin (2*pi*changingfreqs(2,i)*t(i)+ phi2);
    Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i); % Composite Excitation Signal
end

SNRdB = 20;
noise = - A1 + awgn(zeros(size(t))+A1,SNRdB,'measured');
Y =Y + noise;   % infecting the signal with additive noise
figure,plot(t, Y, t,FirstFrequencySignal, t,SecondFrequencySignal,t, noise, '-.', 'LineWidth', linewidth)
axis([0.5-.03 0.5+.03 -3  3])
ylabel ('Amplitude (units)')
xlabel('Time(seconds)')
noise = strcat('Noise,SNR=', num2str(SNRdB), ' dB')
legend('Composite','1^{st} Sinusoid','2^{nd} Sinusoid',noise)
grid on
%------------------------------------------------------------------------%

%---------------------------- ANF parameters -----------------------------%
gamma = 0.01;
xi  = .15;
%------------------------------------------------------------------------%

initialFreq1 =  2*pi*60-12;
initialFreq2 =  2*pi*120+24;

%----------------------------------Cascade Configuration------------------%
NewY = Y;
[X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
NewY = Y - 2*xi*X2./X3;                                   % Cascaded
[X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);

figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),'g',  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth)
% axis([a b 50 200])
% title ('Cascade Configuration')
title ('(a)')
Ax = legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
Ax.FontSize = 14;
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on
axis tight
%-------------------------------------------------------------------------%

%------------------------Parallel Configuration---------------------------%
NewY = Y;
[X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
NewY = Y;                                              % Identical Input
[X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);

figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),'g',  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth+.3)
% axis([a b 50 200])
% title ('Parallel Configuration')
Ax = legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
Ax.FontSize = 14;
title ('(b)')
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on
axis tight
%-------------------------------------------------------------------------%

% %----------------------Prefiltered Configuration -------------------------%
% r = 0.94;
% omega0 = initialFreq1*T;                  % Computing Notch at specific Freqs
% b1 = [1 -2*cos(omega0) 1];
% a1 = [1 -2*r*cos(omega0) r^2];
% NewY= filter(b1-a1,a1,Y);                           % Prefiltering
% [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
% 
% r = 0.93;
% omega0 = initialFreq2*T;                  % Computing Notch at specific Freqs
% b2 = [1 -2*cos(omega0) 1];
% a2 = [1 -2*r*cos(omega0) r^2];
% NewY= filter(b2-a2,a2,Y);                           % Prefiltering
% [X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);
% 
% 
% h = fvtool (b1-a1,a1,b2-a2,a2, 'Fs', 1/T);
% set(h,'DesignMask','off'); % Turn off design mask
% hchildren = get(h,'children');
% haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
% hline = get(haxes,'children');
% set(hline,'linewidth',linewidth)
% str1 = strcat('Bandpass, f_1 = ', num2str(initialFreq1/(2*pi)), 'Hz');
% str2 = strcat('Bandpass, f_2 = ', num2str(initialFreq2/(2*pi)), 'Hz');
% legend (str1,str2)
% 
% 
% figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth)
% axis([a b 50 200])
% % title ('Prefiltered Configuration')
% legend('True 1st Freq','True 2nd Freq','Estimate 1','Estimate 2' )
% ylabel('Frequency (Hz)')
% xlabel('Time (seconds)')
% grid on
%-------------------------------------------------------------------------%