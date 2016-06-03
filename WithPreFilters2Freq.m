clear all;
close all;
clc;
% 
% %----------------------------Few Initializations-----------------------%
% T = 0.001;
% a = 0;
% b = 1;
% t= [a:T:b];
% linewidth = 2.5;
% %----------------------------------------------------------------------%
% 
% %---------------------Preparing the signal with 2 sinusoids---------------%
% N = length(t);
% Y = zeros(1,N);
% FirstFrequencySignal = zeros(1,N);
% SecondFrequencySignal = zeros(1,N);
% changingfreqs = zeros(2,N);
% A1 = 1;(4/pi);
% phi1 = pi/3;
% A2 = 1;(4/pi);
% phi2 = pi/7;
% 
% for i =1:length(t)
%     if (t(i) <= .33)
%         changingfreqs(1,i)  = 72;
%         changingfreqs(2,i)  = 130;
%     end
%     if (t(i) > .33 && t(i) <= .66)
%         changingfreqs(1,i)  = 57;
%         changingfreqs(2,i)  = 150;
%     end
%     if (t(i) > .66)
%         changingfreqs(1,i)  = 80;
%         changingfreqs(2,i)  = 140;
%     end
%     %-------Make it static
%     changingfreqs(1,i) = 60;
%     changingfreqs(2,i)  = 120;
%     FirstFrequencySignal(i) =  A1 * sin (2*pi*changingfreqs(1,i)*t(i)+ phi1);
%     SecondFrequencySignal(i) =  A2 * sin (2*pi*changingfreqs(2,i)*t(i)+ phi2);
%     Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i) ;
% %     + A1 * sin (2*pi*170*t(i)+ phi1)+ A1 * sin (2*pi*260*t(i)+ phi1)+ A1 * sin (2*pi*359*t(i)+ phi1); % Composite Excitation Signal
% end
% 
% SNRdB = 0;
% noise = - 2*A1 + awgn(zeros(size(t))+2*A1,SNRdB,'measured');
% Y =Y + noise;   % infecting the signal with additive noise
% figure,plot(t, Y, t,FirstFrequencySignal, t,SecondFrequencySignal,t, noise, '-.', 'LineWidth', linewidth)
% axis([0.5-.03 0.5+.03 -3  3])
% ylabel ('Amplitude (units)')
% xlabel('Time(seconds)')
% noise = strcat('Noise,SNR=', num2str(SNRdB), ' dB')
% legend('Composite','1^{st} Sinusoid','2^{nd} Sinusoid',noise)
% grid on
% %------------------------------------------------------------------------%
% 
% %---------------------------- ANF parameters -----------------------------%
% gamma = 0.001;
% xi  = .15;
% %------------------------------------------------------------------------%
% 
% initialFreq1 =  2*pi*66;
% initialFreq2 =  2*pi*132;
% 
% %----------------------------------Cascade Configuration Special------------------%
% NewY = Y;
% 
% % [X1, X2,X3] = Special4thOderBP_NonUniform4thOrderANFFixedBlock (NewY, t, initialFreq1, gamma, xi, initialFreq1, xi+.5) 
% % 
% % % [X1, X2,X3] = SpecialBP_NonUniform4thOrderANFFixedBlock (NewY, t, initialFreq1, gamma, xi, initialFreq1, xi+.2,  0);
% % 
% % % [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
% % % NewY = Y - 2*xi*X2./X3;                                   % Cascaded
% % [X21, X22,X23]= Special4thOderBP_NonUniform4thOrderANFFixedBlock (NewY, t, initialFreq2, gamma, xi, initialFreq2, xi+.5) 
% % 
% % % [X21, X22,X23]= SpecialBP_NonUniform4thOrderANFFixedBlock (NewY, t, initialFreq2, gamma, xi, initialFreq2, xi+.2,  0); 
% % % NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);
% % 
% % figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth)
% % axis([a b 50 170])
% % %  title ('Band Pass Special combo')
% % Ax = legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
% % Ax.FontSize = 14;
% % ylabel('Frequency (Hz)')
% % xlabel('Time (seconds)')
% % grid on 
%   
% %----------------------------------Cascade Configuration------------------%
% NewY = Y;
% [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
% NewY = Y - 2*xi*X2./X3;                                    % Cascaded
% [X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);
% 
% figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth)
% axis([a b 50 170])
% %  title ('Cascade Configuration') 
% Ax = legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
% Ax.FontSize = 14;
% ylabel('Frequency (Hz)')
% xlabel('Time (seconds)')
% grid on
% 
% %-------------------------------------------------------------------------%
%  
% %------------------------Parallel Configuration---------------------------%
% % NewY = Y;
% % [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
% % NewY = Y;                                              % Identical Input
% % [X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);
% % 
% % figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth+.3)
% % axis([a b 50 200])
% %  title ('Parallel Configuration')
% % legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
% % ylabel('Frequency (Hz)')
% % xlabel('Time (seconds)') 
% % grid on
% %-------------------------------------------------------------------------%
% 
% %----------------------Prefiltered Configuration -------------------------%
% 
% r = 0.94;
% omega0 = initialFreq1*T;                  % Computing Notch at specific Freqs
% b1 = [1 -2*cos(omega0) 1];
% a1 = [1 -2*r*cos(omega0) r^2];
% NewY= filter(b1-a1,a1,Y);                           % Prefiltering
% 
% [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
% 
% r = 0.95;
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
%  title ('Prefiltered Configuration')
% legend('True 1st Freq','True 2nd Freq','Estimate 1','Estimate 2' )
% ylabel('Frequency (Hz)')
% xlabel('Time (seconds)')
% grid on
% %-------------------------------------------------------------------------%
% 





%%        Previous Version

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
A1 = 1; %(4/pi);
phi1 = pi/3;
phi1=pi/2;
A2 = 1;2; %(4/pi);
phi2 = pi/7;
% 
% for i =1:length(t)
%     if (t(i) <= .33)
%         changingfreqs(1,i)  = 72;
%         changingfreqs(2,i)  = 130;
%     end
%     if (t(i) > .33 && t(i) <= .66)
%         changingfreqs(1,i)  = 57;
%         changingfreqs(2,i)  = 150;
%     end
%     if (t(i) > .66)
%         changingfreqs(1,i)  = 80;
%         changingfreqs(2,i)  = 140;
%     end
%     FirstFrequencySignal(i) =  A1 * sin (2*pi*changingfreqs(1,i)*t(i)+ phi1);
%     SecondFrequencySignal(i) =  A2 * sin (2*pi*changingfreqs(2,i)*t(i)+ phi2);
%     Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i); % Composite Excitation Signal
% end


for i =1:length(t)
    if (t(i) <= .33)
        changingfreqs(1,i)  = 120;
        changingfreqs(2,i)  = 60;
    end
    if (t(i) > .33 && t(i) <= .66)
        changingfreqs(1,i)  = 120;
        changingfreqs(2,i)  = 60;
    end
    if (t(i) > .66)
        changingfreqs(1,i)  = 120;
        changingfreqs(2,i)  = 60;
    end
    FirstFrequencySignal(i) =  A1 * sin (2*pi*changingfreqs(1,i)*t(i)+ phi1);
    SecondFrequencySignal(i) =  A2 * sin (2*pi*changingfreqs(2,i)*t(i)+ phi2);
    Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i); % Composite Excitation Signal
end



SNRdB = 20;
noise = - 2 + awgn(zeros(size(t))+2,SNRdB,'measured');
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
% gamma = 0.05;
gamma=0.001
% xi  = .3;
xi=0.15
%------------------------------------------------------------------------%

initialFreq1 =  2*pi*56;
initialFreq2 =  2*pi*125;

%----------------------------------Cascade Configuration------------------%
NewY = Y;
[X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
NewY = Y - 2*xi*X2./X3;                                   % Cascaded
[X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);

figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth)
axis([a b 50 150])
%title ('Cascade Configuration')
legend('1st true frequency','2nd true frequency','1st frequency estimate','2nd frequency estimation' )
%legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on

%-------------------------------------------------------------------------%

% %------------------------Parallel Configuration---------------------------%
% NewY = Y;
% [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);
% NewY = Y;                                              % Identical Input
% [X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);
% 
% figure, plot(t, changingfreqs(1,:),t, changingfreqs(2,:),  t,X3/(2*pi), '-.k', t,X23/(2*pi), '-.r', 'LineWidth', linewidth+.3)
% %axis([a b 50 130])
% % title ('Parallel Configuration')
% legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
% ylabel('Frequency (Hz)')
% xlabel('Time (seconds)')
% grid on
% %-------------------------------------------------------------------------%

%----------------------Prefiltered Configuration -------------------------%
r = 0.94;
omega0 = initialFreq1*T;                  % Computing Notch at specific Freqs
b1 = [1 -2*cos(omega0) 1];
a1 = [1 -2*r*cos(omega0) r^2];
NewY= filter(b1-a1,a1,Y);                           % Prefiltering
[Z1, Z2,Z3]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq1, gamma, xi);

r = 0.93;
omega0 = initialFreq2*T;                  % Computing Notch at specific Freqs
b2 = [1 -2*cos(omega0) 1];
a2 = [1 -2*r*cos(omega0) r^2];
NewY= filter(b2-a2,a2,Y);                           % Prefiltering
[Z21, Z22,Z23]= NonUniform4thOrderANFFixedBlock(NewY, t, initialFreq2, gamma, xi);


h = fvtool (b1-a1,a1,b2-a2,a2, 'Fs', 1/T);
set(h,'DesignMask','off'); % Turn off design mask
hchildren = get(h,'children');
haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
hline = get(haxes,'children');
set(hline,'linewidth',linewidth)
str1 = strcat('Bandpass, f_1 = ', num2str(initialFreq1/(2*pi)), 'Hz');
str2 = strcat('Bandpass, f_2 = ', num2str(initialFreq2/(2*pi)), 'Hz');
Ax = legend (str1,str2);
Ax.FontSize = 14;
 title ('(b)')

figure, plot(t, changingfreqs(1,:), 'g', t, changingfreqs(2,:), 'g',  t,X3/(2*pi), 'b', t,X23/(2*pi), 'b', t, Z3/(2*pi), '-.r', t,Z23/(2*pi), '-.r', 'LineWidth', linewidth) %, LineWidth'), %linewidth) %, linewidth)
title('(a)')
axis([a b 50 130])
% title ('Prefiltered Configuration')
Ax = legend('True 1^{st} Freq',' True 2^{nd} Freq','Cascade 1^{st} Freq Estimate','Cascade 2^{nd} Freq Estimate', 'Prefiltered 1^{st} Freq Estimate','Prefiltered 2^{nd} Freq Estimate' )
Ax.FontSize = 14;
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on
%-------------------------------------------------------------------------%
