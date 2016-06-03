clear all;
 close all;
clc;
%----------------------------Few Initializations--------------------------%
T = 0.001;
a = 0;
b = 3;
t= [a:T:b];
linewidth = 2.5;
fontSize = 14;
%-------------------------------------------------------------------------%

N = length(t);
modifiedt = zeros (1, N);

min = -5;                   % Min value
max = 5;                    % Max value  
r = (max-min).*rand(1,N) + min;  % Randomly distributed Values

for i = 2:N                               % Jittered model with additive non-uniform sampling
   modifiedt(i) = modifiedt(i-1) + T+ r(i)*T*10^-1;
end

%-----------------------------Tuned ANF parameters------------------------%
% gamma = 0.001;
gamma=0.00001;
% xi  = .15;
xi=0.015;
%-------------------------------------------------------------------------%
A1 = 1; 
phi1 = pi/2;
A2 = 1;
phi2 = pi/7;

Y = zeros(size(modifiedt));
changingfreqs = zeros (2, length(modifiedt));
for i =1:length(modifiedt)
    if (modifiedt(i) <= .33)
        changingfreqs(1,i)  = 61;
        changingfreqs(2,i)  = 65;
    end
    if (modifiedt(i) > .33 && modifiedt(i) <= .66)
        changingfreqs(1,i)  = 61;
        changingfreqs(2,i)  = 65;
    end
    if (modifiedt(i) > .66)
        changingfreqs(1,i)  = 61;
        changingfreqs(2,i)  = 65;
    end
    changingfreqs(1,i)  = 61;
    changingfreqs(2,i)  = 65;
    FirstFrequencySignal =  A1 * sin (2*pi*changingfreqs(1,i)*modifiedt(i)+ phi1);
    SecondFrequencySignal=  A2 * sin (2*pi*changingfreqs(2,i)*modifiedt(i)+ phi2);
    Y(i) = FirstFrequencySignal + SecondFrequencySignal + A1 * sin (2*pi*110*modifiedt(i)+ phi2)+ + A1 * sin (2*pi*510*modifiedt(i)+ phi2) + A1 * sin (2*pi*310*modifiedt(i)+ phi2); % Composite Excitation Signal
end
%--------------------------Preparing the signal with 2 Sines--------------%
% Y = zeros(size(t));
% FirstFrequencySignal = zeros(size(t));
% SecondFrequencySignal = zeros(size(t));
freqs = [61 65];

% for i =1:length(t)
%     FirstFrequencySignal(i) =  A1 * sin (2*pi*freqs(1,1)*t(i)+ phi1);
%     SecondFrequencySignal(i) =  A2 * sin (2*pi*freqs(1,2)*t(i)+ phi2);
%     Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i); % Composite Excitation Signal
% end
SNRdB = 0;
PureSignal = Y;
noise = - 2*A1 + awgn(zeros(size(t))+2*A1,SNRdB,'measured');

Y = Y + noise; awgn(Y,SNRdB,'measured');   % infecting the signal with noise
% figure,plot(t, Y, t,FirstFrequencySignal, 'g', t,SecondFrequencySignal, 'r',t , Y - PureSignal, '-.', 'LineWidth', linewidth)
% axis([0.5-.03 0.5+.03 -3  3])
% ylabel ('Amplitude (units)')
% xlabel('Time(seconds)')
% noise = strcat('Noise,SNR=', num2str(SNRdB), ' dB')
% legend('Composite','1^{st} Signal, 61 Hz','2^{nd} Signal, 65 Hz',noise)
% grid on
initialFreq1 =  2*pi*63;
%initialFreq1 =  2*pi*58;
initialFreq2 =  2*pi*67;
%----------------------------------Cascade Configuration------------------%
NewY = Y;
% [ X1,X2,X3] = Special4thOderBP_NonUniform4thOrderANFFixedBlock (NewY, modifiedt, initialFreq1, gamma, xi, initialFreq1, 0.25);
% [ X1,X2,X3] = SpecialBP_NonUniform4thOrderANFFixedBlock (NewY, modifiedt, initialFreq1, gamma, xi, initialFreq1, xi,  0);
[X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq2, gamma, xi);
NewY = Y - 2*xi*X2./X3(N);                                   % Cascaded
% [X21, X22,X23]= Special4thOderBP_NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq2, gamma, xi, initialFreq1, 0.25);

[X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq1, gamma, xi);


 figure, plot(modifiedt, changingfreqs(1,:), 'g', modifiedt, changingfreqs(2,:), 'b',  modifiedt,X3/(2*pi), '-.k', modifiedt,X23/(2*pi), '-.r', 'LineWidth', linewidth)
 axis([a b 60 70])
% title ('0 dB')
Ax = legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
Ax.FontSize = fontSize;
%legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on


% NewY = Y;
% % [ X1,X2,X3] = Special4thOderBP_NonUniform4thOrderANFFixedBlock (NewY, modifiedt, initialFreq1, gamma, xi, initialFreq1, 0.25);
% % [ X1,X2,X3] = SpecialBP_NonUniform4thOrderANFFixedBlock (NewY, modifiedt, initialFreq1, gamma, xi, initialFreq1, xi,  0);
% [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq1, gamma, xi);
% NewY = Y - 2*xi*X2./X3(N);                                   % Cascaded
% % [X21, X22,X23]= Special4thOderBP_NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq2, gamma, xi, initialFreq1, 0.25);
% 
% [X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq2, gamma, xi);
% 
% 
%  figure, plot(modifiedt, changingfreqs(1,:), 'g', modifiedt, changingfreqs(2,:), 'b',  modifiedt,X3/(2*pi), '-.k', modifiedt,X23/(2*pi), '-.r', 'LineWidth', linewidth)
%  axis([a b 60 70])
% % title ('0 dB')
% Ax = legend('1^{st} true frequency','2^{nd} true frequency','1^{st} frequency estimate','2^{nd} frequency estimate' )
% Ax.FontSize = fontSize;
% %legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
% ylabel('Frequency (Hz)')
% xlabel('Time (seconds)')
% grid on