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


%------------------------------Creating Jittered Sampling---------------%
Non_Uni_Samples = 940;
modifiedt = zeros (1, Non_Uni_Samples);
for i = 2:Non_Uni_Samples                               % Jittered model with additive non-uniform sampling
   modifiedt(i) = modifiedt(i-1) + T+ rand()*10^-4;
end
figure, 
hold on
stem(t,ones(size(t)), 'b','LineWidth',2.3)
stem(modifiedt,ones(size(modifiedt)),'r','LineWidth',2.3)
axis([0.02 0.04  -.5 1.5])
Ax = legend( 'Uniform Sampling', 'Additive Nonuniform Sampling')
leg = findobj(Ax,'type','text')
set(leg,'fontsize', 18)
grid on
hold off
%------------------------------------------------------------------------%


%---------------------Prep 2 sinusoids------------------------------%

Y = zeros(size(modifiedt));
FirstFrequencySignal = zeros(size(modifiedt));
SecondFrequencySignal = zeros(size(modifiedt));
changingfreqs = zeros(2,size(modifiedt));
A1 = 1;(4/pi);
phi1 = pi/3;
A2 = 1;(4/pi);
phi2 = pi/7;

cont_changingfreqs = zeros(2,size(modifiedt));

for i =1:length(t)
    if (t(i) <= .33)
        cont_changingfreqs(1,i)  = 72;
        cont_changingfreqs(2,i)  = 130;
    end
    if (t(i) > .33 && t(i) <= .66)
        cont_changingfreqs(1,i)  = 57;
        cont_changingfreqs(2,i)  = 150;
    end
    if (t(i) > .66)
        cont_changingfreqs(1,i)  = 80;
        cont_changingfreqs(2,i)  = 140;
    end
end

for i =1:length(modifiedt)
    if (modifiedt(i) <= .33)
        changingfreqs(1,i)  = 72;
        changingfreqs(2,i)  = 130;
    end
    if (modifiedt(i) > .33 && modifiedt(i) <= .66)
        changingfreqs(1,i)  = 57;
        changingfreqs(2,i)  = 150;
    end
    if (modifiedt(i) > .66)
        changingfreqs(1,i)  = 80;
        changingfreqs(2,i)  = 140;
    end
    FirstFrequencySignal(i) =  A1 * sin (2*pi*changingfreqs(1,i)*modifiedt(i)+ phi1);
    SecondFrequencySignal(i) =  A2 * sin (2*pi*changingfreqs(2,i)*modifiedt(i)+ phi2);
    Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i); % Composite Excitation Signal
end

SNRdB = 20;
noise = - A1 + awgn(zeros(size(modifiedt))+A1,SNRdB,'measured');
Y =Y + noise;   % infecting the signal with additive noise
figure,plot(modifiedt, Y, modifiedt,FirstFrequencySignal, modifiedt,SecondFrequencySignal,modifiedt, noise, '-.', 'LineWidth', linewidth)
axis([0.5-.03 0.5+.03 -3  3])
ylabel ('Amplitude (units)')
xlabel('Time(seconds)')
noise = strcat('Noise,SNR=', num2str(SNRdB), ' dB')
legend('Composite','1^{st} Sinusoid','2^{nd} Sinusoid',noise)
grid on

%-----------------------------------------------------------------------%

%---------------------Preparing the signal with 2 sinusoids---------------%

% A1 = 1;(4/pi);
% phi1 = pi/3;
% A2 = 1;(4/pi);
% phi2 = pi/7;
% freq1 = 72;
% freq2 = 120;
% 
% FirstFrequencySignal = A1* sin( 2*pi*freq1* modifiedt + phi1);
% SecondFrequencySignal = A2* sin( 2*pi*freq2* modifiedt + phi2);
% Y = FirstFrequencySignal + SecondFrequencySignal ;
% 
% 
% SNRdB = 0;
% noise = - A1 + awgn(zeros(size(modifiedt))+A1,SNRdB,'measured');
% Y =Y + noise;   % infecting the signal with additive noise
% figure,plot(modifiedt, Y,modifiedt,FirstFrequencySignal, modifiedt,SecondFrequencySignal,modifiedt, noise, '-.', 'LineWidth', linewidth)
% axis([0.5-.03 0.5+.03 -3  3])
% ylabel ('Amplitude (units)')
% xlabel('Time(seconds)')
% noise = strcat('Noise,SNR=', num2str(SNRdB), ' dB')
% legend('Composite','1^{st} Sinusoid','2^{nd} Sinusoid',noise)
% grid on
%------------------------------------------------------------------------%

%---------------------------- ANF parameters -----------------------------%
gamma = 0.001;
xi  = .15;

gamma = 0.01;
xi  = .2;
%------------------------------------------------------------------------%

initialFreq1 =  2*pi*60;
initialFreq2 =  2*pi*132;

%----------------------------------Cascade Configuration------------------%
NewY = Y;
[X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq1, gamma, xi);
NewY = Y - 2*xi*X2./X3;                                   % Cascaded
[X21, X22,X23]= NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq2, gamma, xi);

figure,
stem(modifiedt,X23/(2*pi),'r', 'LineWidth', linewidth-.5)
hold on
stem(t,cont_changingfreqs(2,:), 'g', 'LineWidth', linewidth-.5)
hold off
axis([.5 .52 149.5 151.5])
% title ('Cascade Configuration')
legend('Estimate 2', 'True 2^{nd} Freq' )
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on

 
figure, 
plot(modifiedt, changingfreqs(1,:),modifiedt, changingfreqs(2,:),  modifiedt,X3/(2*pi), '-.k', modifiedt,X23/(2*pi), '-.r', 'LineWidth', linewidth)

% scatter(modifiedt,changingfreqs(1,:))
% scatter(modifiedt,X3/(2*pi))
% hold off

% axis([a b changingfreqs(1,1)-15 changingfreqs(2,1)+15])
axis([a b 50 200])
% title ('Cascade Configuration')
legend('True 1^{st} Freq','True 2^{nd} Freq','Estimate 1','Estimate 2' )
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on
%-------------------------------------------------------------------------%