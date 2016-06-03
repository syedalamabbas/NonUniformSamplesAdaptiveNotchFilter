clear all;
close all;
clc;

%--------------------------Selecting tuning parameters externally---------%
gamma = .01;
xi = .15;


T = 0.001;
a = 0;
b = 1;
t= [a:T:b];
N = length(t);

fontSize = 14;


%-----------------------Jittered Model-----------------------%
Non_Uni_Samples = length(t);
modifiedt = zeros (1, Non_Uni_Samples);

min = -1;                   % Min value
max = 1;                    % Max value  
r = (max-min).*rand(1,Non_Uni_Samples) + min;  % Randomly distributed Values

for i = 2:Non_Uni_Samples                               % Jittered model with additive non-uniform sampling
   modifiedt(i) = modifiedt(i-1) + T+ r(i)*T*10^-1;
end
%---------------------------------------------------------------------%
modifiedt = t;     % Making it all uniform


Y = zeros(1,N);
FirstFrequencySignal = zeros(zeros(1,N));
SecondFrequencySignal = zeros(zeros(1,N));
ThirdFrequencySignal = zeros(zeros(1,N));
FourthFrequencySignal = zeros(zeros(1,N));
FifthFrequencySignal = zeros(zeros(1,N));

changingfreqs = zeros(5,N);
A1 = 1;(4/pi);
phi1 = pi/3;
A2 = .5* A1;
phi2 = pi/14;
A3 = .33* A1;
phi3 = pi/7;
A4 = .25* A1;
phi4 = pi/19;
A5 = .2* A1;
phi5 = pi/5;

fund = 72;
for i =1:length(t)
    
    if(t(i) > .33)
        fund = 60;
    end
    
     if(t(i) > .66)
        fund = 80;
    end
    
    changingfreqs(1,i)  = fund;
    changingfreqs(2,i)  = 2*fund;
    changingfreqs(3,i)  = 3*fund;
    changingfreqs(4,i)  = 4*fund;
    changingfreqs(5,i)  = 5*fund;
    FirstFrequencySignal(i) =  A1 * sin (2*pi*changingfreqs(1,i)*t(i)+ phi1);
    SecondFrequencySignal(i) =  A2 * sin (2*pi*changingfreqs(2,i)*t(i)+ phi2);
    ThirdFrequencySignal(i) =  A3 * sin (2*pi*changingfreqs(3,i)*t(i)+ phi3);
    FourthFrequencySignal(i) =  A4 * sin (2*pi*changingfreqs(4,i)*t(i)+ phi4);
    FifthFrequencySignal(i) =  A5 * sin (2*pi*changingfreqs(5,i)*t(i)+ phi5);
    
    Y(i) = FirstFrequencySignal(i) + SecondFrequencySignal(i)+ThirdFrequencySignal(i)+FourthFrequencySignal(i)+FifthFrequencySignal(i); % Composite Excitation Signal
end

SNRdB = 20;
PureSignal = Y;
Y =awgn(Y,SNRdB,'measured');   % infecting the signal with noise


figure,plot(t, Y,t,FirstFrequencySignal, t,SecondFrequencySignal,t,ThirdFrequencySignal, t,FourthFrequencySignal,t,FifthFrequencySignal ,t, (Y - PureSignal), 'LineWidth', 2.5)
title('(a)')
axis([0 0.035 -2.5  2.5])
grid on
ylabel ('Amplitude (units)')
xlabel('Time(seconds)')
noise  = strcat ('Noise, SNR=', num2str(SNRdB))
legend('Composite','Fundamental','2nd Harmonic','3rd Harmonic','4th Harmonic','5th Harmonic', noise)
% print -depsc2 FiveFrequencies_Noise
% 
 f_c = 74;
 initialFreq = 2*pi*f_c; 
% r = 0.87;
% omega0 = initialFreq*T;                  % Computing Notch at specific Freqs
% b = [1 -2*cos(omega0) 1];
% a = [1 -2*r*cos(omega0) r^2];
% % fvtool (b-a,a, 'Fs', 1/T);
% % str = strcat('Bandpass, f = ', num2str(f_c), 'Hz');
% % legend (str)
% % print -depsc2 BandPass_PreFilter
% 
% NewY= filter(b-a,a,Y);  
 NewY = Y; % Prefiltering
 [X1, X2,X3]= NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq, gamma, xi);
 
% [X1, X2,X3]= Special4thOderBP_NonUniform4thOrderANFFixedBlock (Y, modifiedt, initialFreq, gamma, xi, initialFreq, xi+.2);

% NonUniform4thOrderANFFixedBlock(NewY, modifiedt, initialFreq, gamma, xi);

% [X1, X2,X3]= FourthOrderANFFixedBlock(NewY, T, initialFreq, gamma, xi);
figure, plot3(X1,X2,X3/(2*pi), 'LineWidth', 1.4)
title('(b)')
% title('3D state space plot from Agent based ODE solver Hsu ANF 1')
grid on
xlabel('x_1') 
ylabel('x_2')
zlabel('\theta (Hz)')
axis tight
% print -depsc2 StateEvolutionFund

figure, plot(t, changingfreqs(1,:), modifiedt,X3/(2*pi), '-.r','LineWidth', 2.5)
title('(a)')
axis([a b 57 83]) 
Ax = legend('True Value','Estimated Value' )
Ax.FontSize = fontSize;
ylabel('Frequency (Hz)') 
xlabel('Time (seconds)')
grid on
% print -depsc2 FundFreqTracking

M = (changingfreqs(1,:)- X3/(2*pi))/N;                  % Mean Square Error
figure, semilogy( t, M.^2, '-.r','LineWidth', 2.5)
title('(b)')
% axis([a b 10^(-5) 10^3]) 
Ax = legend( 'Mean Squared Error' )
Ax.FontSize = fontSize
ylabel('Mean Squared Error in Frequency Estimates')
xlabel('Time (seconds)')
 grid on
% print -depsc2 EstimateFundError
