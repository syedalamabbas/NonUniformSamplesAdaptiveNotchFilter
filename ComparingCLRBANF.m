clear all;
close all;
clc;

%-------------------------Few Initializations --------------------------------%
T = 0.001;                      % Sampling interval
a = 0;                          % t minimum
b = 1;                          % t maximum
t= [a+T/2:T:b];                 % Sampled interval
N = length(t);
jitter2=(2*rand(1, N)-1)*T/8;    %jitter in [-T/3, T/3]

% figure, 
% 
% hold on
% stem(t,ones(size(t)), 'b','LineWidth',2.3)
% t=t+jitter2;
% stem(t,ones(size(t)), 'g','LineWidth',2.3)
% hold off
%-----------------------------ANF parameters---------------------------%
gamma = 0.001;
xi  = 0.15;
%----------------------------------------------------------------------%

%-----------------------------CLRB-------------------------------------%
F = [1:1:250]*0.001;
M = length(F);
CLRB = zeros(1,M);
for i = 1:M
    temp1=0;
    temp2 = 0;
    temp3=0;
    temp4=0;
    for j=0:N-1
        temp1 = temp1 + (t(j+1))^2;
        temp2 = temp2 + (t(j+1))^2*(cos(4*pi*F(i)*t(j+1)));
        temp3 = temp3 + (t(j+1))^2*(sin(4*pi*F(i)*t(j+1)));
        temp4=temp4+sin(4*pi*F(i)*t(j+1));
    end
    Temp(i)=temp3;
    temp=temp1+(temp2*temp2+temp3*temp3)^(1/2);
    CLRB(i) = 1/temp;
end
CLRB=100*CLRB;
%----------------------------------------------------------------------%

%--------------------------Speficied Parameters for the signal----------%
A =  1;
phi = pi/2;

stable = ( (gamma / (4*xi)) < 1)
f =  [1:1:250];
K = length(f);
f_Est = zeros(size(f));
PureSignal =  zeros (1,N);
SNRdB =[20]; %[80]; % [20];
sig_square = 10^(-6);
E = 100;
MatF = zeros(K,E);
varF = zeros(K,1);

for m = 1:length(SNRdB)
    for i =1: K
        for j=1:E
            Y =A * sin (2*pi*f(i)*t + phi ) - 1 + awgn(zeros(1,N)+1,SNRdB(m),'measured');   % infecting the signal with noise
            flag = (rand() < .5) ;
            if (flag == 0 )
                freqError = +10;
            else
                freqError = -10;
            end
%             [X1, X2,X13]=SpecialComplementary_NonUniform4thOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi); 
            [X1, X2,X13]=NonUniform4thOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);
            f_Est (i) = X13(N)/(2*pi);
            MatF (i,j) = f_Est (i);
        end
        varF(i) = (2*pi)^2* var(MatF(i,:));
        meanF(i)=mean(MatF(i,:));
    end
end

figure
plot(F*1000, meanF-F*1000,'b', 'LineWidth', 1.8)
title ('(a)')
xlabel('Frequency Value (Hz)')
ylabel('Mean Value Error of Frequency Estimation')
grid on

figure,
plot(F*1000, log(abs(varF/CLRB'))/log(10),'b', 'LineWidth', 1.8)%, '-',  F*1000, 10*log(CLRB)/log(10), '--')
title ('(b)')
xlabel('Frequency Value (Hz)')
ylabel('Logarithm of Relative Frequency Variance Error') % log_{10) Var(\theta)')
grid on
% % print -depsc2 CLRBForN100
% %print -depsc2 CLRBForN500
%------------------------------------------------------------%