close all;
clear all;
clc;

%----------------Few Initializations-----------------------------------%
% T = 0.0001;                   % Sampling interval high

a = 0;                          % t minimum
b = 1;                          % t maximum 


SamplingRatesT = [.00001:.00001:.001];
gamma = 0.001;
xi  = .15; 
A =  1;
phi = pi/2;
SNRdB = 5;
f =  170* ones(size(SamplingRatesT)); %Hz

f_Est = zeros(size(SamplingRatesT));
f_4Est = zeros(size(SamplingRatesT));
f_2Est = zeros(size(SamplingRatesT));

for i = 1: length(SamplingRatesT)
    T = SamplingRatesT(i);          % New Sampling Rate
    t= [a:T:b];                     % Uniform spacing
    N = length(t);
    Y = A * sin (2*pi*f(i)*t + phi ) - 1+ awgn(zeros(size(t))+1,SNRdB,'measured');     % infecting the signal with noise
    flag = (rand() < .5) ;
    if (flag == 0 )
        freqError = +10;
    else
        freqError = -10;
    end
    
    [H1, H2,H3]=NonUniformSecondOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);         % 2nd order Estimation
    f_2Est (i) = H3(N)/(2*pi);
    
    [X1, X2,X13]=NonUniformThirdOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);        % 3rd order Estimation
    f_Est (i) = X13(N)/(2*pi);
    
    [B1, B2,B3]= NonUniform4thOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);        % 4th order Estimation
    f_4Est (i) = B3(N)/(2*pi);
end

figure, 
semilogx (SamplingRatesT,abs(f_2Est-f),SamplingRatesT,abs(f_Est-f),SamplingRatesT,abs(f_4Est-f), 'Linewidth', 1.5 )
xlabel('Sampling Rate T (sec)\rightarrow ')
ylabel ('Absolute error in frequency  (Hz)')
legend( '2^{nd} Order Approximation','3^{rd} Order Approximation','4^{th} Order Approximation')
grid on
axis tight



% T = 0.001;                      % Sampling interval low
% 
% 
% t= [a:T:b];                     % Uniform spacing
% N = length(t);
% % jitter=(2*rand(1, N)-1)*T/4; %jitter in [-T/3, T/3]
% % t=t+jitter;
% 
% % jitter2=(2*rand(1, N)-1)*T/8; %jitter in [-T/3, T/3]           % Generating Jitter
% % t=t+jitter2;
% %----------------------------------------------------------------------%
% 
% %-------------------------Specified parameters for the ANF----------------%
% gamma = 0.001;
% xi  = .15; 
% %-------------------------------------------------------------------------%
% 
% %--------------------------Speficied Parameters for the signal------------%
% A =  1;
% phi = pi/2;
% 
% stable = ( (gamma / (4*xi)) < 1)            % Checking the stability criteria
% f =  [2.5:5:1000];                           % Set of Frequencies to be evaluated
% K = length(f);
% f_Est = zeros(size(f));
% f_4Est = zeros(size(f));
% f_2Est = zeros(size(f));
% 
% % SNRdB = 20;
% SNRdB = 5;
% %-------------------------------------------------------------------------%
% 
% 
% %----------------------------Processing-----------------------------------%
% 
% for i =1: K
%     Y = A * sin (2*pi*f(i)*t + phi ) - 1+ awgn(zeros(size(t))+1,SNRdB,'measured');     % infecting the signal with noise
%     flag = (rand() < .5) ;
%     if (flag == 0 )
%         freqError = +10;
%     else
%         freqError = -10;
%     end
%     
%     [H1, H2,H3]=NonUniformSecondOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);         % 2nd order Estimation
%     f_2Est (i) = H3(N)/(2*pi);
%     
%     [X1, X2,X13]=NonUniformThirdOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);        % 3rd order Estimation
%     f_Est (i) = X13(N)/(2*pi);
%     
%     [B1, B2,B3]= NonUniform4thOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);        % 4th order Estimation
%     f_4Est (i) = B3(N)/(2*pi);
% end 
% 
% figure, 
% grid on
% hold on
% scatter(f, abs(f_2Est-f),  'MarkerFaceColor','m')
% scatter(f, abs(f_Est-f),  'MarkerFaceColor','g')
% scatter(f, abs(f_4Est-f),'d', 'MarkerFaceColor','b')
% hold off
% xlabel('Frequency Value (Hz)')
% ylabel ('Absolute error in frequency  (Hz)')
% legend( '2^{nd} Order Approximation','3^{rd} Order Approximation','4^{th} Order Approximation')
% axis ([0 250 0 30])
% % axis ([0 1000 0 30])
% %  print -depsc2 OrderComparison  
% %-------------------------------------------------------------------------%

