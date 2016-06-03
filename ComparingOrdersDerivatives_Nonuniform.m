 close all;
clear all;
clc;

%----------------Few Initializations-----------------------------------%
%  T = 0.0001;                   % Sampling interval high
  T = 0.001;                      % Sampling interval low
% a = 0;                          % t minimum
% b = 1;                          % t maximum 
fontSize = 14;
% t= [a:T:b];                     % Uniform spacing
% N = length(t);
% jitter=(2*rand(1, N)-1)*T/4; %jitter in [-T/3, T/3]
% t=t+jitter;

% jitter2=(2*rand(1, N)-1)*T/8; %jitter in [-T/3, T/3]           % Generating Jitter
% t=t+jitter2;
%----------------------------------------------------------------------%




%-------------------------Specified parameters for the ANF----------------%
gamma = 0.001;
xi  = .15; 
%-------------------------------------------------------------------------%

%--------------------------Speficied Parameters for the signal------------%
A =  1;
phi = pi/2;

stable = ( (gamma / (4*xi)) < 1)            % Checking the stability criteria
f =  [2.5:5:1000];                           % Set of Frequencies to be evaluated
K = length(f);
f_Est = zeros(size(f));
f_4Est = zeros(size(f));
f_5Est = zeros(size(f));
f_2Est = zeros(size(f));

N = K*5*2;
t = zeros (1, N);

%----------------------------------Introducing Jitter---------------%
min = -5;                   % Min value
max = 5;                    % Max value  
r = (max-min).*rand(1,N) + min;  % Randomly distributed Values

for i = 2:N                               % Jittered model with additive non-uniform sampling
   t(i) = t(i-1) + T+ r(i)*T*10^-1;
end
% SNRdB = 20;
SNRdB = 5;
%-------------------------------------------------------------------------%

%----------------------------Processing-----------------------------------%
for i =1: K
    noise = awgn(zeros(size(t))+1,SNRdB,'measured');
    Y = A * sin (2*pi*f(i)*t + phi ) - 1 + noise;     % infecting the signal with noise
    flag = (rand() < .5) ;
    if (flag == 0 )
        freqError = +10;
    else
        freqError = -10;
    end
    
    initialFreq = 2*pi*(f(i)+ freqError);
    
    [H1, H2,H3]=NonUniformSecondOrderANFFixedBlock(Y, t, initialFreq, gamma, xi);         % 2nd order Estimation
    f_2Est (i) = H3(N)/(2*pi);
    
    [X1, X2,X13]=NonUniformThirdOrderANFFixedBlock(Y, t, initialFreq, gamma, xi);        % 3rd order Estimation
    f_Est (i) = X13(N)/(2*pi);
    
   
%      [ B1,B2,B3] = SpecialBP_NonUniform4thOrderANFFixedBlock (Y, t, initialFreq, gamma, xi, initialFreq, xi + .2, 0);
%      [ B1,B2,B3] = Special4thOderBP_NonUniform4thOrderANFFixedBlock(Y, t, initialFreq, gamma, xi, initialFreq, xi+.5);
[B1, B2,B3]= NonUniform4thOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);        % 4th order Estimation
%       [B1, B2,B3]= SpecialComplementary_NonUniform4thOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);        % 4th order Estimation
    f_4Est (i) = B3(N)/(2*pi);
    
%     if (i == 50)
%      Y_derivative = A *2*pi*f(i)* cos (2*pi*f(i)*t + phi );
%      Y_Estimated_Derivatives = B1 .* B3 *(-2*xi);
%      figure, plot (t, Y_derivative, t, Y_Estimated_Derivatives)
%      legend ('True', 'Approximation')
%     end

%     [F1, F2,F3]= NonUniform5thOrderANFFixedBlock(Y, t, 2*pi*(f(i)+ freqError), gamma, xi);        % 5th order Estimation
%     f_5Est (i) = F3(N)/(2*pi);
end 

figure, 
grid on
hold on
scatter(f, abs(f_2Est-f),  'MarkerFaceColor','m')
scatter(f, abs(f_Est-f),  'MarkerFaceColor','g')
scatter(f, abs(f_4Est-f),'d', 'MarkerFaceColor','b', 'MarkerEdgeColor','b')
%  scatter(f, abs(f_5Est-f),'s', 'MarkerFaceColor','r')
hold off
xlabel('Frequency Value (Hz)')
ylabel ('Absolute error in frequency  (Hz)')
Ax = legend( '2^{nd} Order Approximation','3^{rd} Order Approximation','4^{th} Order Approximation' )
Ax.FontSize = fontSize;
% ,'5^{th} Order Approximation')
%   axis ([0 700 0 160])
 axis ([0 250 0 30])
% axis ([0 1000 0 30])
%  print -depsc2 OrderComparison  
%-------------------------------------------------------------------------%

