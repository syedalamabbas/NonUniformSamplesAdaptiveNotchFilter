 close all;
clear all;
clc;

%----------------Few Initializations-----------------------------------%
T = 0.001;                      % Sampling interval
a = 0;                          % t minimum
b = 1;                          % t maximum        
t_cont =  [a:T*.001:b];         % Dense 
t= [a:T:b];                     % Uniformly spaced
fontSize = 14;
lineWidth = 2;
%----------------------------------------------------------------------%


%------------------------------Creating Jittered Sampling---------------%
Non_Uni_Samples = length(t);
modifiedt = zeros (1, Non_Uni_Samples);

min = -5;                   % Min value
max = 5;                    % Max value  
r = (max-min).*rand(1,Non_Uni_Samples) + min;  % Randomly distributed Values

for i = 2:Non_Uni_Samples                               % Jittered model with additive non-uniform sampling
   modifiedt(i) = modifiedt(i-1) + T + r(i)*T*10^-1;
end
 
 
% figure, 
% hold on
% stem(t,ones(size(t)), 'b','LineWidth',2)
% stem(modifiedt,ones(size(modifiedt)),'r','LineWidth',2)
% axis([0.02 0.04  -.5 1.5])
% Ax = legend( 'Uniform Sampling', 'Additive Nonuniform Sampling')
% leg = findobj(Ax,'type','text')
% set(leg,'fontsize', fontSize)
% grid on
% hold off
%------------------------------------------------------------------------%


%-------------------------Specified parameters for the ANF----------------%
gamma = 0.001;
xi  = .15; 
%-------------------------------------------------------------------------%

%--------------------------Speficied Parameters for the signal----------%
N = length(t);
A =  1;
phi = pi/2;

stable = ( (gamma / (4*xi)) < 1)
%     f =  60;
 f =  170;
f_Est = zeros(size(f));
 
SNRdB = 5;
%------------------------------------------------------------------------%

%-----------------------Generating the Signals ---------------------------%
Y_Uniform = A * sin (2*pi*f*t + phi ) - 1+ awgn(zeros(size(t))+1,SNRdB,'measured');                        % infecting the signal with additive noise
Y_NonUniform = A * sin (2*pi*f*modifiedt + phi ) -1 + awgn(zeros(size(modifiedt))+1,SNRdB,'measured');     % infecting the signal with additive noise

% Y_NonUniform = A* sin (2*pi*f*modifiedt + phi )  + A * sin (2*pi*230*modifiedt + phi )+ A * sin (2*pi*207.8*modifiedt + phi )+ + A * sin (2*pi*507.8*modifiedt + phi )+ + A * sin (2*pi*407.8*modifiedt + phi ) -5 + awgn(zeros(size(modifiedt))+5,SNRdB,'measured');     % infecting the signal with additive noise
% 
% - 1+ awgn(zeros(size(modifiedt))+1,SNRdB,'measured'); 

% + A * sin (2*pi*230*modifiedt + phi )+ A * sin (2*pi*207.8*modifiedt + phi )+ + A * sin (2*pi*507.8*modifiedt + phi )+ + A * sin (2*pi*407.8*modifiedt + phi ) -5 + awgn(zeros(size(modifiedt))+5,SNRdB,'measured');     % infecting the signal with additive noise


Y_Continuous =  A * sin (2*pi*f*t_cont + phi );



figure, 
hold on
plot (t_cont,Y_Continuous, 'g','LineWidth',lineWidth)
stem(t,Y_Uniform, 'b','LineWidth',2)
stem(modifiedt,Y_NonUniform,'r','LineWidth',lineWidth)
hold off
axis([0.02 0.04  -2 2])
Ax = legend( 'Continuous (Noiseless)','Uniform (SNR = 5dB)', 'Nonuniform (SNR = 5dB)')
Ax.FontSize = fontSize;
% leg = findobj(Ax,'type','text')
% set(leg,'FontSize', fontSize) 
ylabel ('Amplitude (units)')
xlabel ('Time(s)')
grid on 
% print -depsc2 SignalUnifNonUnif
%-----------------------------------------------------------------------%

 
%------------------------Processing Frequency Estimation-----------------%
flag = (rand() < .5) ;
% if (flag == 0 )
       freqError = +17;
%         freqError = +6;
% else
%        freqError = -17;
%        freqError = -6;
% end
initialFreq = 2*pi*(f+ freqError);
%   initialFreq  = 2*pi*205; 

[S1, S2,S3]= FourthOrderANFFixedBlock(Y_Uniform, T, initialFreq, gamma, xi);                      % This is a uniformly sampled system

[H1, H2,H3]= NonUniform4thOrderANFFixedBlock(Y_NonUniform, modifiedt, initialFreq, gamma, xi);    % This is a non-uniformly sampled system

% %  bandPassFreq = initialFreq;
% %  xi_c = xi ; 
% [ M1,M2,M3] = SpecialBP_NonUniform4thOrderANFFixedBlock (Y_NonUniform, modifiedt, initialFreq, gamma, xi, initialFreq, xi+.2,  0);
% 
% [ X1,X2,X3] = Special4thOderBP_NonUniform4thOrderANFFixedBlock (Y_NonUniform, modifiedt, initialFreq, gamma, xi, initialFreq, xi+.2);
% 
% [G1,G2,G3] = SpecialComplementary_NonUniform4thOrderANFFixedBlock(Y_NonUniform, modifiedt, initialFreq, gamma, xi); 


% figure, plot (modifiedt,Y_NonUniform, modifiedt,BP1)
% legend('True/Desired', 'BP output')
 
figure, plot(t_cont,f*ones(size(t_cont)),'g', 'LineWidth', lineWidth)
title('(a)')
hold on
plot (t, S3/(2*pi),'b' ,'LineWidth', lineWidth)
plot (modifiedt, H3/(2*pi),'r', 'LineWidth', lineWidth)
% plot (modifiedt, M3/(2*pi),'k', 'LineWidth', lineWidth) 
% plot (modifiedt, X3/(2*pi),'m', 'LineWidth', lineWidth)
% plot (modifiedt, G3/(2*pi),'c', 'LineWidth', lineWidth)
hold off
str1 = 'Estimated Uniform';
str2 = 'Estimated Nonuniform';
% str3 = 'Estimated with BP Intergrated'
% str4 = 'Estimated with BP 4th order Intergrated'
% str5 = 'Estimated with Complementary Intergrated'
% Ax = legend('True Value',str1, str2,str3,str4, 'Location','North')
Ax = legend('True Value',str1, str2, 'Location','North')
Ax.FontSize = fontSize;
% leg = findobj(Ax,'type','text')
% set(leg,'fontsize', fontSize)
ylabel ('Frequency (Hz)')
xlabel ('Time(s)')
% axis([a b 155 190])
grid on  
axis tight

% print -depsc2 FreqConverUnifNonUnif    % Higher Freq initialization, 170
% print -depsc2 FreqConverUnifNonUnif1     % Lower Freq initialization

% print -depsc2 FreqConverUnifNonUnif2    % Higher Freq initialization, 60
% print -depsc2 FreqConverUnifNonUnif3     % Lower Freq initialization
%------------------------------------------------------------------------%


%------------Plotting the first derivative of the approximation--------%

%  Y_derivative = A *2*pi*f* cos (2*pi*f*modifiedt + phi );
%  Y_Estimated_Derivatives = H1 .* H3 *(-2*xi);
%  figure, plot (modifiedt, Y_derivative, modifiedt, Y_Estimated_Derivatives, 'LineWidth' , 2.3)
% Ax = legend ('True 1^{st} Derivative of y', 'Discrete ANF Approx, Dy')
% Ax.FontSize = fontSize;
%  axis([.4 .6 -600 600])
% % axis([.4 .5 -1200 1200])
% ylabel ('Amplitude (units)')
% xlabel ('Time(s)')
% grid on
