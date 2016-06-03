close all;
clear;
clc;

%------------------------------Few Initializations------------------------%
T = 0.001;
a = 0;
b = 1;
t= [a:T:b];                     % Uniformly spaced
N = length(t);
fontSize = 14;
%-------------------------------------------------------------------------%


%--------------------------------ANF Parameters---------------------------%
gamma = 0.001;
xi  = .15; 
%-------------------------------------------------------------------------%

%--------------------------Speficied Parameters for the signal------------%
A =  1;
phi = pi/2;
stable = ( (gamma / (4*xi)) < 1)
%  f =  60;
 f =  170;
f_Est = zeros(size(f));
SNRdB = 5;

Y_Uniform = A * sin (2*pi*f*t + phi ) - A + awgn(zeros(size(t))+A ,SNRdB,'measured');   % infecting the signal with noise
%-------------------------------------------------------------------------%

%-----------------------Processing----------------------------------------%
   freqError = +17;  % 10%
%   freqError = +6;  % 10%

[H1, H2,H3]= NonUniform4thOrderANFFixedBlock(Y_Uniform, t, 2*pi*(f+ freqError), gamma, xi);

%-------------------------------------------------------------------------%

%Initial Values

f_true = ones(size(t))*f;
X0 = zeros(1,3);
X0(3) = 2*pi*(f+ freqError);



%--------------------Direct R-K Method--------------------------%
SolutionY = RungeKutta4thOrderANF (t, X0,f,SNRdB, xi, gamma);
%---------------------------------------------------------------%

%--------------------Adomian Decomposition Method------------------------%
SolutionADM = FunctionalANFAdomianSpecial (t,X0, Y_Uniform, xi, gamma);
%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
% t, SolutionADM(:,3)/(2*pi), '-.', 'ADM',
figure, plot(t,f_true, 'g',t, SolutionY(:,3)/(2*pi),  t, H3/(2*pi),'--r','LineWidth' , 2.2 )
title ('(a)')
Ax = legend('True Value','Direct R-K Method', 'Proposed Method')
% Ax = legend('true value','MATLAB ode45(Hsu ANF)','Proposed Discrete Filter','Location','northeast' )
leg = findobj(Ax,'type','text')
Ax.FontSize = fontSize;
% set(leg,'fontsize', 14)

% axis([a b 170 200]) 
xlabel('Time (sec)')
% ylabel('Solution, \theta (t)/(2\pi)')
grid on
axis tight


%  t, abs(f_true-SolutionADM(:,3)'/(2*pi))
% figure, semilogy(t, abs(f_true-SolutionY(:,3)'/(2*pi)), t, abs(f_true-H3/(2*pi)), 'LineWidth' , 2.2)
% grid on
% legend('Direct R-K Method', 'Proposed Method')
% xlabel('Time (sec)')
% ylabel('Absolute Error in solution, \theta(t)/(2\pi).')
% axis tight
%-------------------------------------------------------------------------%

% print -depsc2 Convergence1
% print -depsc2 Convergence2
% print -depsc2 Convergence3

