clear;
clc;
close all;
a = 0;
b = 10;
T = .001;
tr=a:T:b;

f = 60;
f_initial = 66;   % 10% error in initialization

initialvalues=[0.5 .5 .5 (2*f_initial*pi)^2];
[t,Xia]=ode45('FunctionXiaAdaptive',tr,initialvalues);

initialvalues=[0.5 .5 .5 (2*f_initial*pi)^2];
[t,Marino]=ode45('FunctionMarinoAdaptive',tr,initialvalues);


trueFreq  =2*60*pi; 
freqFactor = 2*pi;

gamma = 0.001;
xi  = .15; 
A = 1;

y = A*sin(trueFreq*tr + pi/3);     % Pure sinusoid y(t)

% SNRdB = 5;
% y = y - A + awgn(zeros(size(tr))+ A ,SNRdB,'measured')  ;      % No Additive noise

[H1, H2,H3]= NonUniformThirdOrderANFFixedBlock(y, tr, 2*pi*(f_initial), gamma, xi);

[M1,M2,M3] = ThirdOrderMojiriANFBlock( y, T, 2*pi*(f_initial), gamma, xi+.05 );


figure, plot(tr,trueFreq*ones(size(tr))/freqFactor, 'g',tr, sqrt(Xia(:,4))/(freqFactor),'-.m',tr, sqrt(Marino(:,4))/(freqFactor),'--b',tr,M3/freqFactor , '-.k', tr,H3/freqFactor , '-.r','LineWidth', 2.5 )
% title ('Freq Estimation') 
xlabel('Time (s)')  
ylabel('Frequency (Hz)') 
Ax = legend('True frequency value','Xia Observer', 'Marino-Tomei Observer','Mojiri ANF', 'Proposed Discrete Method') 
Ax.FontSize = 14;
grid on  
  
