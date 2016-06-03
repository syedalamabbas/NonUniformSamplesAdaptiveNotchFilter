clear;
clc;
close all;

%------------------------Few Initializations------------------------------%
T = 0.001;                              % Sampling Rate
a = 0;                                  % t minimum
b = 10;                                 % t maximum
t= [a:T:b];                             % Uniformly spaced
%-------------------------------------------------------------------------%


%--------------------------ANF parameters---------------------------------%
gamma = 0.001;
xi  = .15; 
%-------------------------------------------------------------------------%

%------------------------------Signal parameters--------------------------%
A =  1;
phi = pi/2;

stable = ( (gamma / (4*xi)) < 1)

f =  450;
f_Est = zeros(size(f));
 
SNRdB = 5;
Y_Uniform = A * sin (2*pi*f*t + phi ) - A + awgn(zeros(size(t))+A,SNRdB,'measured');  % infecting the signal with noise
%-------------------------------------------------------------------------%


%----------------------------Modulation-----------------------------------%
Deltaf =  40;
f_ci = f + Deltaf;  % Modulator Frequency

Y_Sine  = Y_Uniform .* sin(2*pi*f_ci* t);           % Sine Modulate
Y_Cosine  = Y_Uniform .* cos(2*pi*f_ci* t);         % Cosine Modulate
%-------------------------------------------------------------------------%


%------------------------Low Passing Stage--------------------------------%
Fs = 1/T; % 1 kHz sampling frequency
d=fdesign.lowpass('Fp,Fst,Ap,Ast', 50, 75,1,90, Fs); 
designmethods(d);
Hd = design(d,'equiripple');
h = fvtool(Hd);   % get the handle to the graphics
set(h,'DesignMask','off'); % Turn off design mask
hchildren = get(h,'children');
haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
hline = get(haxes,'children');
set(hline,'LineWidth', 1.8)
legend('lowpass filter')
title ('(a)')


% figure, plot (t,Y_Sine)                       % Can check before
hold on
Y_Sine = filter(Hd,Y_Sine);                     % Apply Lowpass 1 stage
hold off
% figure, plot (t,Y_Sine)                       % Can check after

% figure, plot (t,Y_Cosine)                     % Can check before
hold on
Y_Cosine = filter(Hd,Y_Cosine);                 % Apply Lowpass 2 stage
hold off
% figure, plot (t,Y_Cosine)                     % Can check after
%-------------------------------------------------------------------------%


%--------------------------Amplitude Estimation---------------------------%
A_Estimated = sum( sqrt( (2*Y_Sine).^2 + (2* Y_Cosine).^2) )/ length(t);  
%-------------------------------------------------------------------------% 

%-------------------------Boosting Amplitude------------------------------%
Y_Sine =  Y_Sine / A_Estimated * 1.4;
%-------------------------------------------------------------------------%

%----------------------Final Offset Estimation----------------------------%
[H1, H2,H3]= FourthOrderANFFixedBlock(Y_Sine, T, 2*pi*(Deltaf+15), gamma, xi);
%-------------------------------------------------------------------------%


%------------------------ Plotting Solution --------------------------------------%
figure, plot(t,(f_ci-f)*ones(size(t)),'g', 'LineWidth', 2.3)
title ('(b)')
hold on
plot (t, H3/(2*pi),'r', 'LineWidth', 2.3)
hold off
str1 = 'Estimated frequency offset';
% str1 = strcat('Estimated Uniform,','N=', num2str(length(t)) );
% str2 = strcat('Estimated Nonuniform,','N=', num2str(length(modifiedt)) );
Ax = legend('True frequency offset',str1,  'Location','North')
Ax.FontSize = 14;
% leg = findobj(Ax,'type','text')
% set(leg,'fontsize', 18)
ylabel ('Frequency (Hz)')
xlabel ('time(s)')
grid on  
%-------------------------------------------------------------------------%