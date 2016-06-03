clear all;
close all;
clc;


%---------------------------Few initializations---------------------------%
T = 0.001;
a = 0;
b = 1;
t= [a:T:b];
%-------------------------------------------------------------------------%

%--------------------------Speficied Parameters for the signal----------%
N = length(t);
K = 7;
A =  .5 + (1)*rand(1,K);
phi = pi* rand(1,K);

f = [78, 137, 177, 234, 288, 325, 437];

PureSignal =  zeros (1,N);
Signals = zeros (K,N);
for i =1: K
    Signals(i,:) = A(i)*sin(2*pi*f(i)*t + phi(i));  
    PureSignal = PureSignal +  Signals(i,:);
end
SNRdB = 10;
Y =awgn(PureSignal,SNRdB,'measured');   % infecting the signal with noise


linewidth = 2.5;
figure,
hold on
plot(t, Y,'-.g','LineWidth',  linewidth-.5)
plot(t, Y - PureSignal, '-.', 'LineWidth', linewidth-1)
hold off
ylabel ('Amplitude (units)')
xlabel('Time(seconds)')
noise = strcat('Noise,SNR=', num2str(SNRdB), ' dB')
legend('Composite',noise)
grid on
% print -depsc2 SevenSignals_Noise

colors = ['b','g','r','c','m','y','k'];
figure, 
hold on
for i =1:K 
    plot (t,Signals(i,:), colors(i), 'LineWidth', linewidth);
end
hold off
axis ([.0 .03 -1.9 1.9])
str1 = strcat(num2str(f(1)), ' Hz');
str2 = strcat(num2str(f(2)), ' Hz');
str3 = strcat(num2str(f(3)), ' Hz');
str4 = strcat(num2str(f(4)), ' Hz');
str5 = strcat(num2str(f(5)), ' Hz');
str6 = strcat(num2str(f(6)), ' Hz');
str7 = strcat(num2str(f(7)), ' Hz');
ylabel ('Amplitude (units)')
xlabel('Time(seconds)')
legend( str1,str2,str3,str4,str5,str6,str7)
grid on 
% print -depsc2 DecomSevenSignals_Noise


figure, stem (f, A, 'LineWidth', 2.3)
xlabel('Frequency (Hz)')
ylabel ('Amplitude (units)')
grid on


initialFreqs = [70, 130, 180, 240, 270, 330, 450];     % Narrow BandPass center initial Frequencies

%----------------------Prefiltered Configuration -------------------------%

% Passband  +- 10% from initial center frequency %
% Stopband  +- 15% from initial center frequency %

Fs = 1/T; % 1 kHz sampling frequency as defined in the signal

d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' , initialFreqs(1)-.15* initialFreqs(1) , initialFreqs(1)-.1* initialFreqs(1) ,initialFreqs(1) + .12* initialFreqs(1), initialFreqs(1) + .18* initialFreqs(1), 90,2,90, Fs); 
% designmethods(d);
Hd1 = design(d,'equiripple');


d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' , initialFreqs(2)-.15* initialFreqs(2) , initialFreqs(2)-.1* initialFreqs(2) ,initialFreqs(2) + .1* initialFreqs(2), initialFreqs(2) + .15* initialFreqs(2), 90,2,90, Fs); 
Hd2 = design(d,'equiripple');

d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' , initialFreqs(3)-.15* initialFreqs(3) , initialFreqs(3)-.1* initialFreqs(3) ,initialFreqs(3) + .1* initialFreqs(3), initialFreqs(3) + .15* initialFreqs(3), 90,2,90, Fs); 
Hd3 = design(d,'equiripple');

d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' , initialFreqs(4)-.15* initialFreqs(4) , initialFreqs(4)-.1* initialFreqs(4) ,initialFreqs(4) + .1* initialFreqs(4), initialFreqs(4) + .15* initialFreqs(4), 90,2,90, Fs); 
Hd4 = design(d,'equiripple');

d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' , initialFreqs(5)-.1* initialFreqs(5) , initialFreqs(5)-.08* initialFreqs(5) ,initialFreqs(5) + .08* initialFreqs(5), initialFreqs(5) + .15* initialFreqs(5), 90,2,90, Fs); 
Hd5 = design(d,'equiripple');

d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' , initialFreqs(6)-.04* initialFreqs(6) , initialFreqs(6)-.02* initialFreqs(6) ,initialFreqs(6) + .02* initialFreqs(6), initialFreqs(6) + .04* initialFreqs(6), 90,2,90, Fs); 
Hd6 = design(d,'equiripple');

d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2' , initialFreqs(7)-.1* initialFreqs(7) , initialFreqs(7)-.04* initialFreqs(7) ,initialFreqs(7) + .04* initialFreqs(7), initialFreqs(7) + .1* initialFreqs(7), 90,2,90, Fs); 
Hd7 = design(d,'equiripple');

h = fvtool(Hd1,Hd2, Hd3, Hd4, Hd5, Hd6, Hd7 );   % get the handle to the graphics
set(h,'DesignMask','off'); % Turn off design mask
hchildren = get(h,'children');
haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
hline = get(haxes,'children');
set(hline,'LineWidth', 1.8)
str1 = strcat('Bandpass, f_1(center) = ', num2str(initialFreqs(1)), 'Hz');
str2 = strcat('Bandpass, f_2(center) = ', num2str(initialFreqs(2)), 'Hz');
str3 = strcat('Bandpass, f_3(center) = ', num2str(initialFreqs(3)), 'Hz');
str4 = strcat('Bandpass, f_4(center) = ', num2str(initialFreqs(4)), 'Hz');
str5 = strcat('Bandpass, f_5(center) = ', num2str(initialFreqs(5)), 'Hz');
str6 = strcat('Bandpass, f_6(center) = ', num2str(initialFreqs(6)), 'Hz');
str7 = strcat('Bandpass, f_7(center) = ', num2str(initialFreqs(7)), 'Hz');
legend (str1,str2, str3, str4, str5, str6, str7);


NewY= filter(Hd1,Y);                            % Prefiltering
% figure, plot(t, Y, t , NewY)
% xlabel('t(sec)')
% ylabel('A(units)')
% legend('Signal', 'Filtered in 70 Hz range')

gamma = 0.001;
xi  = .15;
[X1, X2,X13]= NonUniform4thOrderANFFixedBlock(NewY, t, 2*pi*initialFreqs(1), gamma, xi);


NewY= filter(Hd2,Y);                            % Prefiltering

gamma = 0.001;
xi  = .15;
[X1, X2,X23]= NonUniform4thOrderANFFixedBlock(NewY, t, 2*pi*initialFreqs(2), gamma, xi);


NewY= filter(Hd3,Y);                            % Prefiltering

gamma = 0.001;
xi  = .15;
[X1, X2,X33]= NonUniform4thOrderANFFixedBlock(NewY, t, 2*pi*initialFreqs(3), gamma, xi);

NewY= filter(Hd4,Y);                           % Prefiltering

gamma = 0.001;
xi  = .015;
[X1, X2,X43]= NonUniform4thOrderANFFixedBlock(NewY, t, 2*pi* initialFreqs(4), gamma, xi);


%----------------------Using High Frequency estimation Trick for the other 3 freq--------------%

%-----------------------------Estimating f(5)---------------------------%
% figure, plot(t, Y);
NewY= filter(Hd5,Y);                          % Prefiltering
% figure, plot(t, NewY); 

deltaf = 70;                                % select offset 
f_c1 = initialFreqs(5)+deltaf;
Y_Sine  = NewY .* sin(2*pi*f_c1* t);           % Sine Modulate

Y_Cosine  = NewY .* cos(2*pi*f_c1* t);         % Cosine Modulate


d=fdesign.lowpass('Fp,Fst,Ap,Ast', deltaf, deltaf+10,1,90, Fs); 
Hd = design(d,'equiripple');
h = fvtool(Hd);   % get the handle to the graphics
set(h,'DesignMask','off'); % Turn off design mask
hchildren = get(h,'children');
haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
hline = get(haxes,'children');
set(hline,'LineWidth', 1.8)
legend('lowpass filter')

% figure, plot(t, Y_Sine)
Y_Sine = filter(Hd,Y_Sine);                     % Apply Lowpass 1 stage
% figure, plot(t, Y_Sine)
Y_Cosine = filter(Hd,Y_Cosine);                 % Apply Lowpass 2 stage

A_Estimated = sum( sqrt( (2*Y_Sine).^2 + (2* Y_Cosine).^2) )/ length(t);  
Y_Sine =  Y_Sine / A_Estimated * 1.5;

% figure, plot(t,Y_Sine)

gamma = 0.001;
xi  = .15; 
[X1, X2,X53]= NonUniform4thOrderANFFixedBlock(Y_Sine, t, 2*pi*deltaf, gamma, xi);

%------------------------------------------------------------------------%


%--------------------------Estimating f(6)-----------------------------%
NewY= filter(Hd6,Y);                          % Prefiltering
% figure, plot(t, NewY); 

deltaf = 80;                                % select offset 
f_c2 = initialFreqs(6)+deltaf;
Y_Sine  = NewY .* sin(2*pi*f_c2* t);           % Sine Modulate

Y_Cosine  = NewY .* cos(2*pi*f_c2* t);         % Cosine Modulate


d=fdesign.lowpass('Fp,Fst,Ap,Ast', deltaf, deltaf+10,1,90, Fs); 
Hd = design(d,'equiripple');
h = fvtool(Hd);   % get the handle to the graphics
set(h,'DesignMask','off'); % Turn off design mask
hchildren = get(h,'children');
haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
hline = get(haxes,'children');
set(hline,'LineWidth', 1.8)
legend('lowpass filter')

% figure, plot(t, Y_Sine)
Y_Sine = filter(Hd,Y_Sine);                     % Apply Lowpass 1 stage
% figure, plot(t, Y_Sine)
Y_Cosine = filter(Hd,Y_Cosine);                 % Apply Lowpass 2 stage

A_Estimated = sum( sqrt( (2*Y_Sine).^2 + (2* Y_Cosine).^2) )/ length(t);  
Y_Sine =  Y_Sine / A_Estimated * 2;

% figure, plot(t,Y_Sine)

gamma = 0.001;
xi  = .15; 
[X1, X2,X63]= NonUniform4thOrderANFFixedBlock(Y_Sine, t, 2*pi*deltaf, gamma, xi);
%-------------------------------------------------------------------------%

%--------------------------Estimating f(7)-----------------------------%
% figure, plot(t, Y)                          % Observe Before
NewY= filter(Hd7,Y);                          % Prefiltering
% figure, plot(t, NewY)                       % Observe After


deltaf = 100;                                % select offset 
f_c3 = initialFreqs(7)+deltaf;
Y_Sine  = NewY .* sin(2*pi*f_c3* t);           % Sine Modulate

Y_Cosine  = NewY .* cos(2*pi*f_c3* t);         % Cosine Modulate


% d=fdesign.lowpass('Fp,Fst,Ap,Ast', deltaf, deltaf+10,1,90, Fs);   % Dont use lowpass because of aliased frequency about 12 Hz

% Hd = design(d,'equiripple');
% h = fvtool(Hd);   % get the handle to the graphics
% set(h,'DesignMask','off'); % Turn off design mask
% hchildren = get(h,'children');
% haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
% hline = get(haxes,'children');
% set(hline,'LineWidth', 1.8)
% legend('lowpass filter')
% 
% % figure, plot(t, Y_Sine)
% Y_Sine = filter(Hd,Y_Sine);                     % Apply Lowpass 1 stage
% % figure, plot(t, Y_Sine)
% Y_Cosine = filter(Hd,Y_Cosine);                 % Apply Lowpass 2 stage

A_Estimated = sum( sqrt( (2*Y_Sine).^2 + (2* Y_Cosine).^2) )/ length(t);  
Y_Sine =  Y_Sine / A_Estimated * 1.5;

% figure, plot(t,Y_Sine)

gamma = 0.001;
xi  = .015; 
[X1, X2,X73]= NonUniform4thOrderANFFixedBlock(Y_Sine, t, 2*pi*deltaf, gamma, xi);
%-------------------------------------------------------------------------%


unitStep = ones(1,N);
figure, plot(t, unitStep* f(1),t, unitStep* f(2),t, unitStep* f(3),t, unitStep* f(4),  t,X13/(2*pi), '-.', t,X23/(2*pi), '-.', t,X33/(2*pi), '-.',  t,X43/(2*pi), '-.', 'LineWidth', linewidth)
title('(a)')
% axis tight 
% title ('Prefiltered Configuration')
Ax = legend('True freq 1','True freq 2','True freq 3','True freq 4', 'Estimate 1','Estimate 2','Estimate 3','Estimate 4' )
Ax.FontSize = 14;
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on
% print -depsc2 EstimatesForFiveArbitrary
%-------------------------------------------------------------------------%


figure, plot (t, unitStep* (f_c1 - f(5)),t, unitStep* (f_c2 - f(6)),t, unitStep* (f_c3 - f(7)), t,X53/(2*pi), '-.', t,X63/(2*pi),'-.',t,X73/(2*pi),'-.',  'LineWidth', linewidth)
title('(b)')
ylabel ('Frequency (Hz)')
xlabel('Time (seconds)')
Ax = legend('True freq offset 1','True freq offset 2', 'True freq offset 3', 'Estimate offset 1','Estimate offset 2', 'Estimate offset 3')
Ax.FontSize = 14;
grid on
