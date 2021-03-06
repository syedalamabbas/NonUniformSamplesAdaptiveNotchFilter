clear;
clc;
close all;

%-------------------------Few Initializations --------------------------------%
T = 0.001;                                  % Sampling interval
a = 0;                                      % t minimum
b = 1;                                      % t maximum
t= [a:T:b];                                 % Uniformly spaced
%-----------------------------------------------------------------------%


%--------------------------ANF Parameter----------------------------------%
gamma = 0.001;
xi  = .15; 
%-------------------------------------------------------------------------%


%----------------------------Signal Parameter-----------------------------%
A1 = 1;
A2 =  2; 
A3 =  4;
phi = pi/2;

f =  60;
f_Est = zeros(size(f));
 
SNRdB = 5;

y2 =  A2 * sin (2*pi*f*t + phi ) - A2 + awgn(zeros(size(t))+A2,SNRdB,'measured');    % infecting the signal with noise proportional to the amplitude
y1 =  A1 * sin (2*pi*f*t + phi ) - A1 + awgn(zeros(size(t))+A1,SNRdB,'measured');    % infecting the signal with noise proportional to the amplitude
y3 =  A3 * sin (2*pi*f*t + phi ) - A3 + awgn(zeros(size(t))+A3,SNRdB,'measured');    % infecting the signal with noise proportional to the amplitude

figure, plot (t, A2 * sin (2*pi*f*t + phi ), t, y2, 'LineWidth', 1.5)
xlabel('Time(sec)')
ylabel('Amplitude(units)')

axis([0.4 0.5  -7.5 7.5 ])
% legend ('SNR = 5 dB')
legend ('Pure Signal', strcat('Noisy Signal, SNR =', num2str(SNRdB), 'dB'))
grid on

freq_offset = 18;
[H1, H2,H13]= FourthOrderANFFixedBlock(y1, T, 2*pi*(f+freq_offset), gamma, xi);


freq_offset = 18;
[H1, H2,H23]= FourthOrderANFFixedBlock(y2, T, 2*pi*(f+freq_offset), gamma, xi);


freq_offset = 18;
[H1, H2,H33]= FourthOrderANFFixedBlock(y3, T, 2*pi*(f+freq_offset), gamma, xi);


figure , plot3 (H1, H2, H23/ (2*pi), 'LineWidth', 2.3)
xlabel ('x_1')
ylabel ('x_2')
zlabel ('\theta (Hz)')
grid on  

hold on
xImage = [-40 -40; 40 40];   %# The x data for the image corners
yImage = [0 0; 0 0];             %# The y data for the image corners
zImage = [50 90; 50 90 ];   %# The z data for the image corners
surf(xImage,yImage,zImage, 'LineWidth', 2.5);
xlabel ('x_1')
ylabel ('x_2')
zlabel ('\theta (in Hz)')
colormap hsv 
alpha(.2)
hold off
title('(a)')

Indexes = find(H2 < 10^-3);     % locating point lying on the plane x_2 = 0 
 
figure, scatter (H1(Indexes ),H23 (Indexes)/ (2*pi) , '.')
 greater58 =  find(58 < H23 /(2*pi));
 lesser62 = find(62 > H23 /(2*pi));
 markerIndexes = intersect(greater58,lesser62);
 hold on
 scatter (H1(markerIndexes ),H23 (markerIndexes)/ (2*pi) , '.r')
 hold off
% title('Poincare Section')
title('(b)')
xlabel ('x_1')
ylabel ('\theta')
grid on
 
figure, plot (t, f*ones(size(t)), t, H13/ (2*pi), t, H23/ (2*pi), t, H33/ (2*pi), 'LineWidth', 2.5)
title('(a)')
legend ('True Value', 'Estimate, A = 1', 'Estimate, A = 2', 'Estimate, A = 4')
ylabel ('Frequency (Hz)')
xlabel ('time(s)')
grid on  

