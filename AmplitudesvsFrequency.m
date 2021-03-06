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

sqrt(4* xi / gamma)
%----------------------------Signal Parameter-----------------------------%
%  A = 2;
f =  60;
freq_offset = 28;
phi = pi/2;
f_Est = zeros(size(f));

SNR = 5;

Amplitudes = [.1:.5: 500];
FreqError = ones(size(Amplitudes ));
%Amplitudes(i) * sin (2*pi*56*t + phi )  Harmonic pollution

for i=1:length(Amplitudes) 
    y =   Amplitudes(i) * sin (2*pi*f*t + phi ) - Amplitudes(i) + awgn(zeros(size(t))+Amplitudes(i),SNR,'measured');    % infecting the signal with noise proportional to the amplitude
    
    [H1, H2,H13]= FourthOrderANFFixedBlock(y, T, 2*pi*(f+freq_offset), gamma, xi);
    
    if (Amplitudes(i)==2)
        figure, plot (t, Amplitudes(i) * sin (2*pi*f*t + phi ), t, y, 'LineWidth', 1.5)
        xlabel('Time(sec)')
        ylabel('Amplitude(units)')
        axis([0.4 0.5  -7.5 7.5 ])
        legend ('Pure Signal', strcat('Noisy Signal, SNR =', num2str(SNR), 'dB'))
        grid on
        title('(a)')


        figure, plot (t, f*ones(size(t)), t, H13/ (2*pi),  'LineWidth', 2.5)
        legend ('True Value', strcat('Estimate, SNR= 5 dB, A =', num2str(Amplitudes(i))) )
        ylabel ('Frequency (Hz)')
        xlabel ('time(s)')
        grid on  
    end
    
    FreqError(i) = (f -  H13 (length(t))/ (2*pi));
    
    
end

figure, plot(Amplitudes, abs(FreqError), 'LineWidth',1.5) 
title('(a)')
xlabel('Amplitude (units) \rightarrow') 
ylabel('Absolute Frequency Error (Hz)')
axis tight 
grid on
 
% A1 = .5; 
% A2 = 67;
% A3 = 879; 

A1 = 1; 
A2 = 2;
A3 = 4; 


y2 =  A2 * sin (2*pi*f*t + phi ) - A2 + awgn(zeros(size(t))+A2,SNR,'measured');    % infecting the signal with noise proportional to the amplitude
y1 =  A1 * sin (2*pi*f*t + phi ) - A1 + awgn(zeros(size(t))+A1,SNR,'measured');    % infecting the signal with noise proportional to the amplitude
y3 =  A3 * sin (2*pi*f*t + phi ) - A3 + awgn(zeros(size(t))+A3,SNR,'measured');    % infecting the signal with noise proportional to the amplitude

figure, plot (t, A2 * sin (2*pi*f*t + phi ), t, y2, 'LineWidth', 1.5)
xlabel('Time(sec)')
ylabel('Amplitude(units)')

axis([0.4 0.5  -7.5 7.5 ])
% legend ('SNR = 5 dB')
legend ('Pure Signal', strcat('Noisy Signal, SNR =', num2str(SNR), 'dB'))
grid on

freq_offset = 18;
[H1, H2,H13]= FourthOrderANFFixedBlock(y1, T, 2*pi*(f+freq_offset), gamma, xi);



freq_offset = 18;
[H1, H2,H23]= FourthOrderANFFixedBlock(y2, T, 2*pi*(f+freq_offset), gamma, xi);


freq_offset = 18;
[H1, H2,H33]= FourthOrderANFFixedBlock(y3, T, 2*pi*(f+freq_offset), gamma, xi);


figure , plot3 (H1, H2, H13/ (2*pi), 'LineWidth', 2.3)
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
zlabel ('\theta')
colormap hsv
alpha(.2)
hold off

Indexes = find(H2 < 10^-3);     % locating point lying on the plane x_2 = 0 
 
figure, scatter (H1(Indexes ),H13 (Indexes)/ (2*pi) , '.')
% title('Poincare Section')
xlabel ('x_1')
ylabel ('\theta')
grid on

figure, plot (t, f*ones(size(t)), t, H13/ (2*pi), t, H23/ (2*pi), t, H33/ (2*pi), 'LineWidth', 2.5)
title('(b)')
legend ('True Value', strcat('Estimate, A =', num2str(A1)), strcat('Estimate, A =', num2str(A2)), strcat('Estimate, A =', num2str(A3)))
ylabel ('Frequency (Hz)')
xlabel ('time(s)')
grid on  

