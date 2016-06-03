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
A = 2;
f =  60;
freq_offset = 18;
phi = pi/2;
f_Est = zeros(size(f));
 
SNRdBs = [0:.1:50];
FrequencyError = ones(size(SNRdBs ));

for i=1:length(SNRdBs) 
    y =  A * sin (2*pi*f*t + phi ) - A + awgn(zeros(size(t))+A,SNRdBs(i),'measured');    % infecting the signal with noise proportional to the amplitude
    
    [H1, H2,H13]= FourthOrderANFFixedBlock(y, T, 2*pi*(f+freq_offset), gamma, xi);
    
    if (i==1)
        figure, plot (t, A * sin (2*pi*f*t + phi ), t, y,t,2*xi*H2./H13,'LineWidth', 1.5)
        title ('(b)')
        xlabel('Time(sec)')
        ylabel('Amplitude(units)')
        axis([0.4 0.5  -5.5 5.5 ])
        legend ('Pure Signal, 60 Hz', strcat('Noisy @ SNR =', num2str(SNRdBs(i)), 'dB'), 'Reconstructed Signal')
        grid on
        
       


        figure, plot (t, f*ones(size(t)), t, H13/ (2*pi),  'LineWidth', 2.5)
        title ('(a)')
        legend ('True Value', strcat('Estimate, A = 2, SNR =', num2str(SNRdBs(i)), ' dB') )
        ylabel ('Frequency (Hz)')
        xlabel ('time(s)')
        grid on  
    end
    
    FrequencyError(i) = (f - H13 (length(t))/ (2*pi));
end

figure, scatter(SNRdBs, FrequencyError, 'MarkerFaceColor','m')
 title ('(a)')
xlabel('SNRdB \rightarrow')
ylabel('Frequency Error (Hz)')
axis ([0 50  -1.5 1.5 ])
grid on
 


% Find Mean value in Freq error vs SNRs
MeanFreqError = zeros(size(SNRdBs));
N = 500;

for i=1:length(SNRdBs) 
    NErrors = zeros (1,N);
    for j = 1:N
        y =  A * sin (2*pi*f*t + phi ) - A + awgn(zeros(size(t))+A,SNRdBs(i),'measured');    % infecting the signal with noise proportional to the amplitude
        
        [H1, H2,H13]= FourthOrderANFFixedBlock(y, T, 2*pi*(f+freq_offset), gamma, xi);
        
        NErrors(j) = (f - H13 (length(t))/ (2*pi));
    end
    MeanFreqError(i) = mean(NErrors);
end

figure, plot(SNRdBs, MeanFreqError, 'LineWidth',1.5)
 title ('(b)')
xlabel('SNRdB \rightarrow')
ylabel('Mean Frequency Error (Hz)')
axis tight
grid on


