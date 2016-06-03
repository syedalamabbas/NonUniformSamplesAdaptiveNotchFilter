
close all;

%% Whistle TEST
[wave,fs]=audioread('TRAINWHI.wav'); % read file into memory */
t=0:1/fs:(length(wave)-1)/fs; % and get sampling frequency */
N = length(t);

A = 3.5; 
wave = A * wave ;
sound(wave,fs); % see what it sounds like */
figure, plot(t,wave); 
title('Train whistle wave file.');
ylabel('Amplitude');
          xlabel('Length (in seconds)');
          grid on

whistleFreq =  330;         
myWave =  wave(1:N)' + cos(2*pi*whistleFreq*t); % +  sin(2*pi*12127.35*t)
SNRdB = 120;
 myWave =  myWave - (A+1) + awgn(zeros(size(t))+(A+1) ,SNRdB,'measured');
 
wave = myWave;
pause(2)
sound(cos(2*pi*whistleFreq*t),fs)
pause(4)
sound(wave,fs); % see what it sounds like */
pause(2)


% plot(t,wave);
figure, 
          plot(t,wave);
          title(['Superimposed train wave file with a ', num2str(whistleFreq),' Hz whistle.']);
          ylabel('Amplitude');
          xlabel('Length (in seconds)');
          grid on
          
%            axis ([ .4 .45 -6 6])
 % graph it – try zooming while its up…not much visible until you do*/
% n=length(wave)-1; 
% f=0:fs/n:fs;
% wavefft=abs(fft(wave)); % perform Fourier Transform *
% 
% figure,
%           plot(f,wavefft); % plot Fourier Transform */
%           xlabel('Frequency in Hz');
%           ylabel('Magnitude');
%           title('Frequency Spectrum');
%           grid on

          
          %-------------------------Specified parameters for the ANF----------------%
gamma = 0.00005;
xi  = .015; 
initialFreq = 2*pi* 390;
[H1, H2,H3]= NonUniform4thOrderANFFixedBlock(wave, t,initialFreq, gamma, xi);    % This is a non-uniformly sampled system
figure, plot (t,  ones(size(t))*whistleFreq,t, H3/(2*pi),'r', 'LineWidth', 2.1)
ylabel('Frequency in Hz');
xlabel('Length (in seconds)');
grid on
legend('True Value','Estimated Whistle Frequency')



% %----------------------------MATLAB Audio example check-------------------------------------%
% load mtlb
% fs = 7418;
% segmentlen = 100;
% noverlap = 90;
% NFFT = 128;
% 
% spectrogram(mtlb,segmentlen,noverlap,NFFT,fs,'yaxis')
% sound(mtlb,fs)
% dt = 1/Fs;
% I0 = round(0.1/dt);
% Iend = round(0.25/dt);
% x = mtlb(I0:Iend);
% c = cceps(x);
% t = 0:dt:length(x)*dt-dt;
% 
% 
% gamma = 0.01;
% xi  = .15; 
% initialFreq = 2*pi* 200;
% 
% 
% [H1, H2,H3]= NonUniform4thOrderANFFixedBlock(c, t,initialFreq, gamma, xi);    % This is a non-uniformly sampled system
% figure, plot (t, H3/(2*pi),'r', 'LineWidth', 2.1)
% 
% %-------------------------------------------------------------------------%
