clear 
clc
close all

theta = 2*pi* 500; xi =1.4;
s=tf('s');

%% Band Pass Filter
H = (2*xi*theta*s)/(s^2+ 2*xi*theta*s+ theta^2);
figure, 
h = bodeplot(H)
p=getoptions(h);    % Create a plot options handle p.
p.FreqUnits = 'Hz'; % Modify frequency units.
setoptions(h,p);    % Apply plot options to the Bode plot and render.
grid on

 
%% Band Reject Filter
H = (s^2+ theta^2)/(s^2+ 2*xi*theta*s+ theta^2);
figure, 
h1 = bodeplot(H)
p=getoptions(h1);    % Create a plot options handle p.
p.FreqUnits = 'Hz'; % Modify frequency units.
setoptions(h1,p);    % Apply plot options to the Bode plot and 
                    % render.
grid on


%% Low Pass Filter
H = (theta^2)/(s^2+ 2*xi*theta*s+ theta^2);
figure, 
h2 = bodeplot(H)
p=getoptions(h2);    % Create a plot options handle p.
p.FreqUnits = 'Hz'; % Modify frequency units.
setoptions(h2,p);    % Apply plot options to the Bode plot and 
                    % render.
grid on



%% High Pass Filter
H = (s^2)/(s^2+ 2*xi*theta*s+ theta^2);
figure, 
h3 = bodeplot(H)
p=getoptions(h3);    % Create a plot options handle p.
p.FreqUnits = 'Hz'; % Modify frequency units.
setoptions(h3,p);    % Apply plot options to the Bode plot and 
                    % render.
grid on