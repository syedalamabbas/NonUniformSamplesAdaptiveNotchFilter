clear;
clc;
close all; 

%----------------------------Few Initializations-----------------------%
T = 0.001;                      % Sampling interval
a = 0;                          % t minimum
b = 1;                          % t maximum
t= [a:T:b];                     % total time interval
%----------------------------------------------------------------------%

%--------------------------ANF parameters-----------------------------%
% gamma = .9;             % This was the setting for the previous figure
% xi = 1.4;

gamma = .01;
xi = .15;
%---------------------------------------------------------------------%


%----------------------Defining the signal---------------------------%

Non_Uni_Samples = length(t);
modifiedt = zeros (1, Non_Uni_Samples);

min = -5;                   % Min value
max = 5;                    % Max value  
r = (max-min).*rand(1,Non_Uni_Samples) + min;  % Randomly distributed Values

for i = 2:Non_Uni_Samples                               % Jittered model with additive non-uniform sampling
   modifiedt(i) = modifiedt(i-1) + T+ r(i)*T*10^-1;
end

modifiedt = t;




% Non_Uni_Samples = length(t);
% modifiedt = zeros (1, Non_Uni_Samples);
% for i = 2:Non_Uni_Samples                               % Jittered model with additive non-uniform sampling
%    modifiedt(i) = modifiedt(i-1) + T+ rand()*10^-4;
% end


Y = zeros(size(t));
changingfreq = zeros(size(t));
A = 1;
for i =1:Non_Uni_Samples
    if (t(i) <= .33)
        changingfreq(i)  = 72;
    end
    if (t(i) > .33 && t(i) <= .66)
        changingfreq(i)  = 60 + 0*(t(i)-2);
    end
    if (t(i) > .66)
        changingfreq(i)  = 80 + 0*(modifiedt(i)-4);
    end
    Y(i) = A * sin (2*pi*changingfreq(i)*modifiedt(i)+ pi/3);         % Excitation Signal
end
SNRdB = 120;
Y = Y - A + awgn(zeros(size(t))+A,SNRdB,'measured');
% figure,plot(t, Y)
% xlabel ('Time (seconds)')
% ylabel ('Amplitude (units)')
% grid on
% axis tight
% print -depsc2 HopSignal

%----------------------------------------------------------------------%



%----------------------------Processing Solutions---------------------%
initialFreq = 2*pi*55;       % initialization error

%---------------------Hsu Globally convergent Estimator-----------------%
% [ X1,X2,X3] = Special4thOderBP_NonUniform4thOrderANFFixedBlock (Y, modifiedt, initialFreq, gamma, xi, initialFreq,xi+2);

  [X1,X2,X3] = NonUniform4thOrderANFFixedBlock(Y, modifiedt, initialFreq, gamma, xi);
%  [X1,X2,X3] = SpecialComplementary_NonUniform4thOrderANFFixedBlock(Y, modifiedt, initialFreq, gamma, xi+.2);
%----------------------------------------------------------------------%

%[MX1, MX2, MX3] = NonUniform5thOrderANFFixedBlock(Y, t, initialFreq, gamma, xi);

%--------------------------Mojiri 1 ANF------------------------------------%

% [MX1, MX2, MX3] = ThirdOrderMojiriANFBlock( Y, T, initialFreq, gamma, xi);
%----------------------------------------------------------------------%
%-----------------------------MATLAB Ode45 --------------------------------%
initialvalues=[0 0 initialFreq];
[t1,x1]=ode45('FunctionANF',t,initialvalues);
%-------------------------------------------------------------------------%
%------------------------------------------------------------------------%


%-------------------------Plotting Solutions------------------------------%

linewidth = 2.3;

% figure, plot3(x1(:,1),x1(:,2),x1(:,3)/(2*pi), 'LineWidth', linewidth)
% axis([-.4 .4 -200 200 30 80])
% % title('3D state space plot from MATLAB ode45 solver')
% % grid on
% xlabel('x1')
% ylabel('x2')
% zlabel('\theta (Hz)')
% grid on
% print -depsc2 MATLAB3DState

figure, plot3(X1,X2,X3/(2*pi),'r','LineWidth', 1)
% axis([-.4 .4 -200 200 30 80])
% title('3D state space plot from Agent based ODE solver Hsu ANF')
% grid on
xlabel('x_1')
ylabel('x_2')
zlabel('\theta (Hz)')
grid on
% print -depsc2 HsuANF3D

Indexes = find(X2 < 10^-3);     % locating point lying on the plane x_2 = 0 

figure, scatter (X1(Indexes ),X3 (Indexes)/ (2*pi) , '.')
% title('Poincare Section')
xlabel ('x_1')
ylabel ('\theta (Hz)')
grid on

% figure, plot3(MX1,MX2,MX3/(2*pi),'-.g','LineWidth', linewidth)
% axis([-.5 1.5 -150 800 30 80])
% % title('3D state space plot from Agent based ODE solver Mojiri ANF')
% % grid on
% xlabel('x1')
% ylabel('x2')
% zlabel('\theta (Hz)')
% grid on
% print -depsc2 MojiriANF3D

 figure, plot(t, changingfreq,'g', t, x1(:,3)'/(2*pi), ':r' , t,X3/(2*pi), '-.b',  'LineWidth',2.3)
%  t,MX3/(2*pi), '-.g',
% figure, plot(t, changingfreq,':k','LineWidth', linewidth)
hold on
% plot(t,x1(:,3)'/(2*pi),'b','LineWidth', linewidth-1) 
% plot(t,X3/(2*pi),'-.g','LineWidth', linewidth) 
% plot(t,MX3/(2*pi),'-.r' )
hold off
% grid on
axis([a b 55 90])
Ax = legend('True value','MATLAB ode45 solution','Proposed solution','Location','northeast' )
% Ax = legend('true value','MATLAB ode45(Hsu ANF)','Proposed Discrete Filter','Location','northeast' )
Ax.FontSize = 14;
% leg = findobj(Ax,'type','text')
% set(leg,'fontsize', 14)
ylabel('Frequency (Hz)')
xlabel('Time (seconds)')
grid on
% print -depsc2 SingleFreqComparison

% figure, semilogy( t, abs(changingfreq- X3/(2*pi)), '-.b',t, abs(changingfreq- MX3/(2*pi)), '-.g')
 figure, semilogy( t,abs(changingfreq- x1(:,3)'/(2*pi)), '-.r', t, abs(changingfreq- X3/(2*pi)), '-.b', 'LineWidth', linewidth)
% figure, semilogy( t,abs(changingfreq- x1(:,3)'/(2*pi)), '-.r', t, abs(changingfreq- X3/(2*pi)), '-.b', 'LineWidth', linewidth-1)
axis([a b 10^(-8) 10^4])
 Ax = legend('MATLAB ode45 solution', 'Proposed solution (4^{th} order)', 'Location','southeast')
% Ax = legend('MATLAB ode45(Hsu ANF)', 'Proposed Discrete Filter', 'Location','southeast')
leg = findobj(Ax,'type','text')
set(leg,'fontsize', 14)
% legend('Hsu ANF','Mojiri ANF' )
ylabel('Absolute Error in Frequency Estimates (Hz)')
xlabel('Time (seconds)')
grid on
% print -depsc2  SingleFreqErrorComparison

%------------------------------------------------------------------------%


% iter = length(t);
% insamples = Y;
% xst = zeros(2,1); % state of adaptive filter.
% temp = zeros(2,1); % intermediate signals.
% pihalf = 0.5*pi; % pi/2.
% theta = 0; % initial value for notch frequency parameter.
% sth = sin(theta);
% cth = cos(theta);
% bw = 0.0020*pi; % bandwidth parameter for notch filter.
% sth2 = sin(bw);
% cth2 = cos(bw);
% mu = 0.0005; % adaptive filter step size.
% %
% % Run adaptive lattice notch filter:
% %
% figure,
% plot(t, changingfreq,':k','LineWidth', linewidth+1)
% hold on
% freqold = 70*T; %initial notch frequency.
% for kk=1:iter
%  insig = insamples(kk);
%  temp = [cth2 -sth2;sth2 cth2]*[insig; xst(2)];
%  error = mu*(insig + temp(2)); %notch filter output times step size.
%  theta = theta - error*xst(1); %coefficient update.
%  freqnew = 0.5*acos(cos(theta+pihalf))/pi; %instantaneous freq. estimate.
%  if kk>1
%    plot([kk-1 kk]*T, [freqold freqnew]/T,'LineWidth', linewidth+1 );
%  end
%  freqold = freqnew;
%  sth = sin(theta);
%  cth = cos(theta);
%  xst = [cth -sth;sth cth] * [temp(1);xst(1)];
% end
% % title('Frequency Estimates of Lattice Adaptive Notch Filter');
% Ax = legend ('true value', 'Lattice Adaptive Notch Filter(Regalia)')
% leg = findobj(Ax,'type','text')
% set(leg,'fontsize', 18)
% axis([a b 50 90])
% hold off
% grid on
% print -depsc2 RegaliaFilterComp
