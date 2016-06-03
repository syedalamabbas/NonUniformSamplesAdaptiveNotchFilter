function [varargout]=harmonic_est(x,varargin)
%
%[theta,theta_curve,b,a]=harmonic_est(x,P,Fs,theta,mu,r)
%
% 
% Implements harmonic frequency tracking algorithm as described in 
% IEEE Signal Processing Magazine (189) 11/09 by Tan and Jiang.
% Parameters are:
% 
% Input:
%     x     - (Nx1, required) Signal to be tracked
%     P     - (1x1, required) Order of how many harmonics (including fundamental)
%              to be estimated.
%     Fs    - (1x1, optional) Sampling frequency (Hz). If present output is returned
%             in Hertz, if absent output is returned in radians.
%     theta - (1x1, optional) Initial guess of the fundamental frequency capture range.
%     mu    - (1x1, optional) Step size of the LMS algorithm.
%     r     - (1x1, optinoal) Magnitude of the pole for the harmonic comb filter
%             (0<r<1).
%
% Output:
%     theta -       (1x1,required) Final estimate of the fundamental frequency.
%     theta_curve - (Nx1,optional) Tracking curve of the instantenous fundamental
%                                  frequency.
%     b-            (1xP,optional) b coefficients of the comb filter.
%     a-            (1xP,optional) a coefficients fo the comb filter.
% 
% 
% %Example
% clear all;close all;clc
% N=1200;
% Fs=8000;
% tm=[0:1/Fs:(N-1)/Fs]';
% t=[0:400:N]+1;
% snr=10^(-18/20);
% F=[1000 1075 975];
% x=[];
% true=[];
% for i=1:3
%     T=2*pi*F(i).*tm(t(i):t(i+1)-1);
%     sig=sin(T)+0.5*cos(T*2)+0.25*cos(T*3)+randn(N/3,1).*snr;
%     x=[x;sig];
%     true=[true;ones(N/3,1)*F(i)];
% end
% 
% [theta,theta_curve,b,a]=harmonic_est(x,3,Fs);
% subplot(211)
% plot(tm,theta_curve)
% hold on
% plot(tm,true,'r--','LineWidth',3)
% grid on
% xlabel('Time')
% ylabel('Fundamental Frequency Estimate')
% legend('Tracking','True')
% subplot(212)
% [H,F]=freqz(b,a,N,Fs);
% plot(F,log10(abs(H)))
% title('Final Comb Filter')
% xlabel('Frequency')
% ylabel('Magnitude')
% 
%
%
%
%%%  Written by Ikaro Silva, 2010


[N,M]=size(x);
step=500;
r= 0.85;
Fs=[];
theta=[];
mu=10^-3;


P=varargin{1}; %required
if(nargin>2 && ~isempty(varargin{2}) )
    Fs=varargin{2};
end
if(nargin>3 && ~isempty(varargin{3}) )
    theta=varargin{3};
end
if(nargin>4 && ~isempty(varargin{4}) )
    mu=varargin{4};
end
if(nargin>5 && ~isempty(varargin{5}) )
    r=varargin{5};
end

THETA=linspace(0,pi/P,step);
F=THETA;
Ntheta=length(THETA);
MSE=zeros(Ntheta,1);
MSE1=zeros(Ntheta,1);

%Step 1- Get frequency capture range
if(isempty(theta))
    for n=1:Ntheta
        C=2*cos([1:P].*THETA(n));
        B=[ones(P,1) -C' ones(P,1)];
        A=[ones(P,1) -r.*C' (r^2).*ones(P,1)];
        [b,a]=sos2tf([B A]);
        ym=filter(b,a,x);
        MSE(n)=mean(ym.^2);
        y1=filter(B(1,:),A(1,:),x);
        MSE1(n)=mean(y1.^2);
    end
    if(~isempty(Fs))
        F=F*Fs/(2*pi);
    end
    DC_MSE=mean(MSE);
    %Get initial value for theta
    [trash,theta]=max(DC_MSE-MSE1);
    theta=THETA(theta);
end

%Step 2 - Apply LMS to optimize Theta
beta=zeros(P+1,1);
ym=zeros(P+1,1);
theta_curve=zeros(N,1)+NaN;
ym_old=zeros(P+1,2);
beta_old=zeros(P+1,2);
ym_const=zeros(P+1,2);
ym_old(1,:)=[0 0];

for n=1:N

       ym=rec_step(ym_old,ym_const,x(n),theta,P+1,r);
       beta=rec_step(beta_old,ym_old(:,1),0,theta,P+1,r);
       
       ym_old=[ym ym_old(:,1)];
       beta_old=[beta beta_old(:,1)];
       
       theta_curve(n)=theta;       
       theta= theta - 2*mu*ym(end)*beta(end);
end

if(~isempty(Fs))
    theta=theta*Fs/(2*pi);
    theta_curve=theta_curve*Fs/(2*pi);
end

if(nargout > 0)
    varargout(1)={theta};
end
if(nargout> 1)
    varargout(2)={theta_curve};
end

if(nargout > 2)
    if(~isempty(Fs))
        w=theta*2*pi/Fs;
    else
        w=theta;
    end
        
    C=2*cos([1:P].*w);
    B=[ones(P,1) -C' ones(P,1)];
    A=[ones(P,1) -r.*C' (r^2).*ones(P,1)];
    [b,a]=sos2tf([B A]);
    varargout(3)={b};
    varargout(4)={a};
end



%%%%%%%%% end of main %%%%%%%%%%%%

function out = rec_step(out_old,const,init,theta,P,r)

out=zeros(P,1);
out(1)=init;

for p=2:P
    
    w=(p-1)*theta;
    out(p)= out(p-1) - 2*cos(w)*out_old(p-1,1) + ...
        2*(p-1)*sin(w)*const(p-1) + out_old(p-1,2) + ...
        2*r*cos(w)*out_old(p,1) - (r^2)*out_old(p,2) - ...
        2*r*(p-1)*sin(w)*const(p);
end



