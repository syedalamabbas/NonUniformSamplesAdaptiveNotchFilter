%---------------------------------Generating CLRB---------------------%
sigma=1;
A=1;

N=1000;
deltaT=1/N; %0.001;

theta=1:1:N/4;
theta=theta*2*pi;

a=0;
b=0;
c=0;
a2=0;
b2=0;
c2=0;

SNRdB = 20;

for n=1:N
    temp=(rand-1/2)/4;
  a2=a2+(n+temp)^2*deltaT^2  
  b2=b2+(n+temp)^2*deltaT^2*cos(2*((n+temp)*deltaT*theta));
  c2=c2+(n+temp)^2*deltaT^2*sin(2*((n+temp)*deltaT*theta));
   a=a+(n-1)^2*deltaT^2;
   b=b+(n-1)^2*deltaT^2*cos(2*((n-1)*deltaT*theta));
   c=c+(n-1)^2*deltaT^2*sin(2*((n-1)*deltaT*theta));
end

%f2=a2+b2;
f=a+(b.^2+c.^2).^(1/2);
f2=a2+(b2.^2+c2.^2).^(1/2);
f=f/10^(SNRdB/10);
f2=f2/10^(SNRdB/10);
figure, plot(theta/(2*pi), 2*(f2.^(-1)),'b', 'LineWidth', 1.8)
%plot(theta/(2*pi), 2*(f2.^(-1)))
xlabel('Frequency Value (Hz)')
ylabel ('Cramer-Rao Lower Bound')
axis ([0 250 0.46 0.62])
grid on

%-------------------------------------------------------------------------%