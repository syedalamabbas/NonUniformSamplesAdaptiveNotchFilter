function [X1,X2,M]  = NonUniform4thOrderANFFixedBlock (Y, modifiedt, initialFreq, gamma, xi)

%----------Y = 1x N input signal with unknown Frequency---------%
%----------T = uniform sampling interval------------------------%

N = length(Y);
%---State of Agent 1
X1 = zeros(1,N);        %Output 1
DX1 = 0;
D2X1 = 0;
D3X1 = 0;
D4X1 = 0;

%---State of Agent 2
X2 = zeros(1,N);        %Output 2
DX2 = 0;
D2X2 = 0;
D3X2 = 0;
D4X2 = 0;

%---State of Agent 3
M = zeros(1,N);        %Output 3
M(1)  = initialFreq;
DM = 0;
D2M = 0;
D3M = 0;
D4M = 0;


X3 = 0; 
DX3 = 0;
D2X3 = 0;
D3X3 = 0;
D4X3 = 0;

X4 = 0; 
DX4 = 0;
D2X4 = 0;
D3X4 = 0;
D4X4 = 0;

X5 = 0; 
DX5 = 0;
D2X5 = 0;
D3X5 = 0;
D4X5 = 0;

X6 = 0; 
DX6 = 0;
D2X6 = 0;
D3X6 = 0;
D4X6 = 0;


DY = 0;
D2Y = 0;
D3Y = 0;
D4Y = 0;
D5Y = 0;

%-------------------------------Hsu Globally convergent Estimator-------%
for i =2:N
    %      DY = [Y(i)-Y(i-1)+T^2*M(i-1)^2*Y(i-1)/2 + T^4*M(i-1)^4*Y(i-1)/24]/(T - M(i-1)^2*T^3/6);      %  Using Backward Differentiation Formula
    
    T = double(modifiedt(i)- modifiedt(i-1));
%     T = 0.001;

    DY = X1(i-1)*M(i-1)*(-2*xi);
%     DY = (Y(i)-Y(i-1)*(1 - T^2*M(i-1)^2/2 + T^4*M(i-1)^4/24) )/(T-M(i-1)^2*T^3/6 + M(i-1)^4*T^5/120 );
    
    X3 = M(i-1)^2;
    X4 = X3*Y(i-1);
    X5 = 2*xi*M(i-1)*X2(i-1);
    X6 = X1(i-1)*X3;
    
    D2Y= -X3*Y(i-1);
    D3Y = -X3*DY;
    D4Y = X3^2*Y(i-1);
    D5Y = X3^2*DY;
    
    
    %----------------------------1st order--------------------------------%
    DX1 = X2(i-1);
    DX2 = X4 - X5 - X6;
    DM = -gamma*(X4 - X5)*X1(i-1);
    DX3 = 2*M(i-1)*DM;
    DX4 = X3*DY+Y(i-1)*DX3;
    DX5 = 2*xi*(M(i-1)*DX2+ X2(i-1)*DM);
    DX6 = X1(i-1)*DX3 + X3*DX1; 
    %---------------------------------------------------------------------%
    
    %-----------------------2nd Order-------------------------------------%
    D2X1 = DX2;
    D2X2 = DX4 - DX5 - DX6;
    D2M = -gamma*( (X4 - X5)*DX1 + (DX4 - DX5)*X1(i-1));
    D2X3 = 2*(M(i-1)*D2M+ DM^2);
    D2X4 = -X3^2*Y(i-1) + 2*DX3*DY+Y(i-1)*D2X3;
    D2X5 = 2*xi*(M(i-1)*D2X2+ 2*DM*DX2 + X2(i-1)*D2M);
    D2X6 = X1(i-1)*D2X3 + 2*DX3*DX1 +X3*D2X1; 
    %---------------------------------------------------------------------%
    
    %-----------------------3rd Order-------------------------------------%
    D3X1 = D2X2;
    D3X2 = D2X4 - D2X5 - D2X6;
    D3M = -gamma*( (X4 - X5)*D2X1 + 2*(DX4 - DX5)*DX1 + (D2X4 - D2X5)*X1(i-1));
    D3X3 = 2*(M(i-1)*D3M+ 3*DM*D2M);
    D3X4 = -X3*(M(i-1)^2*DY+2*M(i-1)*DM*Y(i-1)) - 3*X4*DX3+ 3*D2X3*DY+Y(i-1)*D3X3;
    D3X5 = 2*xi*(M(i-1)*D3X2+ 3*D2M*DX2 +3*DM*D2X2 + X2(i-1)*D3M);
    D3X6 = X1(i-1)*D3X3 + 3*D2X3*DX1+ 3*DX3*D2X1 +X3*D3X1; 
    %---------------------------------------------------------------------%
    
    
    %-----------------------4th Order-------------------------------------%
    D4X1 = D3X2;
    D4X2 = D3X4 -D3X5 - D3X6;
    D4M = -gamma*((X4 - X5)*D3X1 + 3*(D2X4 - D2X5)*DX1 + 3*(DX4 - DX5)*D2X1 + (D3X4 - D3X5)*X1(i-1) );
    D4X3 = 2*(M(i-1)*D3M+ 3*DM*D2M);
    D4X4 = -X3*(M(i-1)^2*DY+2*M(i-1)*DM*Y(i-1)) - 3*X4*DX3+ 3*D2X3*DY+Y(i-1)*D3X3;
    D4X5 = 2*xi*(M(i-1)*D3X2+ 3*D2M*DX2 +3*DM*D2X2 + X2(i-1)*D3M);
    D4X6 = X1(i-1)*D3X3 + 3*D2X3*DX1+ 3*DX3*D2X1 +X3*D3X1; 
    %---------------------------------------------------------------------%
    
    
    %------------------5th order---------------------------%
    D5X1=D4X2;
    D5X2= 12*DM^2*D2Y - 12*DM^2*D2X1 - M(i-1)^2*D4X1 + (M(i-1))^2*D4Y - 6*X1(i-1)*D2M^2 + 6*Y(i-1)*D2M^2 - 24*DM*DX1*D2M + 24*DM*DY*D2M - 8*xi*DM*D3X2 - 8*xi*DX2*D3M - 12*xi*D2M*D2X2 - 2*M(i-1)*X1(i-1)*D4M + 2*M(i-1)*Y(i-1)*D4M - 8*M(i-1)*DM*D3X1 - 8*M(i-1)*DX1*D3M - 12*M(i-1)*D2M*D2X1 - 8*X1(i-1)*DM*D3M + 8*M(i-1)*DM*D3Y + 8*M(i-1)*DY*D3M + 12*M(i-1)*D2M*D2Y + 8*Y(i-1)*DM*D3M - 2*xi*M(i-1)*D4X2 - 2*xi*X2(i-1)*D4M;
    D5M=- 6*gamma*D2X1*(2*Y(i-1)*DM^2 + M(i-1)^2*D2Y - 4*xi*DM*DX2 + 2*M(i-1)*Y(i-1)*D2M + 4*M(i-1)*DM*DY - 2*xi*M(i-1)*D2X2 - 2*xi*X2(i-1)*D2M) - 4*gamma*D3X1*(M(i-1)^2*DY - 2*xi*M(i-1)*DX2 - 2*xi*X2(i-1)*DM + 2*M(i-1)*Y(i-1)*DM) - gamma*X1(i-1)*(12*DM^2*D2Y + M(i-1)^2*D4Y + 6*Y(i-1)*D2M^2 + 24*DM*DY*D2M - 8*xi*DM*D3X2 - 8*xi*DX2*D3M - 12*xi*D2M*D2X2 + 2*M(i-1)*Y(i-1)*D4M + 8*M(i-1)*DM*D3Y + 8*M(i-1)*DY*D3M + 12*M(i-1)*D2M*D2Y + 8*Y(i-1)*DM*D3M - 2*xi*M(i-1)*D4X2 - 2*xi*X2(i-1)*D4M) - gamma*(M(i-1)^2*Y(i-1) - 2*xi*M(i-1)*X2(i-1))*D4X1 - 4*gamma*DX1*(6*DM^2*DY + M(i-1)^2*D3Y + 6*M(i-1)*DM*D2Y + 6*M(i-1)*DY*D2M + 6*Y(i-1)*DM*D2M - 2*xi*M(i-1)*D3X2 - 2*xi*X2(i-1)*D3M - 6*xi*DM*D2X2 - 6*xi*DX2*D2M + 2*M(i-1)*Y(i-1)*D3M);

    
%     D5X1 = D4X2;
%     D5X2 = -12*DM ^2*X4 - 12*DM^2*D2X1 - X3*D4X1 + X3^3*Y(i-1) - 6*X1(i-1)*D2M^2 + 6*Y(i-1)*D2M^2 - 24*DM*DX1*D2M + 24*DM*DY*D2M - 8*xi*DM*D3X2 - 8*xi*DX2*D3M - 12*xi*D2M*D2X2 - 2*M(i-1)*X1(i-1)*D4M + 2*M(i-1)*Y(i-1)*D4M - 8*M(i-1)*DM*D3X1 - 8*M(i-1)*DX1*D3M - 12*M(i-1)*D2M*D2X1 - 8*X1(i-1)*DM*D3M - 8*M(i-1)*DM*X3*DY + 8*M(i-1)*DY*D3M - 12*M(i-1)*D2M*X3*Y(i-1) + 8*Y(i-1)*DM*D3M - 2*xi*M(i-1)*D4X2 - 2*xi*X2(i-1)*D4M;
%     D5M  = - 6*gamma*D2X1*(2*Y(i-1)*DM^2 - X3^2*Y(i-1) - 4*xi*DM*DX2 + 2*M(i-1)*Y(i-1)*D2M + 4*M(i-1)*DM*DY - 2*xi*M(i-1)*D2X2 - 2*xi*X2(i-1)*D2M - 4*gamma*D3X1*((M(i-1))^2*DY - 2*xi*M(i-1)*DX2 - 2*xi*X2(i-1)*DM + 2*M(i-1)*Y(i-1)*DM) - gamma*X1(i-1)*(-12*DM^2*X3*Y(i-1) + X3^2*Y(i-1) + 6*Y(i-1)*D2M^2 + 24*DM*DY*D2M - 8*xi*DM*D3X2 - 8*xi*DX2*D3M - 12*xi*D2M*D2X2 + 2*M(i-1)*Y(i-1)*D4M - 8*M(i-1)*DM*X3*DY + 8*M(i-1)*DY*D3M - 12*M(i-1)*D2M*X3*Y(i-1) + 8*Y(i-1)*DM*D3M - 2*xi*M(i-1)*D4X2 - 2*xi*X2(i-1)*D4M) - gamma*((M(i-1))^2*Y(i-1) - 2*xi*M(i-1)*X2(i-1))*D4X1 - 4*gamma*DX1*(6*DM^2*DY - (M(i-1))^4*DY - 6*M(i-1)^3*DM*Y(i-1) + 6*M(i-1)*DY*D2M + 6*Y(i-1)*DM*D2M - 2*xi*M(i-1)*D3X2 - 2*xi*X2(i-1)*D3M - 6*xi*DM*D2X2 - 6*xi*DX2*D2M + 2*M(i-1)*Y(i-1)*D3M ));
%     
    %-------------------------------------------------------------%
    
    
    %-----------Agents- Solving n Linear Differential Equations------%
    X1(i) = X1(i-1) + T*DX1 + T^2/2*D2X1 + T^3/6*D3X1+ T^4/24*D4X1+ T^5/120*D5X1;
    X2(i) = X2(i-1) + T*DX2 + T^2/2*D2X2 + T^3/6*D3X2+ T^4/24*D4X2 + T^5/120*D5X2;
    M(i) = M(i-1) + T*DM + T^2/2*D2M + T^3/6*D3M+ T^4/24*D4M+ + T^5/120*D5M;

    %------------------------------------------------------------------%
end
%----------------------------------------------------------------------%
end
