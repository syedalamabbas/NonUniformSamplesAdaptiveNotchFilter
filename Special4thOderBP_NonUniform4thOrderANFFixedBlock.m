function [ X5,X6,Theta]  = Special4thOderBP_NonUniform4thOrderANFFixedBlock (Y, modifiedt, initialFreq, gamma, xi, bandPassFreq, xi_c)

%----------Y = 1x N input signal with unknown Frequency---------%
%----------T = uniform sampling interval------------------------%

N = length(Y);
%---State of Agent 1

X1 = zeros(1,N);        %Output 1
DX1 = 0;
D2X1 = 0;
D3X1 = 0;
D4X1 = 0;


X2 = zeros(1,N);        %Output 1
DX2 = 0;
D2X2 = 0;
D3X2 = 0;
D4X2 = 0;


X3 = zeros(1,N);        %Output 1
DX3 = 0;
D2X3 = 0;
D3X3 = 0;
D4X3 = 0;

%---State of Agent 2
X4 = zeros(1,N);        %Output 2
DX4 = 0 ;
D2X4 = 0;
D3X4 = 0;
D4X4 = 0;


X5 =  zeros(1,N);
DX5 = 0;
D2X5 = 0;
D3X5 = 0;
D4X5 = 0;

X6 =  zeros(1,N);
DX6 = 0;
D2X6 = 0;
D3X6 = 0;
D4X6 = 0;

%---State of Agent 5
Theta = zeros(1,N);        %Output 5
Theta(1)  = initialFreq;
DTheta = 0;
D2Theta = 0;
D3Theta = 0;
D4Theta = 0;


DY = 0;
D2Y = 0;
D3Y = 0;
D4Y = 0;

% bandPassFreq = initialFreq-10;
% xi_c = xi+.2;

% Y =  * Y;

%-------------------------------Hsu Globally convergent Estimator-------%
for i =2:N
    
    T = double(modifiedt(i)- modifiedt(i-1));
    
    %-------------------------Signal Derivatives-------------------------%
%      bandPassFreq = Theta(i-1);
    DY = X5(i-1)*Theta(i-1)*(-2*xi);
    D2Y = -Theta(i-1)^2*Y(i-1);
    D3Y = -Theta(i-1)^2*DY;
    D4Y = Theta(i-1)^4*Y(i-1);
    
    
    %----------------------------1st order--------------------------------%
    DX1 = X2(i-1);
    DX2 = 2*bandPassFreq*(xi_c)*DY - bandPassFreq^2* X1(i-1) - 2*bandPassFreq*(xi_c )* X2(i-1);
    DX3 = X4(i-1);
    DX4 = 2*bandPassFreq*(xi_c)*DX1 - bandPassFreq^2* X3(i-1) - 2*bandPassFreq*(xi_c )* X4(i-1);
    DX5 = X6(i-1);
    DX6 = Theta(i-1)^2*X3(i-1)- Theta(i-1)^2*X5(i-1) - 2*xi*Theta(i-1)*X6(i-1);
    DTheta = -gamma*X5(i-1)*(DX6+Theta(i-1)^2*X5(i-1));
    %---------------------------------------------------------------------%
    
    %-----------------------2nd Order-------------------------------------%
    D2X1 = DX2;
    D2X2 = 2*bandPassFreq*(xi_c)*D2Y - bandPassFreq^2* DX1 - 2*bandPassFreq*(xi_c )* DX2;
    D2X3 = DX4;
    D2X4 = 2*bandPassFreq*(xi_c)*D2X1 - bandPassFreq^2* DX3 - 2*bandPassFreq*(xi_c)* DX4;
    D2X5 = DX6;
    D2X6 = Theta(i-1)^2*DX3- Theta(i-1)^2*DX5 - 2*xi*Theta(i-1)*DX6 - 2*xi*X6(i-1)*DTheta + 2*Theta(i-1)*X3(i-1)*DTheta - 2*Theta(i-1)*X5(i-1)*DTheta ;
    D2Theta = -gamma*X5(i-1)*(Theta(i-1)^2*DX5 + D2X6 + 2*Theta(i-1)*DTheta*X5(i-1)) -gamma*DX5*(DX6+Theta(i-1)^2*X5(i-1));
    %---------------------------------------------------------------------%
    
    %-----------------------3rd Order-------------------------------------%
    D3X1 = D2X2;
    D3X2 = 2*bandPassFreq*(xi_c)*D3Y - bandPassFreq^2* D2X1 - 2*bandPassFreq*(xi_c )* D2X2;
    D3X3 = D2X4;
    D3X4 = 2*bandPassFreq*(xi_c)*D3X1 - bandPassFreq^2* D2X3 - 2*bandPassFreq*(xi_c)* D2X4;
    D3X5 = D2X6;
    D3X6 = 2*X3(i-1)*DTheta^2 - 2* X5(i-1)* DTheta^2 + Theta(i-1)^2*D2X3 - Theta(i-1)^2*D2X5 - 4*xi*DTheta*DX6 + 2*Theta(i-1)* X3(i-1)*D2Theta - 2*Theta(i-1)* X5(i-1)*D2Theta + 4*Theta(i-1)*DTheta*DX3 - 4*Theta(i-1)*DTheta*DX5 - 2*xi*Theta(i-1)*D2X6 - 2*xi*X6(i-1)*D2Theta;
    D3Theta = -gamma*X5(i-1)*( 2*X5(i-1)*DTheta^2 + Theta(i-1)^2*D2X5 + D3X6 + 2*Theta(i-1)*D2Theta*X5(i-1) + 4*Theta(i-1)*DTheta*DX5) -gamma*D2X5*(DX6+Theta(i-1)^2*X5(i-1)) - 2*gamma*DX5*(Theta(i-1)^2*DX5 + D2X6+ 2*Theta(i-1)*X5(i-1)*DTheta);
    %---------------------------------------------------------------------%
    
    
    %-----------------------4th Order-------------------------------------%
    D4X1 = D3X2;
    D4X2 = 2*bandPassFreq*(xi_c)*D4Y - bandPassFreq^2* D3X1 - 2*bandPassFreq*(xi_c )* D3X2;
    D4X3 = D3X4;
    D4X4 = 2*bandPassFreq*(xi_c)*D4X1 - bandPassFreq^2* D3X3 - 2*bandPassFreq*(xi_c )* D3X4;
    D4X5 = D3X6;
    D4X6 =  6*DTheta^2*DX3 - 6*DTheta^2*DX5 + Theta(i-1)^2*D3X3 - Theta(i-1)^2*D3X5 + 6*Theta(i-1)*DTheta*D2X3 + 6*Theta(i-1)*DX3* D2Theta + 6*X3(i-1)*DTheta*D2Theta - 6*Theta(i-1)*DTheta*D2X5 - 6*Theta(i-1)*DX5*D2Theta - 6*X5(i-1)*DTheta*D2Theta - 2*xi*Theta(i-1)*D3X6 - 2*xi*X6(i-1)*D3Theta -6*xi*DTheta*D2X6 - 6*xi*DX6*D2Theta +2*Theta(i-1)*X3(i-1)*D3Theta - 2*Theta(i-1)*X5(i-1)*D3Theta;
    D4Theta = -gamma*D3X5*(DX6+Theta(i-1)^2*X5(i-1))-3*gamma*D2X5*( Theta(i-1)^2*DX5 + D2X6 + 2*Theta(i-1)*X5(i-1)*DTheta) - 3*gamma*DX5*( 2*X5(i-1)*DTheta^2 + Theta(i-1)^2*D2X5 + D3X6 + 2*Theta(i-1)*X5(i-1)*D2Theta + 4*Theta(i-1)*DTheta*DX5 )     -gamma*X5(i-1)*( 6*DX5*DTheta^2 + Theta(i-1)^2*D3X5 + 6*Theta(i-1)*DTheta*D2X5 + 6*Theta(i-1)*D2Theta*DX5 + 6*X5(i-1)*DTheta*D2Theta  + 2*Theta(i-1)*D3Theta*X5(i-1) +  D4X6 );
    %---------------------------------------------------------------------%
    
    
    %-----------Agents- Solving n Linear Differential Equations------%
    X1(i) = X1(i-1) + T*DX1 + T^2/2*D2X1 + T^3/6*D3X1+ T^4/24*D4X1; 
    X2(i) = X2(i-1) + T*DX2 + T^2/2*D2X2 + T^3/6*D3X2+ T^4/24*D4X2;
    X3(i) = X3(i-1) + T*DX3 + T^2/2*D2X3 + T^3/6*D3X3+ T^4/24*D4X3; 
    X4(i) = X4(i-1) + T*DX4 + T^2/2*D2X4 + T^3/6*D3X4+ T^4/24*D4X4;
    X5(i) = X5(i-1) + T*DX5 + T^2/2*D2X5 + T^3/6*D3X5+ T^4/24*D4X5;
    X6(i) = X6(i-1) + T*DX6 + T^2/2*D2X6 + T^3/6*D3X6+ T^4/24*D4X6;
    Theta(i) = Theta(i-1) + T*DTheta + T^2/2*D2Theta + T^3/6*D3Theta+ T^4/24*D4Theta;
    
    %------------------------------------------------------------------%
end

end
