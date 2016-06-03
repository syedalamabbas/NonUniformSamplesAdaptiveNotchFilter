function [X1,X2,X3]  = NonUniformThirdOrderANFFixedBlock (Y, modifiedt, initialFreq, gamma, xi)

%----------Y = 1x N input signal with unknown Frequency---------%
%----------T = uniform sampling interval------------------------%



N = length(Y);
%---State of Agent 1
X1 = zeros(1,N);        %Output 1
dotX1 = 0;
ddotX1 = 0;
dddotX1 = 0;

%---State of Agent 2
X2 = zeros(1,N);        %Output 2
dotX2 = 0;
ddotX2 = 0;
dddotX2 = 0;

%---State of Agent 3
X3 = zeros(1,N);        %Output 3
X3(1)  = initialFreq;
dotX3 = 0;
ddotX3 = 0;
dddotX3 = 0;

dotY = 0;
ddotY = 0;

a1 = 0;
a2 = 0;
a3 = 0;

%-------------------------------Hsu Globally convergent Estimator-------%
for i =2:N
    
    T = double(modifiedt(i)- modifiedt(i-1));
    
%     dotY = [Y(i)-Y(i-1)+T^2*X3(i-1)^2*Y(i-1)/2 - T^4*X3(i-1)^4*Y(i-1)/24]/(T - X3(i-1)^2*T^3/6);      %  Using Backward Differentiation Formula
     dotY = X1(i-1)*X3(i-1)*(-2*xi);
    ddotY = -X3(i-1)^2*Y(i-1);
    
    dotX1 = X2(i-1);
    dotX2 = -2*xi*X3(i-1)*X2(i-1) + X3(i-1)^2*[Y(i-1)-X1(i-1)] ;
    dotX3 = -gamma*[X3(i-1)^2*Y(i-1)- 2*xi*X3(i-1)*X2(i-1)]*X1(i-1);
    
    ddotX1 = dotX2;
    ddotX2 = -2*xi*[X3(i-1)*dotX2+dotX3*X2(i-1)]+ X3(i-1)^2* [dotY - dotX1] +2*X3(i-1)*dotX3*[Y(i-1)-X1(i-1)];
    ddotX3 = -gamma*[X3(i-1)^2*Y(i-1)- 2*xi*X3(i-1)*X2(i-1)]*dotX1- gamma*X1(i-1)*[ X3(i-1)^2*dotY + 2*X3(i-1)*dotX3*Y(i-1)-2*xi*X3(i-1)*dotX2-2*xi*dotX3*X2(i-1)];
    
    dddotX1 = ddotX2;
    dddotX2 = -2*xi*[X3(i-1)*ddotX2+ddotX3*X2(i-1)+2*dotX2*dotX3]+ X3(i-1)^2* [ddotY - ddotX1] + 4* X3(i-1)*dotX3 * [dotY - dotX1] + (2*X3(i-1)*ddotX3 + 2*dotX3^2)*[Y(i-1)-X1(i-1)];
    a1 = -gamma*[X3(i-1)^2*Y(i-1)- 2*xi*X3(i-1)*X2(i-1)]*ddotX1;
    a2 = - 2*gamma*dotX1*[ X3(i-1)^2*dotY + 2*X3(i-1)*dotX3*Y(i-1)-2*xi*(X3(i-1)*dotX2+dotX3*X2(i-1))];
    a3 =  - gamma*X1(i-1)*( X3(i-1)^2*ddotY + 4*X3(i-1)*dotX3*Y(i-1) +2*dotX3^2*Y(i-1)+2*X3(i-1)*ddotX3*Y(i-1) -2*xi*( X3(i-1)*ddotX2 + ddotX3*X2(i-1) + 2*dotX2*dotX3 ));
    dddotX3 = a1 + a2 + a3;
    
    %-----------Agent 1- Solving 1st Linear Differential Equation------%
        X1(i) = X1(i-1) + T*dotX1 + T^2/2* ddotX1 + T^3/6*dddotX1;
    %-------------------------------------------------------------------%
    
    %-----------Agent 2- Solving 2nd Linear Differential Equation------%
        X2(i) = X2(i-1) + T*dotX2 + T^2/2* ddotX2+ T^3/6*dddotX2;
    %------------------------------------------------------------------%
    
    %-----------Agent 3- Solving 3rd Linear Differential Equation------%
        X3(i) = X3(i-1) + T*dotX3 + T^2/2* ddotX3+ T^3/6*dddotX3;
    %------------------------------------------------------------------%
end
%----------------------------------------------------------------------%
end
