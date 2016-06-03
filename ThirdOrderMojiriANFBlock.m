function [ MX1,MX2,MX3 ] = ThirdOrderMojiriANFBlock( Y, T, initialFreq, gamma, xi )
%THIRDORDERMOJIRIANFBLOCK  This function simply evaluates the solution of
%the dynamical system proposed by Mojiri (2007) using a 3rd order approximation

%------------------Mojiri 1 ANF----------------------------------------%
MX1 = zeros(size(Y));        %Output 1
MX2 = zeros(size(Y));        %Output 2
MX3 = zeros(size(Y));        %Output 3
MX3(1) = initialFreq;        % Given input in radians
for i =2:length(Y)
    
    dotY = [Y(i)-Y(i-1)+T^2*MX3(i-1)^2*Y(i-1)/2]/T;
    dotdotY = -MX3(i-1)^2*Y(i-1);
    
    dotX1 = MX2(i-1);
    dotX2 = -2*xi*MX3(i-1)*MX2(i-1) + MX3(i-1)^2*[2*xi*Y(i-1)-MX1(i-1)] ;
    dotX3 = -gamma*[MX3(i-1)^2*Y(i-1)- MX3(i-1)*MX2(i-1)]*MX1(i-1);
    
    dotdotX1 = dotX2;
    dotdotX2 = -2*xi*[MX3(i-1)*dotX2+dotX3*MX2(i-1)]+ MX3(i-1)^2* [2*xi*dotY - dotX1] +2*MX3(i-1)*dotX3*[2*xi*Y(i-1)-MX1(i-1)];
    dotdotX3 = -gamma*[MX3(i-1)^2*Y(i-1)- MX3(i-1)*MX2(i-1)]*dotX1- gamma*MX1(i-1)*[ MX3(i-1)^2*dotY + 2*MX3(i-1)*dotX3*Y(i-1)-MX3(i-1)*dotX2-dotX3*MX2(i-1)];
    
    dotdotdotX1 = dotdotX2;
    dotdotdotX2 = -2*xi*[MX3(i-1)*dotdotX2+dotdotX3*MX2(i-1)+2*dotX2*dotX3]+ MX3(i-1)^2* [dotdotY - dotdotX1] + 4* MX3(i-1)*dotX3 * [2*xi*dotY - dotX1] + (2*MX3(i-1)*dotdotX3 + 2*dotX3^2)*[2*xi*Y(i-1)-MX1(i-1)];
    a1 = -gamma*[MX3(i-1)^2*Y(i-1)- MX3(i-1)*MX2(i-1)]*dotdotX1;
    a2 = - 2*gamma*dotX1*[ MX3(i-1)^2*dotY + 2*MX3(i-1)*dotX3*Y(i-1)-(MX3(i-1)*dotX2+dotX3*MX2(i-1))];
    a3 =  - gamma*MX1(i-1)*( MX3(i-1)^2*dotdotY + 4*MX3(i-1)*dotX3*Y(i-1) +2*dotX3^2*Y(i-1)+2*MX3(i-1)*dotdotX3*Y(i-1) -( MX3(i-1)*dotdotX2 + dotdotX3*MX2(i-1) + 2*dotX2*dotX3 ));
    dotdotdotX3 = a1 + a2 + a3;
    
    
    %----------------------Agent 1----------------------%
    %--- Solving 1st Linear Differential Equation
    MX1(i) = MX1(i-1) + T*dotX1 + T^2/2* dotdotX1 + T^3/6*dotdotdotX1;
    %---------------------------------------------------%
    
    %----------------------Agent 2----------------------%
    %---Solving 2nd  Linear Differential Equation
    MX2(i) = MX2(i-1) + T*dotX2 + T^2/2* dotdotX2+ T^3/6*dotdotdotX2;
    %---------------------------------------------------%
    
    %----------------------Agent 3----------------------%
    %---Solving 3rd  Linear Differential Equation
    MX3(i) = MX3(i-1) + T*dotX3 + T^2/2* dotdotX3+ T^3/6*dotdotdotX3;
    %---------------------------------------------------%
end
end

