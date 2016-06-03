function SolutionX = FunctionalANFAdomianSpecial(t,X0, SignalY, xi, gamma)
N = length(t);
SolutionX = zeros(N,length(X0));
SolutionX(1,:) = X0;   % Intial conditions given

for i = 2:N
    
    % Integrate from 0 to h
    t_min = t(i-1); % interval start
    t_max = t(i); % interval end
    
    %------------------ 1st Adomain Coefficients----------------------%
    % A_0 = x_2_0;
    % B_0 =  theta_0*x_2_0;
    % C_0 =  theta_0^2*x_1_0;
    % D_0 = theta_0^2;
    
    A_0 = SolutionX(i-1,2);
    B_0 = SolutionX(i-1,3)*SolutionX(i-1,2);
    C_0 = SolutionX(i-1,3)^2*SolutionX(i-1,1);
    D_0 = SolutionX(i-1,3)^2;
    E_0 = SolutionX(i-1,3)*SolutionX(i-1,1)*SolutionX(i-1,2);
    
    X1_1 = A_0 * (t_max - t_min); 
    X2_1 = (-2*xi*B_0-C_0+SignalY(i-1)*D_0)*(t_max - t_min);
    X3_1 = (gamma*(2*xi*E_0-SignalY(i-1)*C_0))*(t_max - t_min);
    
    
    %----------------- 2nd Adomain Coefficients -------------------------%
    % A_1 = x_2_1;
    % B_1 = (theta_0*x_2_1 + theta_1*x_2_0);
    % C_1 = (x_1_1*theta_0^2 + 2*theta_1*x_1_0*theta_0);
    % D_1 = (2*theta_0*theta_1);
    % E_1 = (theta_0*x_1_0*x_2_1 + theta_0*x_1_1*x_2_0 + theta_1*x_1_0*x_2_0);
    
    A_1 = X2_1;
    B_1 = SolutionX(i-1,3) * X2_1  + X3_1 * SolutionX(i-1,2);
    C_1 = X1_1*SolutionX(i-1,3)^2 + 2 * X3_1 * SolutionX(i-1,1) * SolutionX(i-1,3);
    D_1 = 2 * SolutionX(i-1,3) * X3_1;
    E_1 = SolutionX(i-1,3) * SolutionX(i-1,1)* X2_1 + SolutionX(i-1,3)* X1_1 *  SolutionX(i-1,1) + X3_1 * SolutionX(i-1,1)*SolutionX(i-1,2);
    
    X1_2 = A_1 * (t_max - t_min);  
    X2_2 = (-2*xi*B_1-C_1+SignalY(i-1)*D_1)*(t_max - t_min); 
    X3_2 = (gamma*(2*xi*E_1-SignalY(i-1)*C_1))*(t_max - t_min); 
    
    
    %-----------------3rd Adomain Coefficients --------------------------%
    % A_2  = x_2_2
    % B_2 = (theta_0*x_2_2 + theta_1*x_2_1 + theta_2*x_2_0);
    % C_2 = (x_1_2*theta_0^2 + 2*x_1_1*theta_0*theta_1 + 2*theta_2*x_1_0*theta_0 + x_1_0*theta_1^2);
    % D_2 = (theta_1^2 + 2*theta_0*theta_2);
    % E_2 = (theta_0*x_1_0*x_2_2 + theta_0*x_1_1*x_2_1 + theta_0*x_1_2*x_2_0 + theta_1*x_1_0*x_2_1 + theta_1*x_1_1*x_2_0 + theta_2*x_1_0*x_2_0);
    
    A_2 = X2_2;
    B_2 = SolutionX(i-1,3)*X2_2 + X3_1*X2_1 + X3_2* SolutionX(i-1,2);
    C_2 = X1_2*SolutionX(i-1,3)^2 + 2 * X1_1 * SolutionX(i-1,3)* X3_1 + 2 * X3_2 * SolutionX(i-1,1)*SolutionX(i-1,3) + SolutionX(i-1,1)*X3_1^2;
    D_2 = X3_1^2 + 2* SolutionX(i-1,3)* X3_2;
    E_2 = SolutionX(i-1,3)* SolutionX(i-1,1)* X2_2 + SolutionX(i-1,3)*X1_1* X2_1 + SolutionX(i-1,3)* X1_2 * SolutionX(i-1,2) + X3_1 * SolutionX(i-1,1)*  X2_1 + X3_1 * X2_1 * SolutionX(i-1,2) +  X3_2 * SolutionX(i-1,1) * SolutionX(i-1,2);
    
    
    X1_3 = A_2 * (t_max - t_min);   
    X2_3 = (-2*xi*B_2-C_2+SignalY(i-1)*D_2)*(t_max - t_min); 
    X3_3 = (gamma*(2*xi*E_2-SignalY(i-1)*C_2))*(t_max - t_min);  
    
    %---------------------4th Adomain Coefficients-----------------------%
    %     A_3 = x_2_3;
    %     B_3 = (theta_0*x_2_3 + theta_1*x_2_2 + theta_2*x_2_1 + theta_3*x_2_0);
    %     C_3 = (theta_1^2*x_1_1 + theta_0^2*x_1_3 + 2*theta_0*theta_1*x_1_2 + 2*theta_0*theta_2*x_1_1 + 2*theta_0*theta_3*x_1_0 + 2*theta_1*theta_2*x_1_0);
    %     D_3 = (2*theta_0*theta_3 + 2*theta_1*theta_2);
    %     E_3 = (theta_0*x_1_0*x_2_3 + theta_0*x_1_1*x_2_2 + theta_0*x_1_2*x_2_1 + theta_0*x_1_3*x_2_0 + theta_1*x_1_0*x_2_2 + theta_1*x_1_1*x_2_1 + theta_1*x_1_2*x_2_0 + theta_2*x_1_0*x_2_1 + theta_2*x_1_1*x_2_0 + theta_3*x_1_0*x_2_0);
    
    
    A_3 = X2_3;
    B_3 = SolutionX(i-1,3) * X2_3 + X3_1* X2_2 + X3_2* X2_1 + X3_3 * SolutionX(i-1,2) ;
    C_3 = X3_1^2*X1_1 + SolutionX(i-1,3)^2*X1_3 + 2* SolutionX(i-1,3)* X3_1 * X1_2 + 2* SolutionX(i-1,3)*  X3_2* X1_1 + 2* X3_1* X3_2* SolutionX(i-1,1);
    D_3 = 2* SolutionX(i-1,3)* X3_3 + 2 * X3_1* X3_2;
    E_3 = SolutionX(i-1,3)* SolutionX(i-1,3)*X2_3 + SolutionX(i-1,3)*X1_1*X2_2 +  SolutionX(i-1,3)*X1_2*X2_1 +  SolutionX(i-1,3)*X1_3*SolutionX(i-1,2)+  X3_1*SolutionX(i-1,1)*X2_2 + X3_1*X1_1*X2_1+ X3_1*X1_2*SolutionX(i-1,2)+ X3_3*SolutionX(i-1,1)*SolutionX(i-1,2);
    
    
    X1_4 = A_3 * (t_max - t_min);   
    X2_4 = (-2*xi*B_3-C_3+SignalY(i-1)*D_3)*(t_max - t_min); 
    X3_4 = (gamma*(2*xi*E_3-SignalY(i-1)*C_3))*(t_max - t_min);  
    
    SolutionX(i,1) =  SolutionX(i-1,1) + X1_1 + X1_2 + X1_3 + X1_4;         % Solution using truncated Adomain Polynomials
    SolutionX(i,2) =  SolutionX(i-1,2) + X2_1 + X2_2 + X2_3 + X2_4;
    SolutionX(i,3) =  SolutionX(i-1,3) + X3_1 + X3_2 + X3_3 + X3_4;
    
   
end
end



