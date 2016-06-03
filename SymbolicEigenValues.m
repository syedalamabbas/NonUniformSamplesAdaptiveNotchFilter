
clear
clc

%% Symbolic Computations only
syms theta_0 eta gamma
disp('This is the Jacobian:')

J = [ 0  1  0  0 ;  -theta_0^2  -2*theta_0*eta  theta_0*sin(theta_0)+ theta_0*cos(theta_0)/eta  theta_0^2 ; 0  -gamma*theta_0*cos(theta_0)  gamma*theta_0*sin(2*theta_0)/(4*eta)  gamma*theta_0^2*cos(theta_0)/(2*eta) ;  -2*theta_0*eta  0   cos(theta_0)  0   ];

J

disp('Doing the Eigen value analysis:')

simp_Eig = simplify(eig(J))

disp('Doing the substiution and simplification further:')

syms alpha 

disp('alpha = 2^(1/2)*(128*eta^4 - gamma^2*cos(4*theta_0) + gamma^2 - 64*eta*gamma*cos(theta_0)^2 - 32*eta^2*gamma*sin(2*theta_0))^(1/2)');

subs_Eig = subs(simp_Eig,  2^(1/2)*(128*eta^4 - gamma^2*cos(4*theta_0) + gamma^2 - 64*eta*gamma*cos(theta_0)^2 - 32*eta^2*gamma*sin(2*theta_0))^(1/2) ,alpha)

%---------------Alpha solving---------------%
disp('Solving alpha for gamma')

alpha = 2^(1/2)*(128*eta^4 - gamma^2*cos(4*theta_0) + gamma^2 - 64*eta*gamma*cos(theta_0)^2 - 32*eta^2*gamma*sin(2*theta_0))^(1/2);
simplify(solve(alpha, gamma))


%% All about Substitution
f_0 = [30 : 3: 300];

lambda_3 = simp_Eig(3);
lambda_4 = simp_Eig(4);

disp('Simplifying the product of two eigen values ')
simplify(lambda_3*lambda_4)
% theta_0
lambda_3 = subs(lambda_3,theta_0, 2*pi*f_0);  
lambda_4 = subs(lambda_4,theta_0, 2*pi*f_0); 

% xi
lambda_3 = subs(lambda_3,eta, .1);  
lambda_4 = subs(lambda_4, eta, .1); 

% gamma
lambda_3 = subs(lambda_3, gamma, .001);  
lambda_4 = subs(lambda_4, gamma, .001);

%% Plotting the Eigen values !
figure, plot(lambda_4, lambda_3, '+' )
xlabel('\lambda_3')
ylabel('\lambda_4')
