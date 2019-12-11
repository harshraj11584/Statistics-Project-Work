function y = trace_psi(design,theta,rho)
    p = 3; %periods
%     rho = 0.5; %Cross-Correlation Coefficient

    % Xj for j= 1,2,3,4,5,6,7,8 corresponds to X_AAA, X_AAB, X_ABA, X_ABB, X_BAA, X_BAB, X_BBA, X_BBB respectively
    X1 = [ ones(p,1) [0;1;0] [0;0;1] [1;1;1] [0;1;1] ]; %X_AAA
    X2 = [ ones(p,1) [0;1;0] [0;0;1] [1;1;-1] [0;1;1] ]; %X_AAB
    X3 = [ ones(p,1) [0;1;0] [0;0;1] [1;-1;1] [0;1;-1] ]; %X_ABA
    X4 = [ ones(p,1) [0;1;0] [0;0;1] [1;-1;-1] [0;1;-1] ]; %X_ABB
    X5 = [ ones(p,1) [0;1;0] [0;0;1] [-1;1;1] [0;-1;1] ]; %X_BAA
    X6 = [ ones(p,1) [0;1;0] [0;0;1] [-1;1;-1] [0;-1;1] ]; %X_BAB
    X7 = [ ones(p,1) [0;1;0] [0;0;1] [-1;-1;1] [0;-1;-1] ]; %X_BBA
    X8 = [ ones(p,1) [0;1;0] [0;0;1] [-1;-1;-1] [0;-1;-1] ]; %X_BBB

    % parameters in one column
    % theta = [mu;beta;tau;gam];
       
    % formula for Linear Predictor : n_j = Xj * theta
    n1 = X1*theta';
    n2 = X2*theta';
    n3 = X3*theta';
    n4 = X4*theta';
    n5 = X5*theta';
    n6 = X6*theta';
    n7 = X7*theta';
    n8 = X8*theta';

    % using inverse logit link to obtain u_ij. Formula : g(u_ij) = n_ij |
    % where u_ij is mean of response variable
    u1 = exp(n1)./(1.+exp(n1));
    u2 = exp(n2)./(1.+exp(n2));
    u3 = exp(n3)./(1.+exp(n3));
    u4 = exp(n4)./(1.+exp(n4));
    u5 = exp(n5)./(1.+exp(n5));
    u6 = exp(n6)./(1.+exp(n6));
    u7 = exp(n7)./(1.+exp(n7));
    u8 = exp(n8)./(1.+exp(n8));

    % Aj = diag ( Var(Y1j), . . . , Var(Ypj) )
    % As Y_ij is Bernoulli Variable for our case, Var(Y_ij) = u_ij * (1 - u_ij)
    % Define Aj in line after eqn 3. For Bernoulli Case, Aj = Dj as defined on pg 38 line  3
    A1 = diag(u1.*(1.-u1));
    A2 = diag(u2.*(1.-u2));
    A3 = diag(u3.*(1.-u3));
    A4 = diag(u4.*(1.-u4));
    A5 = diag(u5.*(1.-u5));
    A6 = diag(u6.*(1.-u6));
    A7 = diag(u7.*(1.-u7));
    A8 = diag(u8.*(1.-u8));

    % Defining Correlation Matrix R(alpha)
    R = rho*ones(p)+(1-rho)*eye(p);
    % Storing Inverse of R(alpha)
    R_inv = inv(R);
    
    % "Length of design"
    design;
    A = design(1)*X1'*sqrt(A1)*R_inv*sqrt(A1)*X1 + design(2)*X2'*sqrt(A2)*R_inv*sqrt(A2)*X2 + design(3)*X3'*sqrt(A3)*R_inv*sqrt(A3)*X3 + design(4)*X4'*sqrt(A4)*R_inv*sqrt(A4)*X4  + design(5)*X5'*sqrt(A5)*R_inv*sqrt(A5)*X5 + design(6)*X6'*sqrt(A6)*R_inv*sqrt(A6)*X6 + design(7)*X7'*sqrt(A7)*R_inv*sqrt(A7)*X7 + design(8)*X8'*sqrt(A8)*R_inv*sqrt(A8)*X8;
   
    
    Q = sqrt(A8)*R_inv*sqrt(A8);
    K = X8'*Q*X8;
    size(A)
    B = eye(5)/(A);
    L = [0 0 0 1 0]';
    C = eye(1)/(L'*B*L);
    %y = log(C);
    y = trace(B*L*C*L'*B*K);
    
    
end