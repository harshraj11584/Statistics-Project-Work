function y = trace_psi(design,theta,rho)
    p = 2; %periods
%     rho = 0.5; %Cross-Correlation Coefficient

    % Xj for j= 1,2,3,4 corresponds to X_AA,X_AB,X_BA,X_BB respectively
    X1 = [ ones(p,1) [1;-1] [1;1]  ];
    X2 = [ ones(p,1) [1;-1] [1;-1]  ];
    X3 = [ ones(p,1) [1;-1] [-1;1]  ];
    X4 = [ ones(p,1) [1;-1] [-1;-1]  ];

    % parameters in one column
    % theta = [mu;beta;tau;gam];
       
    % formula for Linear Predictor : n_j = Xj * theta
    n1 = X1*theta';
    n2 = X2*theta';
    n3 = X3*theta';
    n4 = X4*theta';

    % using inverse logit link to obtain u_ij. Formula : g(u_ij) = n_ij |
    % where u_ij is mean of response variable
    u1 = exp(n1)./(1.+exp(n1));
    u2 = exp(n2)./(1.+exp(n2));
    u3 = exp(n3)./(1.+exp(n3));
    u4 = exp(n4)./(1.+exp(n4));

    % Aj = diag ( Var(Y1j), . . . , Var(Ypj) )
    % As Y_ij is Bernoulli Variable for our case, Var(Y_ij) = u_ij * (1 - u_ij)
    % Define Aj in line after eqn 3. For Bernoulli Case, Aj = Dj as defined on pg 38 line  3
    A1 = diag(u1.*(1.-u1));
    A2 = diag(u2.*(1.-u2));
    A3 = diag(u3.*(1.-u3));
    A4 = diag(u4.*(1.-u4));

    % Defining Correlation Matrix R(alpha)
    R = rho*ones(p)+(1-rho)*eye(p);
    % Storing Inverse of R(alpha)
    R_inv = inv(R);
    
    %print("Length of design")
    design;
    A = design(1)*X1'*sqrt(A1)*R_inv*sqrt(A1)*X1 + design(2)*X2'*sqrt(A2)*R_inv*sqrt(A2)*X2 + design(3)*X3'*sqrt(A3)*R_inv*sqrt(A3)*X3 + design(4)*X4'*sqrt(A4)*R_inv*sqrt(A4)*X4;
    
    
    Q = sqrt(A1)*R_inv*sqrt(A1);
    K = X1'*Q*X1;
    B = eye(3)/(A);
    L = [0 0 1]';
    C = eye(1)/(L'*B*L);
    %y = log(C);
    y = trace(B*L*C*L'*B*K);
    
    
end