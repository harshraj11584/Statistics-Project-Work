function y = trace_psi(design,theta,rho)
    p = 4; %periods
%     rho = 0.5; %Cross-Correlation Coefficient

    % Xj for j= 1,2,...,16 corresponding treatment has been written beside
    X1 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;1;1;1]  ]; %X_AAAA
    X2 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;1;1;-1]  ]; %X_AAAB
    X3 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;1;-1;1]  ]; %X_AABA
    X4 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;1;-1;-1]  ]; %X_AABB
    X5 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;-1;1;1]  ]; %X_ABAA
    X6 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;-1;1;-1]  ]; %X_ABAB
    X7 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;-1;-1;1]  ]; %X_ABBA
    X8 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [1;-1;-1;-1]  ]; %X_ABBB
    X9 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;1;1;1]  ]; %X_BAAA
    X10 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;1;1;-1]  ]; %X_BAAB
    X11 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;1;-1;1]  ]; %X_BABA
    X12 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;1;-1;-1]  ]; %X_BABB
    X13 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;-1;1;1]  ]; %X_BBAA
    X14 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;-1;1;-1]  ]; %X_BBAB
    X15 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;-1;-1;1]  ]; %X_BBBA
    X16 = [ ones(p,1) [0;0;0;1] [0;0;1;0] [0;1;0;0] [-1;-1;-1;-1] ]; %X_BBBB

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
    n9 = X9*theta';
    n10 = X10*theta';
    n11 = X11*theta';
    n12 = X12*theta';
    n13 = X13*theta';
    n14 = X14*theta';
    n15 = X15*theta';
    n16 = X16*theta';

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
    u9 = exp(n9)./(1.+exp(n9));
    u10 = exp(n10)./(1.+exp(n10));
    u11 = exp(n11)./(1.+exp(n11));
    u12 = exp(n12)./(1.+exp(n12));
    u13 = exp(n13)./(1.+exp(n13));
    u14 = exp(n14)./(1.+exp(n14));
    u15 = exp(n15)./(1.+exp(n15));
    u16 = exp(n16)./(1.+exp(n16));

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
    A9 = diag(u9.*(1.-u9));
    A10 = diag(u10.*(1.-u10));
    A11 = diag(u11.*(1.-u11));
    A12 = diag(u12.*(1.-u12));
    A13 = diag(u13.*(1.-u13));
    A14 = diag(u14.*(1.-u14));
    A15 = diag(u15.*(1.-u15));
    A16 = diag(u16.*(1.-u16));

    % Defining Correlation Matrix R(alpha)
    R = rho*ones(p)+(1-rho)*eye(p);
    % Storing Inverse of R(alpha)
    R_inv = inv(R);
    
    % "Length of design"
    design;
    A = design(1)*X1'*sqrt(A1)*R_inv*sqrt(A1)*X1 + design(2)*X2'*sqrt(A2)*R_inv*sqrt(A2)*X2 + design(3)*X3'*sqrt(A3)*R_inv*sqrt(A3)*X3 + design(4)*X4'*sqrt(A4)*R_inv*sqrt(A4)*X4  + design(5)*X5'*sqrt(A5)*R_inv*sqrt(A5)*X5 + design(6)*X6'*sqrt(A6)*R_inv*sqrt(A6)*X6 + design(7)*X7'*sqrt(A7)*R_inv*sqrt(A7)*X7 + design(8)*X8'*sqrt(A8)*R_inv*sqrt(A8)*X8 + design(9)*X9'*sqrt(A9)*R_inv*sqrt(A9)*X9 + design(10)*X10'*sqrt(A10)*R_inv*sqrt(A10)*X10+ design(11)*X11'*sqrt(A11)*R_inv*sqrt(A11)*X11+ design(12)*X12'*sqrt(A12)*R_inv*sqrt(A12)*X12+ design(13)*X13'*sqrt(A13)*R_inv*sqrt(A13)*X13+ design(14)*X14'*sqrt(A14)*R_inv*sqrt(A14)*X14+ design(15)*X15'*sqrt(A15)*R_inv*sqrt(A15)*X15+ design(16)*X16'*sqrt(A16)*R_inv*sqrt(A16)*X16;
 
    
    Q = sqrt(A16)*R_inv*sqrt(A16);
    K = X16'*Q*X16;
    size(A);
    B = eye(5)/(A);
    L = [0 0 0 0 1]';
    C = eye(1)/(L'*B*L);
    %y = log(C);
    y = trace(B*L*C*L'*B*K);    
    
end