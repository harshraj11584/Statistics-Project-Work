ps0 = [ 0.1 0.2 0.3 0.4 ]
res = allocation_ps(ps0)
nseq = 4;
[ps_final,fmin_val] = fmincon(@(prob)allocation_ps(prob), ps0, eye(nseq),ones(nseq,1), ones(1,nseq), [1.0], zeros(nseq,1), ones(nseq,1))


function z = allocation_ps(ps)
%Defining Parameters
ps = [0.1 0.7 0.1 0.1];
mu = [2.0];
beta = [1;3];
tau = [5.0];
gam = [6.0];
rho = 0.3;
res = get_allocation(ps, mu,beta,tau,gam, rho);
z = (res);
end

function y = get_allocation(ps,mu,beta,tau,gam, rho)

p = 2; %periods
t = 2; %treatments

% parameters in one column
theta = [mu;beta;tau;gam];

% Xj for j= 1,2,3,4 corresponds to X_AA,X_AB,X_BA,X_BB respectively
X1 = [ ones(p,1) eye(p) [1;1] [0;1] ];
X2 = [ ones(p,1) eye(p) [1;-1] [0;1] ];
X3 = [ ones(p,1) eye(p) [-1;1] [0;-1] ];
X4 = [ ones(p,1) eye(p) [-1;-1] [0;-1] ];

% formula for Linear Predictor : n_j = Xj * theta
n1 = X1*theta;
n2 = X2*theta;
n3 = X3*theta;
n4 = X4*theta;

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

asymptotic_information_matrix = ps(1)*X1'*sqrt(A1)*R_inv*sqrt(A1)*X1 + ps(2)*X2'*sqrt(A2)*R_inv*sqrt(A2)*X2 + ps(3)*X3'*sqrt(A3)*R_inv*sqrt(A3)*X3 + ps(4)*X4'*sqrt(A4)*R_inv*sqrt(A4)*X4
%n = 1e12;
%asymptotic_information_matrix = n.* asymptotic_information_matrix;
asymptotic_variance = inv(asymptotic_information_matrix)
C = asymptotic_variance(3,3);
y = (C);
end


