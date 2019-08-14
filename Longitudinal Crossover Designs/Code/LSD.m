function y = LSD(ps,mu,alpha,tau,gam,rho,tau1,tau2,tau3,tau4)
p = 4;


mu = 2.0;
alpha = [1.0 2.0 3.0];
tau = [1.0 2.0 3.0];
gam = [1.0 2.0 3.0];
sigmb = 3.0;
rho = 0.2;
ps = [1/4 1/4 1/4 1/4];

tau1 = [0 0 0;0 1 0;1 0 0;0 0 1];
tau2 = [0 0 1;1 0 0;0 0 0;0 1 0];
tau3 = [0 1 0;0 0 0;0 0 1;1 0 0];
tau4 = [1 0 0;0 0 1;0 1 0;0 0 0];


GM = [0 zeros(1,p-1); eye(p-1) zeros(p-1,1)];
gam1 = GM*tau1;
gam2 = GM*tau2;
gam3 = GM*tau3;
gam4 = GM*tau4;
beta = [mu alpha tau gam]';
per = [0 0 0;0 0 1;0 1 0;1 0 0];
X1 = [ones(p,1) per tau1 gam1];
X2 = [ones(p,1) per tau2 gam2];
X3 = [ones(p,1) per tau3 gam3];
X4 = [ones(p,1) per tau4 gam4];


Z1 = X1*beta;
eta11 = exp(Z1)./(1+exp(Z1));
eta12 = 1-eta11;
etad1 = eta12.*eta11;

Z2 = X2*beta;
eta21 = exp(Z2)./(1+exp(Z2));
eta22 = 1-eta21;
etad2 = eta21.*eta22;

Z3 = X3*beta;
eta31 = exp(Z3)./(1+exp(Z3));
eta32 = 1-eta31;
etad3 = eta31.*eta32;

Z4 = X4*beta;
eta41 = exp(Z4)./(1+exp(Z4));
eta42 = 1-eta41;
etad4 = eta41.*eta42;




W1 = diag(sqrt(etad1));
W2 = diag(sqrt(etad2));
W3 = diag(sqrt(etad3));
W4 = diag(sqrt(etad4));

% V1 = sigb*W1*W1;
% V2 = sigb*W2*W2;
% V3 = sigb*W3*W3;


%R = rho*ones(p)+(1-rho)*eye(p);
R = [1 rho rho^2 rho^3;rho 1 rho rho^2;rho^2 rho 1 rho;rho^3 rho^2 rho 1]
R1 = eye(p)/(R);
A1 = W1*R1*W1;
A2 = W2*R1*W2;
A3 = W3*R1*W3;
A4 = W4*R1*W4;

%V = sigmb*eye(p);

A = ps(1)*X1'*A1*X1 + ps(2)*X2'*A2*X2 + ps(3)*X3'*A3*X3 + ps(4)*X4'*A4*X4
B = pinv(A)
C = [B(5:7,5) B(5:7,6) B(5:7,7)]
y = log(det(C));