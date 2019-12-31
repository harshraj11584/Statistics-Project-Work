% H-ALGORITHM : 
warning off;

p=3; %if changing, must change values of p everywhere in code 
rho = 0.5; %Cross-Correlation Coefficient 

global iteration_num
iteration_num=0;

%Step 0

l=1; 
B_pi_l = -inf;
h_grid_space = 0.25;
H = 0.0:h_grid_space:1.0;

%theta = [mu;beta;tau;gam];
"Initial Prior on Theta : "
theta_vals_l = [ ones(1,1+p) ]
prob_vals_l = [1.0] 

% Initial Prior PMF is (theta_vals, prob_vals)
"Starting Design"
design_l = get_optimal_on_the_average_design(theta_vals_l,prob_vals_l,rho)
B_pi_l = B(design_l, theta_vals_l, prob_vals_l,rho);

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,rho);

function step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,rho)
    global iteration_num;
    
    p=3;
    theta_argmax_psi = zeros(1,1+p); psi_max = -inf;
    stopping=true;
    
    for i=1:1:10        
        %initialization of theta_nod according to (-0.5,1.5) limit
        theta_nod = -2*rand(1,1+p)+(1.5);
        options = optimoptions('fmincon','Display','none');
        [this_theta_psi_max, this_psi_max] = fmincon(@(theta) -psi(design_l,theta,rho), theta_nod, [],[], [], [], -0.5*ones(1,1+p), 1.5*ones(1,1+p),[],options);
        this_psi_max = -1.0*this_psi_max;
        if this_psi_max > psi_max
            psi_max = this_psi_max;
            theta_argmax_psi = this_theta_psi_max;
        end
    end
    iteration_num
    iteration_num = iteration_num+1;
    psi_max
    B_pi_l
    if psi_max-B_pi_l>10^-4
        stopping=false;
    end
              
    if stopping==true
        "Found minmax design"
        design_l
        "Least Favourable Distribution"
        "probs"
        prob_vals_l
        "thetas"
        theta_vals_l
    else
        theta_l = theta_argmax_psi;
        step_2(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,theta_l,rho);
    end
end

function step_2(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,theta_l,rho)
    %delta_l = unit mass to theta_l
    step_3(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,theta_l,rho);
end

function step_3(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,theta_l,rho)
    design_l1 = design_l;
    largest_B_pi_t_l1 = -inf;
    prob_vals_l1 = prob_vals_l;
    theta_vals_l1 = theta_vals_l;
    
    theta_vals_l;
    prob_vals_l;
    
    for t_val=2:1:length(H)-1
        theta_vals_t_l1 = [theta_vals_l; theta_l];
        %"old"
        %prob_vals_l1
        prob_vals_t_l1 = [(1-H(t_val)).*prob_vals_l H(t_val)];
        %"H(t_val)"
        %H(t_val)
        %"new"
        %prob_vals_t_l1
        
        design_t_l1 = get_optimal_on_the_average_design(theta_vals_t_l1, prob_vals_t_l1,rho);
        B_pi_t_l1 = B(design_t_l1, theta_vals_t_l1, prob_vals_t_l1,rho);        
        if largest_B_pi_t_l1<B_pi_t_l1
            largest_B_pi_t_l1 = B_pi_t_l1;
            design_l1 = design_t_l1;
            prob_vals_l1 = prob_vals_t_l1;
            theta_vals_l1 = theta_vals_t_l1;
        end
    end
    step_4(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,design_l1, theta_vals_l1, prob_vals_l1, theta_l,rho);    
end

function step_4(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,design_l1, theta_vals_l1, prob_vals_l1, theta_l,rho) 
    B_pi_l1 = B(design_l1, theta_vals_l1, prob_vals_l1,rho);
    if B_pi_l1 > B_pi_l
        B_pi_l = B_pi_l1;
        design_l = design_l1;
        "Assigned l1 to l";
        theta_vals_l = theta_vals_l1;
        prob_vals_l = prob_vals_l1; 
        step_1(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,rho);
    else
        h_grid_space = h_grid_space/2.0;
        H = 0.0:h_grid_space :1.0;
        step_3(H,h_grid_space,theta_vals_l,prob_vals_l,design_l,B_pi_l,theta_l,rho);
    end
end

% function ret_val = get_least_favourable_argument(design)
%     p=3;
%     arg_0 = 10*rand(1,1+p+1);options = optimoptions('fmincon','Display','none');
%     l_f_arg = fmincon(@(theta)-1* psi(design,theta), arg_0, [],[],[],options);
%     ret_val = l_f_arg;
% end

function ret_val = get_optimal_on_the_average_design(theta_vals, prob_vals,rho)
    p=3;
%     des_0 = [ 0.1 0.2 0.3 0.4 ]; 
    des_0 = rand(1,2^p); des_0 = des_0/sum(des_0);
    nseq=2^p;
    options = optimoptions('fmincon','Display','none');
    opt_design = fmincon(@(des)B(des,theta_vals, prob_vals,rho), des_0, [],[], ones(1,nseq), [1.0], zeros(nseq,1), ones(nseq,1),[],options);
    ret_val = opt_design;
end

function ret_val = B (design, theta_vals, prob_vals,rho)
    s = 0;
    for i = 1:length(prob_vals)
        s = s + prob_vals(i)* psi( design, theta_vals(i,:),rho );
    end
    ret_val=s;
end

function y = psi(design,theta,rho)
    p = 3; %periods
    

    % Xj for j= 1,2,3,4,5,6,7,8 corresponds to X_AAA, X_AAB, X_ABA, X_ABB, X_BAA, X_BAB, X_BBA, X_BBB respectively
    X1 = [ ones(p,1) [0;1;0] [0;0;1] [1;1;1]  ]; %X_AAA
    X2 = [ ones(p,1) [0;1;0] [0;0;1] [1;1;-1]  ]; %X_AAB
    X3 = [ ones(p,1) [0;1;0] [0;0;1] [1;-1;1]  ]; %X_ABA
    X4 = [ ones(p,1) [0;1;0] [0;0;1] [1;-1;-1]  ]; %X_ABB
    X5 = [ ones(p,1) [0;1;0] [0;0;1] [-1;1;1]  ]; %X_BAA
    X6 = [ ones(p,1) [0;1;0] [0;0;1] [-1;1;-1]  ]; %X_BAB
    X7 = [ ones(p,1) [0;1;0] [0;0;1] [-1;-1;1]  ]; %X_BBA
    X8 = [ ones(p,1) [0;1;0] [0;0;1] [-1;-1;-1] ]; %X_BBB

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
    
    %print("Length of design")
    design;
    asymptotic_information_matrix = design(1)*X1'*sqrt(A1)*R_inv*sqrt(A1)*X1 + design(2)*X2'*sqrt(A2)*R_inv*sqrt(A2)*X2 + design(3)*X3'*sqrt(A3)*R_inv*sqrt(A3)*X3 + design(4)*X4'*sqrt(A4)*R_inv*sqrt(A4)*X4  + design(5)*X5'*sqrt(A5)*R_inv*sqrt(A5)*X5 + design(6)*X6'*sqrt(A6)*R_inv*sqrt(A6)*X6 + design(7)*X7'*sqrt(A7)*R_inv*sqrt(A7)*X7 + design(8)*X8'*sqrt(A8)*R_inv*sqrt(A8)*X8;
    %n = 1e12;
    %asymptotic_information_matrix = n.* asymptotic_information_matrix;
    asymptotic_variance = inv(asymptotic_information_matrix);
    size(asymptotic_variance);
    C = log(asymptotic_variance(4,4));
    y = (C);
end