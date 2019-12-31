design = [0.25    0.25    0.25    0.25];
% p =[0.6300    0.2700    0.1000];
% theta = [2.7013    1.2000    0.7000    5.4000
%     2.7004    1.2000    0.7000    5.4000
%     5.0000    0.8000    0.7000    1.0000];

% design = [0.4226    0.1745    0.1227    0.2802];
p = [0.0000    0.0000    0.0000    0.0006    0.0001 0.0101    0.0015    0.1852    0.0132    0.3511 0.0181    0.0387    0.1732    0.0255    0.1514 0.0312]

theta = [1.0000    1.0000    1.0000    1.0000
1.5000    1.2449    1.5000    1.5000
0.3477    1.5000    1.5000    1.5000
1.5000    1.5000    1.5000    1.5000
1.5000    0.9571    1.5000    1.5000
1.5000    1.4250    1.5000    1.5000
1.5000    1.5000    1.5000    1.5000
1.5000    1.4304    1.5000    1.5000
1.5000    1.5000    1.5000    1.5000
1.5000    1.4300    1.5000    1.5000
1.5000    1.4999    1.5000    1.5000
1.5000    1.4519    1.5000    1.5000
1.5000    1.4282    1.5000    1.5000
1.5000    1.4606    1.5000    1.5000
1.5000    1.4421    1.5000    1.5000
1.5000    1.4176    1.5000    1.5000];1.0002


rho = 0.5;

l= size(p);
size(theta);
a = [];
for i = 1:l(2)
    a(i) = p(i)*trace_psi(design,theta(i,:),rho);
end
sum(a)