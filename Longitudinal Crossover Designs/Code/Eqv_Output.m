design = [0.2341    0.1172    0.2485    0.4002];
p =[0.6300    0.2700    0.1000];
theta = [2.7013    1.2000    0.7000    5.4000
    2.7004    1.2000    0.7000    5.4000
    5.0000    0.8000    0.7000    1.0000];
l= size(p);
a = [];
for i = 1:l(2)
    a(i) = p(i)*trace_psi(design,theta(i,:))
end
sum(a)