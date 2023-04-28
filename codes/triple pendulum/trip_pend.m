syms alpha1(t) alpha2(t) alpha3(t) Y
syms q [3 1]
syms Lag [3 1]
syms eqns [3 1]

q(1) = alpha1;
q(2) = alpha2;
q(3) = alpha3;
dq = diff(q,t);

m1 = 0.5;
m2 = 0.15;
m3 = 0.05;
l1 = 0.1;
l2 = 0.1;
l3 = 0.1;
g = 9.8;

% REMEBER alpha is the angle with the vertical
x1 = l1*sin(alpha1);
y1 = -l1*cos(alpha1);
x2 = l1*sin(alpha1) + l2*sin(alpha2);
y2 = -l1*cos(alpha1) - l2*cos(alpha2);
x3 = l1*sin(alpha1) + l2*sin(alpha2) + l3*sin(alpha3);
y3 = -l1*cos(alpha1) - l2*cos(alpha2) - l3*cos(alpha3);

T = 0.5*(m1*((diff(x1,t))^2 + (diff(y1,t))^2) + m2*((diff(x2,t))^2 + (diff(y2,t))^2) + m3*((diff(x3,t))^2 + (diff(y3,t))^2));
V = m1*g*y1 + m2*g*y2 + m3*g*y3;
L = T-V;

for i=1:3
    Lag(i) = diff(diff(L,dq(i)),t) - diff(L,q(i)) == 0;
end
for i=1:3
    eqns(i) = simplify(lhs(Lag(i)) - rhs(Lag(i)));
end
eqns == 0 %#ok<EQEFF> 

[VF,Subs] = odeToVectorField(Lag)
odefcn = matlabFunction(VF, 'Vars',{t,Y})
[X,Y] = ode45(odefcn, [0 5], [pi/2 0 pi/2 0 pi 0]);

plot(X, Y)
grid
legend(string(Subs))