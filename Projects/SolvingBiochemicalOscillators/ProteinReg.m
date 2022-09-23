function dy = ProteinReg(t,y)

%Parameters:
v = 0.2;
k1 = 0.045;
k2 = 0.161;
k3 = 0.869;
k7 = 2.174;
Ka = 5.5;
Kb = 15;
m = 3;
n = 2;

%Initial Conditions:
dy = [0;0];

dy(1) = -k1*y(1) + v/(1+(y(2)/Ka)^m);
dy(2) = k2*y(1) + (k7*y(1)*y(2)^n)/(Kb^n+y(2)^n) - k3*y(2);