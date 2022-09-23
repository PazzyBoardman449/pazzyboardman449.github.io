function [tvec,ymat] = ProteinRegSolver(M0,P0)

%Parameters:
v = 3.26;
k1 = 0.045;
k2 = 0.161;
k3 = 0.869;
k7 = 2.174;
Ka = 5.5;
Kb = 15;
m = 3;
n = 2;

%Initial Conditions:

M0 = 0;
P0 = 0;

%ODE Solver Options:
options = odeset('MaxStep',1e-2);

%Computing The Exact Solution:

[tvec,ymat] = ode45(@(t,y) ProteinRegSolverVE(t,y,v,k1,k2,k3,k7,Ka,Kb,m,n))


%Vector Field File for CC2015 Model 2 System

function dydt = ProteinRegSolverVE(t,y,v,k1,k2,k3,k7,Ka,Kb,m,n)
 dydt = zeros(2,1);
 dydt(1) = -k1*y(1) + v/(1+(y(2)/Ka)^m) %dMdt
 dydt(2) = k2*y(1) + (k7*y(1)*y(3)^n)/(Kb^n+P^n) - k3*y(3) %dPdt
 