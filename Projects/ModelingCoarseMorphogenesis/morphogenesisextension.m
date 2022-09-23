%% Solution to the Basic Morphogenisis equation for a flux at both ends.

%Specifying parameters
Dt = 1; %Timestep (s)
Dx = 2.5; %Lattice Step (µm)
D = 1; %Diffusion Constant (µm•m^2•s^-1)
mu = 10e-4; %Degredation Rate (s^-1)
J = 0.3; %Flux at the boundary (µm^-2 s^-1)
L = 500; %Length of the system (m)
Tmax = 50*60; %Total Duration of the System (s)

%Variables
num_steps = Tmax/Dt; %Number of time steps over the full period
num_latt_pts = L/Dx; %Total Number of Lattice Steps

p = zeros(num_steps,num_latt_pts); %Array containig concentration profile across the length of the container for the total duration of the simulation

%Initial Condition:
%{Since concentration is zero everywhere there is no initial condition to be specified as all the aray elements are already zero. %}

%Main For Loop Over Time
for i = 2:num_steps
    
    %Finite Difference Version of the Diffusion Equation with Losses
    for j = 2:num_latt_pts-1
       p(i,j) = p(i-1,j) + D*Dt/Dx^2 *(p(i-1,j-1)-2*p(i-1,j)+p(i-1,j+1))- mu*Dt*p(i-1,j);
    end
    
    %Impose Flux Boundary Condition at x=0
    p(i,1) = p(i,2)+J*Dx/D;
    %Impose Flux Boundary Condition at x=L
    p(i,num_latt_pts) = p(i,num_latt_pts-1)+J*Dx/D;
    
end 

% ===== Plotting =====

%(Comment out the plot which ins't needed - Both are commented out by default.)

%Plot Of Concentration Against Position 

%{
x = (0:(num_latt_pts-1))*Dx;
plot(x,p(100,:),'r','linewidth',3)
hold on
plot(x,p(500,:),'b','linewidth',3)
hold on
plot(x,p(2000,:),'k','linewidth',3)
title('Plot of Morphogen Concentration Against Space','fontsize',16)
xlabel('Position x, µm','fontsize',12)
ylabel('Concentration ρ(x,t), µm^-3','fontsize',12)
legend('t = 100s','t = 500s','t = 2000s','fontsize',12)
%}

%Plot Of Concentration Against Time

%{
t = (0:(num_steps-1))*Dt;
plot(t,p(:,10),'r','linewidth',3)
hold on
plot(t,p(:,100),'b','linewidth',3)
hold on
plot(t,p(:,180),'k','linewidth',3)
title('Plot of Morphogen Concentration Against Time','fontsize',16)
xlabel('Time t, s','fontsize',12)
ylabel('Concentration ρ(x,t), µm^-3','fontsize',12)
legend('x = 25µm','x = 250µm','x = 450µm','fontsize',12,'location','northwest')
%}