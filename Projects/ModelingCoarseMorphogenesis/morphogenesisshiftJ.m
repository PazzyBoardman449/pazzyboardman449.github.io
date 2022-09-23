%% Solution to the Basic Morphogenisis equation, But in the fixed time scenario with different fluxes.

%Specifying parameters
Dt = 1; %Timestep (s)
Dx = 2.5; %Lattice Step (µm)
D = 1; %Diffusion Constant (µm•m^2•s^-1)
mu = 10e-4; %Degredation Rate (s^-1)
J = 0.3; %Flux at the boundary (µm^-2 s^-1)
J2 = 0.5; %Flux at the boundary (µm^-2 s^-1)
J3 = 1; %Flux at the boundary (µm^-2 s^-1)
L = 500; %Length of the system (m)
Tmax = 50*60; %Total Duration of the System (s)
pmax = 5; %Threshold Value For Concentration (µm^-3)

%Variables
num_steps = Tmax/Dt; %Number of time steps over the full period
num_latt_pts = L/Dx; %Total Number of Lattice Steps

p = zeros(num_steps,num_latt_pts); %Array containig concentration profile across the length of the container for the total duration of the simulation
p2 = zeros(num_steps,num_latt_pts); %Array containig concentration profile across the length of the container for the total duration of the simulation
p3 = zeros(num_steps,num_latt_pts); %Array containig concentration profile across the length of the container for the total duration of the simulation

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
    %p(i,num_latt_pts) = p(i,num_latt_pts-1)+J*Dx/D;

end 

for k = 2:num_steps
    
    %Finite Difference Version of the Diffusion Equation with Losses
    for l = 2:num_latt_pts-1
       p2(k,l) = p2(k-1,l) + D*Dt/Dx^2 *(p2(k-1,l-1)-2*p2(k-1,l)+p2(k-1,l+1))- mu*Dt*p2(k-1,l);
           
    end
    
    %Impose Flux Boundary Condition at x=0
    p2(k,1) = p2(k,2)+J2*Dx/D; 
    %p2(k,num_latt_pts) = p(k,num_latt_pts-1)+J2*Dx/D;

end 

for m = 2:num_steps
    
    %Finite Difference Version of the Diffusion Equation with Losses
    for n = 2:num_latt_pts-1
       p3(m,n) = p3(m-1,n) + D*Dt/Dx^2 *(p3(m-1,n-1)-2*p3(m-1,n)+p3(m-1,n+1))- mu*Dt*p3(m-1,n);
           
    end
    
    %Impose Flux Boundary Condition at x=0
    p3(m,1) = p3(m,2)+J3*Dx/D;
    %p3(m,num_latt_pts) = p(m,num_latt_pts-1)+J3*Dx/D;

end 


%Returning Value of x for which Concentration exceeds Threshold:%


% ===== Plotting =====

%Plot Of Concentration Against Position

x = (0:(num_latt_pts-1))*Dx;
pmaxvec = pmax*ones(num_latt_pts,1);

plot(x,p(3000,:),'r','linewidth',1)
hold on
plot(x,p2(3000,:),'b','linewidth',1)
hold on
plot(x,p3(3000,:),'g','linewidth',1)
hold on
plot(x,pmaxvec,'k-','linewidth',1)
%axis([0 250 0 10])
title('Morphogen Concentration at t = 100s for different fluxes','fontsize',16)
xlabel('Position x, µm','fontsize',12)
ylabel('Concentration ρ(x,t), µm^-3','fontsize',12)
legend('ρ(x,t=100s), J = 0.3','ρ(x,t=100s), J = 0.36','ρ(x,t=100s), J = 0.24','Threshold ρ_{max}','fontsize',12)