%Mechanical Parameters
E = 100 %Youngs Modulus/kPa
mu = 30 %Viscosity Pa•s
sigma0 = 100 %Constant Stress
epsilon0 = sigma0/E %Initial Strain

F_Max_Const_Stress = @(t,epsilon) sigma0/mu %Since stress is constant, maxwell equation reduces to this
F_Max_Const_Strain = @(t,epsilon) sigma0/mu %Since stress is constant, maxwell equation reduces to this


t0 = 0 %Initial Time
tfinal = 10 %Final Time
h = tfinal/1000 %Step Size

tout = [t0:h:tfinal];

epsilon = epsilon0;
epsilonout = epsilon;
    for t = t0 : h : tfinal-h
         s1 = F(t,epsilon);
         s2 = F(t+h/2, epsilon+h*s1/2);
         s3 = F(t+h/2, epsilon+h*s2/2);
         s4 = F(t+h, epsilon+h*s3);
         epsilon = epsilon + h*(s1 + 2*s2 + 2*s3 + s4)/6;
         epsilonout = [epsilonout; epsilon];
    end
    
    plot(tout,epsilonout,'b-')
    xlabel('Time, t/s')
    ylabel('Strain, ε(t)')
    title('Constant stress σ_0 = 100kPa')