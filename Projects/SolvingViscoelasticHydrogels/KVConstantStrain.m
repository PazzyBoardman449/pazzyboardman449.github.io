%Mechanical Parameters
E = 100 %Youngs Modulus/kPa
mu = 30 %Viscosity Pa•s
epsilon0 = 2 %Constant Strain
sigma0 = E*epsilon0 %Initial Stress

t0 = 0 %Initial Time
tfinal = 10 %Final Time
h = tfinal/1000 %Step Size

tout = [t0:h:tfinal]
sigma0vec = sigma0*ones(length(tout));

plot(tout,sigma0vec,'r-')
xlabel('Time, t/s')
ylabel('Stress σ/kPa')
title('Constant strain, ε_0 = 2')
%axis([0 5 0 4])