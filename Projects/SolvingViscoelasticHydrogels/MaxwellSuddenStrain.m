%Mechanical Parameters
E = 100 %Youngs Modulus/kPa
mu = 30 %Viscosity Pa•s
epsilon0 = 2 %Constant Strain
sigma0 = E*epsilon0 %Initial Stress

F = @(t,sigma) -(E/mu)*sigma %Since strain is constant, maxwell equation reduces to this

t0 = 0 %Initial Time
tfinal = 10 %Final Time
h = tfinal/1000 %Step Size

tout = [t0:h:tfinal];

sigma = sigma0;
sigmaout = sigma;
    for t = t0 : h : tfinal-h
         s1 = F(t,sigma);
         s2 = F(t+h/2, sigma+h*s1/2);
         s3 = F(t+h/2, sigma+h*s2/2);
         s4 = F(t+h, sigma+h*s3);
         sigma = sigma + h*(s1 + 2*s2 + 2*s3 + s4)/6;
         sigmaout = [sigmaout; sigma];
    end
    
    plot(tout,sigmaout,'r-')
    xlabel('Time, t/s')
    ylabel('Stress, σ/kPa')
    title('Sudden Strain ε_0 = 2')
    axis([0 5 0 200])