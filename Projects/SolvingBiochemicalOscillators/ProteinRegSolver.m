options = odeset('RelTol', 1e-4, 'NonNegative', [1 2]);
[t1,y1] = ode45('ProteinReg', [0 250], [0 0], options);
[t2,y2] = ode45('ProteinReg', [0 250], [10 10], options);
[t3,y3] = ode45('ProteinReg', [0 250], [0 30], options);


plot(t1,y1,'LineWidth',1.5);
legend('mRNA', 'Protein');
xlabel('Time [h]')
ylabel('Concentration [a.u]')
title('Plot of mRNA and Protein Concentration against time')

subplot(3,1,1)
plot(t1,y1,'LineWidth',1.5);
legend('mRNA', 'Protein');
xlabel('Time [h]')
ylabel('Concentration [a.u]')
title('Plot of mRNA and Protein Concentration against time')

subplot(3,1,2)
plot(t2,y2,'LineWidth',1.5);
legend('mRNA', 'Protein');
xlabel('Time [h]')
ylabel('Concentration [a.u]')
title('Plot of mRNA and Protein Concentration against time')

subplot(3,1,3)
plot(t3,y3,'LineWidth',1.5);
legend('mRNA', 'Protein');
xlabel('Time [h]')
ylabel('Concentration [a.u]')
title('Plot of mRNA and Protein Concentration against time')

