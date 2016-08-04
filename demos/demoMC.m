figure;
subplot(2,2,1); 
plot(generateMarkovChain([0;0],eye(2),eye(2),100),'LineWidth',2);
title('$\Sigma=I$','interpreter','latex', 'FontSize',20);
subplot(2,2,2); 
plot(generateMarkovChain([0;0],eye(2),[0.2,0.0;0.0,0.6],100),'LineWidth',2);
title('$\Sigma  = \left[ {\begin{array}{*{20}{c}}{0.2}&0 \\ 0&{0.6}\end{array}} \right]$','interpreter','latex','FontSize',20);
subplot(2,2,3); 
plot(generateMarkovChain([0;0],eye(2),[0.2,0.1;0.1,0.6],100),'LineWidth',2);
title('$\Sigma  = \left[ {\begin{array}{*{20}{c}}{0.2}&0.1 \\ 0.1&{0.6}\end{array}} \right]$','interpreter','latex','FontSize',20);

subplot(2,2,4); 
plot(generateMarkovChain([0;0],eye(2),[0.2,0.3;0.3,0.6],100),'LineWidth',2);
title('$\Sigma  = \left[ {\begin{array}{*{20}{c}}{0.2}&0.3 \\ 0.3 &{0.6}\end{array}} \right]$','interpreter','latex','FontSize',20);
