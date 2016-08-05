
function [] = plotgraphs4(alphas,thetas,residual)

global fignum

% This function plots graphs of the values of alpha, theta and the
% residual.    



% Graph 1: The graph of the residuals.

figure(fignum)

plot(1:1:length(residual),log10(residual),'blue-o','LineWidth',1,...
    'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
axis([1,length(residual),-inf,inf])
xlabel('iteration','FontSize',16)
ylabel('log_{10} residual','FontSize',16)
title('Residual','FontSize',14)

% Increment the figure number.
fignum = fignum + 1;

%%%%%%%%%%%%%%%%%%%%%%%

% Graph 2: The graph of alpha

figure(fignum);
        
xalpha=1:1:length(alphas);
plot(xalpha,alphas,'-bs', 'LineWidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','MarkerSize',10); 
axis([1,length(alphas),-inf,inf])
xlabel('k','FontSize',16)
ylabel('\alpha','FontSize',16)
title('alpha ','FontSize',14);   
      
% Increment the figure number.
fignum = fignum+1;
        
%%%%%%%%%%%%%%%%%%%%%%%

% Graph 3: The graph of theta.

figure(fignum);

xtheta=1:1:length(thetas);
plot(xtheta,thetas,'-bs', 'LineWidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','MarkerSize',10); 
axis([1,length(thetas),-inf,inf])
xlabel('k','FontSize',16)
ylabel('\theta','FontSize',16)
title('theta ','FontSize',14);   

% Increment the figure number.
fignum = fignum+1;


end


    


