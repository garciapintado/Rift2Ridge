% CONVERGENCE PLOT

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, 21-07-2015. Email: 
% mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Input
steps_p = 1:length(nw_it);

% Plot iterations per step
plot(steps_p,nw_it,'-ob')
title('Convergence')
xlabel('Steps')
ylabel('Number of iterations before convergence')
hold on

% Plot remeshing
plot(steps_p(track_remesh==1),nw_it(track_remesh==1),'o', ...
    'MarkerEdgeColor','none','MarkerFaceColor','r')
legend('Iterations','Remesh steps')
hold off