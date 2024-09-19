% CRUSTAL THICKNESS ALONG TIME

%--------------------------------------------------------------------------
% Function written by Elena Ros
%--------------------------------------------------------------------------

time =[dt/ma:dt/ma:dt*istep/ma];
plot(time(1:istep),Crust_thickness(1:istep)/1000,'o','MarkerSize',3) %,'linewidth',30)
xlabel('Time [Ma]','FontSize',9)
ylabel('Crustal thickness [km]','FontSize',9)
title(['Magmatic crustal thickness ',num2str((istep*dt/ma)),' Ma'])
grid on