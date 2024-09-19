% WEAK SEED PLOT

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% % Calculate points of the weak seed
% WSp(1,:) = WS.size*cosd([0:360])+WS.coord(1);
% WSp(2,:) = WS.size*sind([0:360])+WS.coord(2);
hold on
plot(WS.coord(:,1)/1000,WS.coord(:,2)/1000,'.k')
% plot(WSp(1,:)/1000,WSp(2,:)/1000,'k')