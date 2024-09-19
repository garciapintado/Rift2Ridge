% HISTORIC 2ND STRAIN INVARIANT PLOT
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_int     Color of the box        #interfaces x 3 vector  Black
%               and interfaces          with values [0 1]
%
% line_width    Width of the box        #interfaces x 1 vector  1
%               and interfaces          with width values
%
% log_data      Plot logarithmic data   0 no logarithmic data   1
%                                       1 logarithmic data

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
I2_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy      = Nbig'\I2.p(i,:)';
    I2_n(:,i)= Dummy(1:3);
end

% Correct negative E2_n which are the result of interpolation
I2_n = abs(I2_n);

% I2_n = I2(:,1:3)';
% nel = size(ELEM2NODE,2);
% EL2N(1,:) = 1:3;
% for i=1:nel
%     is         = (i-1)*3+1; ie = (i-1)*3+3;
%     GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
%     EL2N(i,:) = is:ie;
% end

%==========================================================================
% PLOT
%==========================================================================
% Plot historic strain
% --------------------
title(['Historic 2nd strain invariant (',num2str(istep*dt/ma-dt/ma),' Myr)'])
xlabel('Distance [km]')
ylabel('Depth [km]')

if exist('log_data','var')
    if log_data
        patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            log10(I2_n(:)),'FaceColor','flat')
    else
        patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            I2_n(:),'FaceColor','flat')
    end
else
    patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
        I2_n(:),'FaceColor','flat')
end
axis tight
shading interp
colorbar
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off