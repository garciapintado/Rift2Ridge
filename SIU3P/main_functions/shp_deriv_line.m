function [N, dNdu] = shp_deriv_line(ipx, nnodel)
    % function [N, dNdu] = shp_deriv_line(ipx, nnodel)
    % +++ purpose +++
    % 1D Canonical Lagrangian shape functions and their derivatives with respect to local coordinates
    % based on known bibliographical values
    % 
    % Author: Javier García_Pintado, MARUM, 2020

    nip  = length(ipx);
    N    = cell(nip,1);
    dNdu = cell(nip,1);
    
    for i=1:nip
        eta  = ipx(i);
        
        switch nnodel
        case 2
            N{i}  = [1/2*(1 - eta);
                     1/2*(1 + eta)];
            dNdu{i} = [-1/2 1/2]';
        case 3  % cuadratic Lagrange polynomials
            N{i} = [1/2*eta*(eta-1);                                          % eta = -1:0.1:1; 
                    1-eta^2;                                                  % figure(); plot(eta,eta.*(eta-1)/2,'-')
                    1/2*eta*(eta+1)];                                         % hold on;  plot(eta,eta.*(eta+1)/2,'-')
            dNdu{i} = [eta-1/2 -2*eta eta+1/2]';
            % hold on;  plot(eta,1 - eta.^2,'-') 
            % for ip=1:nip; xline(ipx(ip)); yline(ipx(ip)); end
        end
    end
end % function

