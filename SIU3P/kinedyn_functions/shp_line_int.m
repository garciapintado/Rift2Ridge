function  [N] = shp_line_int(ipx, nnodel)

nip  = length(ipx);
N    = cell(nip,1);

for i=1:nip
    eta  = ipx(i);

    switch nnodel
        case 2
            SHP   = [0.5*(1 - eta);
                     0.5*(1 + eta)];

        case 3  % cuadratic Lagrange polynomials
            SHP = [eta*(eta-1)/2;                                          % eta = -1:0.1:1; 
                   1-eta^2;                                                % figure(); plot(eta,eta.*(eta-1)/2,'-')
                   eta*(eta+1)/2];                                         % hold on;  plot(eta,eta.*(eta+1)/2,'-')
                                                                           % hold on;  plot(eta,1 - eta.^2,'-') 
                                                                           % for ip=1:nip; xline(ipx(ip)); end; yline(0)
    end
    
    N{i} = SHP;
end