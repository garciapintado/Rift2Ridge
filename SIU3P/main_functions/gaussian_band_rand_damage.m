function RDWS_factorIP = gaussian_band_rand_damage(RDWS, GCOORD, ELEM2NODE, nip)
%
% Get the x coordinates of the central node
% return: RDWS_factorIP cell array, containing two [nel,nip] arrays

X_NODE7 = GCOORD(1,ELEM2NODE(7,:))';         % [nel,1] central nodes                                      
xc = mean(minmax(GCOORD(1,:)));              % x-center of domain

for n = 1:length(RDWS.factor)
    % Calculate the factors for this coordinates using a gaussian function
    RDWS_factor = 1. + RDWS.factor(n) * exp(-(X_NODE7-xc).^2 ./ ((RDWS.sigma/2)^2));       % [] 
    RDWS_factorIP{n} = repmat(RDWS_factor,1,nip);                                          %

    %     hold on
    %     plot(X_NODE7/1000,RDWS_factorIP{n},'.')
    %     plot([-RDWS.sigma RDWS.sigma; -RDWS.sigma RDWS.sigma]/2/1000, ...
    %         [1 1; RDWS.factor(n) RDWS.factor(n)])
    %     hold off
end
