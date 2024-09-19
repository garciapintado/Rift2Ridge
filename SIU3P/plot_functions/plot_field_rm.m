function plot_field_rm(XX,YY,FIELD,n,clim,logc,Topography)
% PLOT_EH_RM(XX,YY,FIELD,N,LIM,LOGC,TOPOGRAPHY) plots a FIELD defined at a
% regular mesh XX, YY using matlab function contourf and cuts the
% topography unsing the TOPOGRAPHY. N defines the number of different
% colors used for the countour fills. CLIM defines the upper and lower
% limits on the colormap (leave it empty for default limits). Set LOGC 1
% for plotting the field in logarithmic (base 10) scale and and LOGC 0 for
% not logarithmic scale.

if logc
    FIELDp = log10(FIELD);
else
    FIELDp = FIELD;
end

% Plot contours
contourf(XX,YY,FIELDp,n,'LineStyle','none')
if ~isempty(clim)
    caxis([clim(1) clim(2)])
end
colorbar
hold on

% Remove areas outside of the domain
patch([Topography(1,:) Topography(1,end) Topography(1,1) Topography(1,1)], ...
    [Topography(2,:) 20 20 Topography(2,1)],'w','LineStyle', ...
    'none')