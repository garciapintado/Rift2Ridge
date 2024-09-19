function c = cmap_sedage(numc,bk_age,pbbk,pabk)
% C = CMAP_SEDAGE(NUMC,PBBK,PABK,BK_AGE) constructs a color scale for
% sediment ages. NUMC is the number of intervals for colors to happen
% (color resolution). BK_AGE is the age of breakup, if not needed make it
% empty []. PBBK are the color points previous to breakup and PABK are the
% color points after breakup. This color points are linearly distributed
% along the time at equidistantly.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, postdoc at University of
% Bremen, 05-11-2018. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Get value limits
val_lim = get(gca,'CLim');
val_lim = [floor(val_lim(1)) ceil(val_lim(2))];
c = get(gcf,'colormap')';

if nargin<2
    bk_age = val_lim(2);
end

if nargin<3
    % Color points to breakup
    pbbk = [0.8 0 0.8; 0.49 0.18 0.56; ...
        0 0 1; ...
        0 1 1;
        1 1 0; 1 0.5 0;
        0.5 0.3 0];
    
    % Color points after breakup
    pabk = [1 1 1; ...
        0.7 0 0];
end

% Calculate colors for sediments pre-breakup sediments
bkp = bk_age*numc/diff(val_lim);
xbbk = linspace(1,bkp,size(pbbk,1));
rbbk = interp1(xbbk,pbbk(:,1),1:floor(bkp),'linear');
gbbk = interp1(xbbk,pbbk(:,2),1:floor(bkp),'linear');
bbbk = interp1(xbbk,pbbk(:,3),1:floor(bkp),'linear');
if floor(bkp)==bkp
    rbbk(end) = pbbk(end,1);
    gbbk(end) = pbbk(end,2);
    bbbk(end) = pbbk(end,3);
end

% Calculate colors for sediments post-breakup sediments
xabk = linspace(floor(bkp)+1,numc,size(pabk,1));
rabk = interp1(xabk,pabk(:,1),floor(bkp)+1:numc,'linear');
gabk = interp1(xabk,pabk(:,2),floor(bkp)+1:numc,'linear');
babk = interp1(xabk,pabk(:,3),floor(bkp)+1:numc,'linear');

% Build up colors
c = [[rbbk rabk]' [gbbk gabk]' [bbbk babk]'];