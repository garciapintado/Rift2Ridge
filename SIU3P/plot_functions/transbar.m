function axt = transbar(color,taxis,raxis,position,orient)
% TRANSPARENT BAR

switch orient
    case 'vertical'
        axt = axes('Color','none','XColor','none','Position',position, ...
            'YAxisLocation','right');
    case 'horizontal'
        axt = axes('Color','none','YColor','none','Position',position);
end

ry = taxis(1):diff(taxis)/10:taxis(2);
sy = length(ry);

X = [zeros(size(ry)) ones(size(ry))]*taxis(2)*0.1;
Y = [ry ry];

C = [[1:sy-1; (1:sy-1)+sy; (2:sy)+sy] [1:sy-1; 2:sy; (2:sy)+sy]];

switch orient
    case 'vertical'
        patch('faces',C','vertices',[X' Y'],'FaceVertexAlphaData', ...
            Y(:),'AlphaDataMapping','scaled','FaceAlpha','interp',...
            'FaceColor',color,'EdgeColor','none')
        hold on
        plot([0 0 1 1 0]*taxis(2)*0.1,[ry(1) ry(end) ry(end) ry(1) ry(1)],'k')
        axt.ALim = raxis;
        axt.YLim = taxis;
        hold off
    case 'horizontal'
        patch('faces',C','vertices',[Y' X'],'FaceVertexAlphaData', ...
            Y(:),'AlphaDataMapping','scaled','FaceAlpha','interp',...
            'FaceColor',color,'EdgeColor','none')
        hold on
        plot([ry(1) ry(end) ry(end) ry(1) ry(1)],[0 0 1 1 0]*taxis(2)*0.1,'k')
        axt.ALim = raxis;
        axt.XLim = taxis;
        hold off
end

axis equal
axis tight