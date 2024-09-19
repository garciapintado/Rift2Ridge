function [layer1x, layer1y, layer2x, layer2y, remesh, layer_corr] = lineintersect (layer1x, layer1y, layer2x, layer2y, remesh, adjust)
% TODO This function is a pile of crap. Improve it using interp1

% figure(22)
% plot(layer1x,layer1y,layer2x,layer2y)

layer_corr = [0 0];

for i = 1:length(layer1x)
    % Find the two neighbourg points in x1
    
    diffx = layer1x(i)-layer2x;
    diffxneg = diffx<0;
    if sum(diffxneg) == 0;
        continue
    end
    postindex = find(min(layer2x(diffxneg))==layer2x);
    diffxpos = diffx>=0;
    if sum(diffxpos) == 0;
        continue
    end
    previndex = find(max(layer2x(diffxpos))==layer2x);
    
    % Interpolates y1 at the x2 point
    layer12y = (layer1x(i)-layer2x(previndex))*((layer2y(postindex)-layer2y(previndex))/(layer2x(postindex)-layer2x(previndex)))+layer2y(previndex);
    
    % Check if 2 is below 1 and adjust in that case
    if layer1y(i) >= layer12y - adjust
        layer1y(i) = layer12y - adjust;
        remesh = 1;
        layer_corr(1) = 1;
    end
end

% figure(12)
% plot(layer1x,layer1y,'o',layer2x,layer2y,'o')
% hold on
% plot(layer1x,layer1y,layer2x,layer2y)
% hold off

for i = 1:length(layer2x)
    % Find the two neighbourg points in x1
    diffx = layer2x(i)-layer1x;
    diffxneg = diffx<0;
    if sum(diffxneg) == 0;
        continue
    end
    postindex = find(min(layer1x(diffxneg))==layer1x);
    diffxpos = diffx>=0;
    if sum(diffxpos) == 0;
        continue
    end
    previndex = find(max(layer1x(diffxpos))==layer1x);
    
    % Interpolates y2 at the x1 point
    layer21y = (layer2x(i)-layer1x(previndex))*((layer1y(postindex)-layer1y(previndex))/(layer1x(postindex)-layer1x(previndex)))+layer1y(previndex);
    
    % Check if 1 is above 2 and adjust in that case
    if layer2y(i) - adjust <= layer21y
        layer1y(previndex) = layer1y(previndex) + layer2y(i) - layer21y - adjust;
        layer1y(postindex) = layer1y(postindex) + layer2y(i) - layer21y - adjust;
        remesh = 1;
        layer_corr(2) = 1;
    end
end

% figure(23)
% plot(layer1x,layer1y,'o',layer2x,layer2y,'o')
% hold on
% plot(layer1x,layer1y,layer2x,layer2y)
% hold off

disp(layer_corr)

end
