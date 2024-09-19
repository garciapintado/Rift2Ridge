function plot_defgrid(TRACKP,E2N_TP,tp_x,tp_y,km,type_p,color)
% plot_defgrid(TRACKP,tp_x,tp_y,km)

switch type_p
    case 'elements'
        for j = 1:size(E2N_TP,2)
            plot(TRACKP(1,E2N_TP([1 2 3 4 1],j))/km,TRACKP(2,E2N_TP([1 2 3 4 1],j))/km,color)
            hold on
        end
    case 'layers'
        Layers = 1:2:length(tp_y);
        Elayers = repmat(Layers*(length(tp_x)-1)-(length(tp_x)-1), ...
            length(tp_x)-1,1) + ...
            repmat((1:length(tp_x)-1)',1,size(Layers,2));
        Elayers = Elayers(:);
        for j = 1:length(Elayers)
            patch(TRACKP(1,E2N_TP(:,Elayers(j)))'/km,TRACKP(2,E2N_TP(:,Elayers(j)))'/km,color)
        end
end

hold off