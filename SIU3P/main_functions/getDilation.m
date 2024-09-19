function dilation = getDilation(Phases,w_ml,vs,x_max_dike,x_max_ml,y_top_ml,y_bot_ml,y_moho, GCOORD,ELEM2NODE, area_el,dila)
    %
    % OUTPUT
    % dilation [nel,1]
    area_el = area_el(Phases==5);

    fprintf(1, 'Dilation:      '); tic;

    dila_dikes = vs*abs(y_top_ml)/(x_max_dike*abs(y_top_ml));
    dila_meltl = w_ml     * vs*abs(y_moho - y_top_ml)/(x_max_ml*abs(y_bot_ml - y_top_ml));
    dila_sill  = (1-w_ml) * vs*abs(y_moho - y_top_ml)/(x_max_ml*abs(y_moho - y_bot_ml));

    nel                 = length(Phases);
    dilation            = zeros(nel,1);
    dilation(Phases==4) = -dila_dikes;
    dilation(Phases==3) = -dila_meltl;

    %sheeted sills- linear decrease of dilation with depth
    if dila=="linear"
        z_el    = mean(reshape(GCOORD(2,ELEM2NODE(1:3,Phases==5)),3,[]));

        wdil = (y_bot_ml-z_el)./(y_bot_ml-y_moho);

        % Normierung
        a1 = x_max_ml*(y_bot_ml-y_moho);
        a2 = sum(area_el.*wdil');
        wdil = wdil .* (a1/a2);
        dilation(Phases==5) = -((1-w_ml) * vs*abs(y_moho - y_top_ml))./(a1.*wdil);

    elseif dila=="deep_ml"
        z_el    = mean(reshape(GCOORD(2,ELEM2NODE(1:3,Phases==5)),3,[]));
        ind = z_el < -4000.;
        Phase5 = find(Phases==5);
        dilation(Phase5(ind)) = -dila_sill;
        
    elseif dila==""
        dilation(Phases==5) = -dila_sill;
    else
        error('distribution of sills is not defined!')
    end

    fprintf(1,'top melt lens: %6.4f dila dike: %6.4e dila ml: %6.4e meltlens_calc: %6.4f \n',...
        y_top_ml/1000, dila_dikes, dila_meltl, y_top_ml/1000);
    end
end % function

