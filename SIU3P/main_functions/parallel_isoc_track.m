function [ISOCHRONS,Basement] = parallel_isoc_track(GCOORD,ELEM2NODE, ...
    Point_id,ISOCHRONS,DISPL,dt,istep,iso_sav,Topography,New_topo, ...
    save_step,save_choice,directory,name)

ISOCHRONS(1:2,:) = track_points(GCOORD,ELEM2NODE,Point_id, ...
    ISOCHRONS(1:2,:),DISPL,dt);
Basement = ISOCHRONS(1:2,ISOCHRONS(3,:)==0);

if floor(istep/iso_sav)==(istep/iso_sav)
    ISOCHRONS = [ISOCHRONS [New_topo; ...
        istep*ones(1,size(Topography,2)); New_topo(2,:)]];
end

% Save ISOCHRONS of the previous time step in their correspondent step
ss.istep = Inf;
c = 0;
while ss.istep ~= istep
    switch save_choice
        case 'cont'
            try
                ss = load([[directory,'/',name],num2str(istep)], ...
                    'istep');
            end
            if ss.istep==istep
                s_name = [[directory,'/',name],num2str(istep)];
                save(s_name,'ISOCHRONS','Basement')
            end
        case 'dis'
            if (istep) == save_step((istep)== ...
                    save_step)
                try
                    ss = load([[directory,'/',name],num2str(istep)], ...
                        'istep');
                end
                if ss.istep==istep
                    s_name = [[directory,'/',name], ...
                        num2str(istep)];
                    save(s_name,'ISOCHRONS','Basement')
                end
            else
                try
                    ss = load([directory,'/',name],'istep');
                end
                if ss.istep==istep
                    save([directory,'/',name],'ISOCHRONS', ...
                        'Basement')
                end
            end
    end
    if ss.istep~=istep
        pause(5)
    end
    c = c+1;
    disp(c)
end