function [MESH]=Te_all(MESH,VAR,PHYSICS,NUMSCALE)

[MESH] = curvature(MESH);

space=0.05; %resolution of new mesh in km over the depth

inod_surface = find(ismember(MESH.PointID,MESH.PointID_top));
inod_bdt=find(ismember(MESH.PointID,MESH.PointID_bdt));
inod_moho=find(ismember(MESH.PointID,MESH.PointID_moho));
inod_lab=find(ismember(MESH.PointID,MESH.PointID_lab));

x_surface=MESH.GCOORD(1,inod_surface);
[x_surface,order_x]=sort(x_surface);
inod_surface = inod_surface(order_x);

Te=length(x_surface);

for i=1:length(x_surface)
    
    % loop over all surface points to estimate Te
    top_z=MESH.GCOORD(2,inod_surface(i));
    [dif1 index_lab] = min( abs(MESH.GCOORD(1,inod_lab) - MESH.GCOORD(1,inod_surface(i))));
    bot_z=MESH.GCOORD(2,inod_lab(index_lab));
    
    x_points=x_surface(i);
    cur_space= (top_z-bot_z)/(floor((top_z-bot_z)/space)+1);
    z_points=bot_z:cur_space:top_z;
   
    [x,z]=meshgrid(x_points,z_points);
    new_MESH=[x(:),z(:)]';
    
    [~,gX_PT] = find_pts_outside_mesh(MESH,new_MESH);
    [els,lc]  = locate_points_2d(MESH.GCOORD,MESH.EL2NOD,gX_PT);

    T_profile    = interp2d_tri367(MESH.EL2NOD(1:6,:),els,lc,VAR.T);
    T_profile=T_profile+273;
    depth_profile = new_MESH(2,:);
    
    

    
%     Er_II_profile=interp2d_tri367(MESH.EL2NOD(1:3,:),els,lc,VAR.Er_II_totl);
%     Er_II_profile=Er_II_profile./NUMSCALE.t0;  
%     rate=mean(Er_II_profile);

    % by now new mesh (mesh for colomn) with resolution close to space is
    % created. T_profile and depth_profile contain mentioned values under
    % the required surface point
    % arrays rev_depth_profile and rev_T_profile contaian the same data but
    % ordered from surface to depth
    
    [rev_depth_profile,order_depth]=sort(depth_profile,'descend');
    rev_T_profile =T_profile(order_depth);
%    rev_Er_II_profile=Er_II_profile(order_depth);
    
    [dif2 index_bdt] = min(abs(MESH.GCOORD(1,inod_bdt) - MESH.GCOORD(1,inod_surface(i))));
    [dif3 nnc] = min(abs(rev_depth_profile - MESH.GCOORD(2,inod_bdt(index_bdt))));

    [dif4 index_moho] = min(abs(MESH.GCOORD(1,inod_moho) - MESH.GCOORD(1,inod_surface(i))));
    [dif5 nncl] = min(abs(rev_depth_profile - MESH.GCOORD(2,inod_moho(index_moho))));
     
    %here we have prepared arrays rev_depth_profile (values mostly negative
    %and ordered to the depth), rev_T_profile (temperature profile
    %corresponding to depth profile) and variables nnc and nncl which are
    %numbers of elements corresponding to bdt and moho respectively. From
    %here we can call Te_Lowry
    
    rev_Litho=[];
    rev_Dens=[];
    for j=1:length(rev_depth_profile)
        if(j<nnc) 
            phaseid=3;
        elseif(j>=nnc & j<nncl)
            phaseid=2;
        else
            phaseid=1;
        end
        
        rev_Dens=[rev_Dens PHYSICS.Dens(phaseid)*(1-PHYSICS.alpha(phaseid)*rev_T_profile(j))];
        
        if(j>1)
            rev_Litho=[rev_Litho rev_Litho(j-1)-rev_Dens(j)*PHYSICS.g*(rev_depth_profile(j-1)-rev_depth_profile(j))];
        else
            rev_Litho=0;
        end
    end
    rev_Litho=rev_Litho.*1.0e3;
    
    rev_Dens=rev_Dens';
    
    rho_uc=mean(rev_Dens(1:nnc));
    rho_lc=mean(rev_Dens((nnc+1):nncl));
    rho_m=mean(rev_Dens((nncl+1):end));


    INPUT_COL.nnc=nnc;
    INPUT_COL.nncl=nncl;
    INPUT_COL.rev_T_profile=rev_T_profile;
    INPUT_COL.space=cur_space;
    INPUT_COL.rev_depth_profile=rev_depth_profile;
    %INPUT_COL.rev_Er_II_profile=rev_Er_II_profile;
    INPUT_COL.rev_Litho=rev_Litho;
    INPUT_COL.curv=MESH.surf_curv(i);
    INPUT_COL.tc_up=nnc*cur_space;
    INPUT_COL.tc_low=nncl*cur_space;
    INPUT_COL.tm=-rev_depth_profile(end);
%    INPUT_COL.rate=rate;

    
    close(figure(666))
    disp(x_surface(i))
    Te(i)=Te_Lowry(MESH,VAR,PHYSICS,INPUT_COL);
%     Te(:)=10; 
    
    
end

while (any(isnan(Te)))

    for i=1:length(x_surface)
        if(isnan(Te(i)))
            if((i-1)>0)
                if(~isnan(Te(i-1)))
                    left=Te(i-1);
                else
                    left=0;
                end
            else
                left=0;
            end
            
            if((i+1)<=length(x_surface))
                if(~isnan(Te(i+1)))
                    right=Te(i+1);
                else
                    right=0;
                end
            else
                right=0;
            end
            
            if(left*right)
                if((x_surface(i)-x_surface(i-1)) < (x_surface(i+1)-x_surface(i)))
                    Te(i)=Te(i-1);
                else
                    Te(i)=Te(i+1);
                end
            elseif(left || right)
                Te(i)=left+right;
            end
                
        end
    end
    
end







if isfield(MESH,'Te') 
    clearvars MESH.Te;
end

MESH.Te=Te;
    
    
    
end









