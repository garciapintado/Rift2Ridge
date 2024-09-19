function var2tec(dirp,diro,VAR) 
% VAR2TEC(DIRP,DIRO,VAR) load the saved time steps from a model in the
% folder DIRP and produce a .plt file for tecplot to visualise the results.
% The .plt file is saved in DIRO under the name of the model folder.
% Variables stored are defined in the cell VAR. This function uses
% mat2tecplot.m function downloaded from tecplot website. The supported
% variables are:
% - E2all:      strain rates
% - Temp:       temperatures

% Author: Miguel Andres-Martinez, University of Bremen
%         andresma@uni-bremen.de

%==========================================================================
%% INPUT
%==========================================================================
% Number of variables to plot
nvar    = length(VAR);

% File names for every step
stepsf  = dir(dirp);
stepsf  = stepsf(3:end);

% Load first file names and saving options
load([dirp,'/',stepsf(1).name],'save_choice','name','sav_every_step', ...
    'surf_proc_s','Point_id')

% Name plt file
bars    = strfind(dirp,'/');
if bars(end)==length(dirp)
    nameplt = dirp(bars(end-1)+1:bars(end)-1);
else
    nameplt = dirp(bars(end)+1:end);
end

% Initialise variables for mat2tecplot
tdata.vformat   = 2*ones(1,nvar+3);
tdata.Nvar      = 3+nvar;
tdata.varnames  = {'Distance','Depth','z'};

% Box plotting
interf = [1 3:3:max(Point_id)-1];
for o = 1:length(interf)
    tdatal(o).title = ['Int',interf(o)];
    tdatal(o).Nvar = 2;
    tdatal(o).varnames = {'x','y'};
end
nint = length(interf);
if surf_proc_s
    tdatal(nint+1).title = 'TopBasement';
    tdatal(nint+1).Nvar = 2;
    tdatal(nint+1).varnames = {'x','y'};
end

c               = 3;
for n = 1:nvar
    if strcmp(VAR{n},'Vel')
        c                   = c+1;
        tdata.varnames{c}   = 'Vx';
        c                   = c+1;
        tdata.varnames{c}   = 'Vy';
        tdata.Nvar          = 4+nvar;
        tdata.vformat       = 2*ones(1,nvar+4);
    else
        c                   = c+1;
        tdata.varnames{c}   = VAR{n};
    end
end

% Build names of files and time steps
nl      = length(name)+1;
step    = [];
laststr = {};
c       = 0;
for n = 1:length(stepsf)
    if ~strcmp(stepsf(n).name,[name,'.mat'])
        if isnan(str2double(stepsf(n).name(2)))
%             step(c)     = str2double(stepsf(n).name(nl:end-5));
%             laststr{c}  = stepsf(n).name(end-5);
        else
            c   = c+1;
            step(c)     = str2double(stepsf(n).name(nl:end-4));
            laststr{c} = [];
        end
        
    end
end
% Sort time steps
[step,stepsort] = sort(step);

% Check for structures
cp = 0;
for n = 1:nvar
    pointn      = strfind(VAR{n},'.');
    if ~isempty(pointn)
        VAR_load{n+cp} = VAR{n}(1:pointn-1);
    elseif strcmp(VAR{n},'ErP')
        VAR_load = {VAR_load{:} 'Gamma' 'YC' 'Yield_T2' 'TAU_xx' 'TAU_xy' ...
            'TAU_yy'};
        cp = cp+5;
    elseif strcmp(VAR{n},'ErC')
        VAR_load = {VAR_load{:} 'TAU_xx' 'TAU_xy' 'TAU_yy' ...
            'Mu_dif_all' 'Mu_dis_all'};
        cp = cp+4;
    else
        VAR_load{n+cp} = VAR{n};
    end
end

% Add default variables for plotting
addvar  = {'GCOORD','ELEM2NODE','nel','dt','istep','ma','GEOMETRY', ...
    'Geo_id','Point_id','Corner_id','Cornin_id','nip'};
nav     = length(addvar);

for n = 1:nav
    VAR_load{nvar+cp+n}    = addvar{n};
end

%==========================================================================
%% STEPS LOOP
%==========================================================================
cc  = 1;
fprintf('MAPPING DATA TO NODES:       0')
for n = 1:length(step)
    % Display percentage
    p   = num2str(floor(n/length(step)*100));
    for pn = 1:cc
        fprintf('\b')
    end
    fprintf(p)
    cc = length(p);
    
    % Load step
    L = load([dirp,'/',name,num2str(step(n)),laststr{stepsort(n)}],VAR_load{:});
    
    e2n     = reshape(1:L.nel*3,3,L.nel)';
    gc      = L.GCOORD(:,L.ELEM2NODE([1 2 3],:));
    
    % Define 
    tdata.FEsurfaces(n).zonename    = 'my zone';
    tdata.FEsurfaces(n).x           = gc(1,:);
    tdata.FEsurfaces(n).y           = gc(2,:);
    tdata.FEsurfaces(n).order       = 3;
    tdata.FEsurfaces(n).e2n         = e2n;
    tdata.FEsurfaces(n).solutiontime= L.dt*L.istep/L.ma;
    
    % Box plotting
    for o = 1:length(interf)
        tdatal(o).lines(n).solutiontime = L.dt*L.istep/L.ma;
    end
    if surf_proc_s
        tdatal(nint+1).lines(n).solutiontime = L.dt*L.istep/L.ma;
    end
    
    %======================================================================
    %% VARIABLE LOOP
    %======================================================================
    c = 0;
    for m = 1:nvar
        % Variable name and values
        vname   = VAR{m};
        pointn  = strfind(vname,'.');
        if ~isempty(pointn)
            vval    = L.(VAR{m}(1:pointn-1)).(VAR{m}(pointn+1:end));
        elseif strcmp(VAR{m},'ErP') || strcmp(VAR{m},'ErC')
            % Do nothing
        else
            vval    = L.(VAR{m});
        end
        
        % Variable switch
        double_var = 0;
        switch vname
            case {'Temp'}
                var3    = vval(L.ELEM2NODE([1 2 3],:));
                var3    = var3(:)';
            case {'Vel'}
                Vx      = vval(1:2:end-1);
                Vy      = vval(2:2:end);
                var3    = Vx(L.ELEM2NODE([1 2 3],:));
                tdata.FEsurfaces(n).v(m+c,:) = var3(:)';
                c       = c+1;
                var3    = Vy(L.ELEM2NODE([1 2 3],:));
                tdata.FEsurfaces(n).v(m+c,:) = var3(:)';
                double_var = 1;
            case {'E2all','Mu_all','I2.f','I2.p','I2.c','RHO'}
                var3 = ip2nodes(vval,L.GCOORD,L.ELEM2NODE,size(vval,2),3);
                var3 = var3(:)';
            case {'ErP'}
                Er_xx   = zeros(size(L.ELEM2NODE,2),L.nip);
                Er_xy   = zeros(size(L.ELEM2NODE,2),L.nip);
                Er_yy   = zeros(size(L.ELEM2NODE,2),L.nip);
                
                Er_xx(L.YC~=0)  = 0.5 * L.Gamma(L.YC~=0)...
                    .*L.TAU_xx(L.YC~= 0)./L.Yield_T2(L.YC~= 0);
                Er_xy(L.YC~=0)  = 0.5 * L.Gamma(L.YC~=0)...
                    .*L.TAU_xy(L.YC~=0)./L.Yield_T2(L.YC~= 0);
                Er_yy(L.YC~=0)  = 0.5 * L.Gamma(L.YC~= 0)...
                    .*L.TAU_yy(L.YC~=0)./L.Yield_T2(L.YC~=0);
                
                ErP     = sqrt(0.5*(Er_xx.^2 + Er_yy.^2) + Er_xy.^2);
                
                var3 = ip2nodes(ErP,L.GCOORD,L.ELEM2NODE,size(vval,2),3);
                var3 = var3(:)';
            case {'ErC'}
                mdin = isnan(L.Mu_dif_all);
                iMu_c = zeros(size(L.ELEM2NODE,2),L.nip);
                iMu_c(~mdin) = (1./L.Mu_dis_all(~mdin) + ...
                    1./L.Mu_dif_all(~mdin));
                iMu_c(mdin) = 1./L.Mu_dis_all(mdin);
                ErC = 0.5.*sqrt(0.5*(L.TAU_xx.^2+L.TAU_yy.^2)+L.TAU_xy.^2) ...
                    .*iMu_c;
                
                var3 = ip2nodes(ErC,L.GCOORD,L.ELEM2NODE,size(vval,2),3);
                var3 = var3(:)';
            otherwise
                error(['Variable "',vname,'" not recognised'])
        end
        
        if ~double_var
            % Set values to variables
            tdata.FEsurfaces(n).v(m+c,:)=var3;
        end
    end
    
    % BOX AND TOP BASEMENT
    L.Point_id(L.Corner_id) = [1 1 max(L.Point_id)-1 max(L.Point_id)-1];
    L.Point_id(L.Cornin_id) = reshape([interf(2:end-1); interf(2:end-1)], ...
        1,length(interf(2:end-1))*2);
    node_edge = ismember(1:size(L.GCOORD,2)-size(L.ELEM2NODE,2), ...
        L.ELEM2NODE(1:3,:));
    for o = 1:nint
        indx_int_edge = find(L.Point_id==interf(o) & node_edge);
        indx_e2n = ismember(L.ELEM2NODE,indx_int_edge);
        indx_2nodes = repmat(sum(indx_e2n)==2,7,1);
        indx_e2n = indx_e2n & indx_2nodes;
        e2ng = reshape(L.ELEM2NODE(indx_e2n),2,sum(indx_e2n(:))/2);
        [~,e2n] = ismember(e2ng,indx_int_edge);
        e2n = unique(sort(e2n)','rows')';
        coord_int = L.GCOORD(:,indx_int_edge);
        cont_int_coord = seg_connect(coord_int,e2n);
        tdatal(o).lines(n).x = cont_int_coord(1,:);
        tdatal(o).lines(n).y = cont_int_coord(2,:);
    end
    if surf_proc_s
        load([dirp,'/',name,num2str(step(n)),laststr{stepsort(n)}], ...
            'Basement');
        tdatal(nint+1).lines(n).x = Basement(1,:);
        tdatal(nint+1).lines(n).y = Basement(2,:);
    end
end

fprintf('\n')
mkdir([diro,nameplt])
% 2D data writting
outp(1) = mat2tecplot(tdata,[diro,nameplt,'/fields.plt']);

% Lines data writting
for o = 1:length(tdatal)
    outp(1+o) = mat2tecplot(tdatal(o),[diro,nameplt,'/lines',num2str(o),'.plt']);
end

if any(outp==-1)
    error('Error in producing the .plt file')
else
    fprintf('COMPLETED')
    fprintf('\n')
end