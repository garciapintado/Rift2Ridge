function [TopoXY,SP] = diff_topo_sealand(TopoXY,GCOORD,SP,dt,ma)
% VERSION OF diff_topo_sealand TO LIMIT SLOPES UNDER THE SEA
% (NOT WORKING!)

% [TOPOXY,SP] = DIFF_TOPO_ARMITAGE(TOPOXY,GCOORD,SP,DT,MA) takes current
% topography TOPOXY and calculates new topographies by diffusion taking
% into account the parameters defined by GCOORD, SP, DT and MA. This
% algorithm uses the finite element diffusion function
% function_nl_erosion_fixedelevation developed by John Armitage.
% ALL units must be in International System.

% TODO test if uplift and pelagic terms should be source terms in the
% diffusion equation.
% TODO due to the adimensionalisation, it is not possible to have no
% hill-slope diffusion. Fix this in the future.

cwide       = 500; % [m]
sres        = 10; % [m]
highan      = 20;
highan_air  = 25;

%==========================================================================
%% SCALING
%==========================================================================
% Longitude of the model [m]
model_long  = max(GCOORD(1,:))-min(GCOORD(1,:));

% Make variables adimensional
if SP.scale
    kappa = SP.kappa/SP.kappa;
    De = SP.c*SP.alpha_sed^SP.nexp*model_long/SP.kappa;
    sealevel = SP.sealevel/model_long;
    Dx_sed = diff(TopoXY(1,:));
    Dx_sed = Dx_sed'/model_long;
    dt_sed = SP.dt*SP.kappa/(model_long*model_long);
    Topography = TopoXY(2,:)';
    Topography = Topography/model_long;
    ksea = SP.ksea*ma/(1e6*SP.kappa);
    kappa_s = SP.kappa_s*ma/(1e6*SP.kappa);
    kdecay = SP.kdecay*model_long;
else
    kappa = SP.kappa;
    De = SP.c*SP.alpha_sed^SP.nexp;
    sealevel = SP.sealevel;
    dt_sed = SP.dt;
    ksea = SP.ksea;
    kappa_s = SP.kappa_s;
    kdecay = SP.kdecay;
end

Old_topo = TopoXY;
MM_SED = [1/3 1/6; 1/6 1/3]; %1/9*([16/5 4/5; 4/5 16/5]+[5/8 5/8; 5/8 5/8]); %[1/3 1/6; 1/6 1/3]; % [1/4 1/4; 1/4 1/4];

% Number of sediment time steps to reach a mechanical time step
ntsteps_sed = dt/SP.dt;
% Calculates the size of the last time step
last_dts = SP.dt*(dt/SP.dt-floor(dt/SP.dt));
if last_dts~=0
    ntsteps_sed = ntsteps_sed+1;
end
% Loop through the erosion/sedimentation time steps
for j = 1:ntsteps_sed
    % Remesh steep slopes in the sea
    %-------------------------------
    if ksea~=0
        % Find submarine nodes
        subm_nodes = TopoXY(2,:)<SP.sealevel;
        % Find the steep slopes
        steep_slopes = ...
            atand(abs(diff(TopoXY(2,:))./diff(TopoXY(1,:))))>=highan;
        steep_slopes = ([steep_slopes 0]+[0 steep_slopes])~=0;
        % Find subaerial steep slopes
        steep_slopes_air = ...
            atand(abs(diff(TopoXY(2,:))./diff(TopoXY(1,:))))>=highan_air;
        steep_slopes_air = ([steep_slopes_air 0]+[0 steep_slopes_air])~=0;
        % Find shallow water and near sea subaerial topographies
        near_sea = TopoXY(2,:)>=sealevel-cwide/2 & ...
            TopoXY(2,:)<=sealevel+cwide/2;
        % Resolution check
        need_remesh = diff(TopoXY(1,:))>2*sres;
        need_remesh = ([need_remesh 0]+[0 need_remesh])~=0;
        % Submarine steep slopes
        sbm_ss = ...
            ((subm_nodes & steep_slopes) | near_sea | steep_slopes_air) ...
            & need_remesh;
        % Calculate steep slope segment indexes
        counter_sbm_ss = 1:length(sbm_ss);
        sbm_ss_seg = sbm_ss*1;
        sbm_ss_seg(sbm_ss) = ...
            [1 cumsum(diff(counter_sbm_ss(sbm_ss))>1)+1];
        % Number of segments
        nseg = unique(sbm_ss_seg(sbm_ss_seg~=0));
        
%         % Plot
%         figure(2); clf
%         plot([TopoXY(1,1) TopoXY(1,end)],[SP.sealevel SP.sealevel],'--')
%         hold on
%         plot(TopoXY(1,:),TopoXY(2,:))
%         plot(TopoXY(1,subm_nodes),TopoXY(2,subm_nodes),'x')
        % Loop through steep slope submarine segments
        for n = 1:size(nseg,2)
            % Indexes of the nodes in the current segment
            seg_nodes = find(sbm_ss_seg==n);
            xleft   = TopoXY(1,seg_nodes(1));
            xright  = TopoXY(1,seg_nodes(end));
            % Remesh x
            xslope = linspace(xleft,xright,round((xright-xleft)/sres));
            % Remesh y
            yslope = interp1(TopoXY(1,seg_nodes),TopoXY(2,seg_nodes),xslope);
            % Add remeshed nodes to topography
            ileft = 1:length(sbm_ss_seg)<seg_nodes(1);
            iright = 1:length(sbm_ss_seg)>seg_nodes(end);
            TopoXY = [TopoXY(:,ileft) [xslope; yslope] TopoXY(:,iright)];
            sbm_ss_seg = [sbm_ss_seg(ileft) (n*ones(length(xslope),1))' ...
                sbm_ss_seg(iright)];
        end
%         plot(TopoXY(1,sbm_ss_seg~=0),TopoXY(2,sbm_ss_seg~=0),'o')
    end
    
    Dx_sed = diff(TopoXY(1,:))';
    Topography = TopoXY(2,:)';
    no_nodes_sed = size(Topography,1);
    Rhs_sed = zeros(size(Topography));
    
    % If it is the last sediment time step and the last_dts is
    % different from 0
    if SP.scale
        if j==ntsteps_sed && last_dts~=0
            dt_sed = last_dts*SP.kappa/(model_long*model_long);
        end
        
        dt_sed_dim = dt_sed*(model_long*model_long)/(SP.kappa);
    else
        if j==ntsteps_sed && last_dts~=0
            dt_sed = last_dts;
        end
        
        dt_sed_dim = dt_sed;
    end
    
    % Pelagic
    tpelagic = pelagic(SP,TopoXY,dt_sed_dim);
    if SP.scale
        tpelagic = tpelagic/model_long;
    end
    Topography = Topography + tpelagic';
    % Calculating new topography
    Topography = function_sealand_k_slim ...
            (Topography,no_nodes_sed,Rhs_sed,MM_SED,De,SP.nexp,dt_sed, ...
            kappa,Dx_sed,sealevel,ksea,kappa_s,kdecay,SP.slim,2,2);
        
    TopoXY = [TopoXY(1,:); Topography'];
end

% Back to old mesh
if ksea~=0
    TopoXY = [Old_topo(1,:); interp1(TopoXY(1,:),TopoXY(2,:),Old_topo(1,:))];
end

% Redimensioning topography
if SP.scale
    Topography = Topography*model_long;
end

% % Plot old and new topographies (uncomment)
% figure(2)
% plot(TopoXY(1,:)/1000,Old_topo/1000,'k')
% hold on
% plot(TopoXY(1,:)/1000,Topography/1000,'r')

function tpelagic = pelagic(SP,TopoXY,dt_sed)
switch SP.pelagic_t
    case 'cnst'
        tpelagic = zeros(1,length(TopoXY));
        tpelagic(TopoXY(2,:)<SP.sealevel) = SP.pelagic_rate*dt_sed;
%         hold on
%         plot(TopoXY(1,:),TopoXY(2,:)+tpelagic,'xr')
    case 'var'
end