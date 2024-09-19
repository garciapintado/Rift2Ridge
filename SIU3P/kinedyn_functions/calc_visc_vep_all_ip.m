function [Visc_els,Visc_v,Visc_e,Visc_p,Gamma,YC] = calc_visc_vep_all_ip...
    (VAR,SETTINGS,PHYSICS,NUMSCALE,GCOORD,EL2NOD,EL2NODP,PhaseID,els,N3_ip,N6_ip,dt)

nel = length(els);
nip = length(N3_ip);
if length(N6_ip)~=nip
    error('3- and 6-node shape functions must be provided for the same number of stress evaluation points');
end

% =========================================================================
% CALCULATE "VISCOUS" VISCOSITY AT INTEGRATION POINT
% =========================================================================
switch SETTINGS.method_eval_visc
    case 'interp_nodal'
        % interpolation from nodal values to integration points
        nnodel = length(N3_ip{1});
        Visc_v = zeros(nel,nip);
        for ip=1:nip
            Visc_v(:,ip) = VAR.Visc( EL2NOD(1:nnodel,els) )'*N3_ip{ip};
        end
        
    case 'mean_nodal'
        % harmonic mean for each element
        nnodel = size(EL2NOD,1);
        if nnodel==7 % no 7th node in viscosity vector
            nnodel = nnodel - 1;
        end
        Visc_v = VAR.Visc( EL2NOD(1:nnodel,els) );
        Visc_v = nnodel./sum(1./Visc_v,1)';
        Visc_v = repmat(Visc_v(:),1,nip);
        
    case 'elem_phases'
        % use the constant viscosities defined in PHYSICS.Visc
        Visc_v = PHYSICS.Visc( PhaseID(els) );
        Visc_v = repmat(Visc_v(:),1,nip);
        
    case 'elem_var'
        if size(VAR.Visc,2)==1
            Visc_v = VAR.Visc(els);
            Visc_v = repmat(Visc_v(:),1,nip);
        elseif size(VAR.Visc,2)==nip
            Visc_v = VAR.Visc(els,:);
        else
            error('VAR.Visc must be of size [nel x 1] or [nel x nip]');
        end
        
    case 'powlaw'
        % -------------------------------------------------------------------------
        % The calculated strain rates and stresses are in units defined by NUMSCALE
        % To convert them to SI units you need to "unscale" them again:
        % strain_SI = VAR.Strain / NUMSCALE.t0;
        % stress_SI = VAR.Stress * (NUMSCALE.Visc0/NUMSCALE.t0);
        % -------------------------------------------------------------------------
        if ~isfield(PHYSICS,'RHEOL')
            error('Powerlaw rheology requires the structure "PHYSICS.RHEOL".');
        end
        RHEOL    = PHYSICS.RHEOL;

        
%         % Parameters for dislocation creep
%         Adis     = RHEOL.Adis(PhaseID(els),:);
%         Ndis     = RHEOL.Ndis(PhaseID(els),:);
%         Qdis     = RHEOL.Qdis(PhaseID(els),:);
%         Vdis     = RHEOL.Vdis(PhaseID(els),:);
%         
%         % Parameters for diffusion creep
%         Adif     = RHEOL.Adif(PhaseID(els),:);
%         Ndif     = RHEOL.Ndif(PhaseID(els),:);
%         Qdif     = RHEOL.Qdif(PhaseID(els),:);
%         Vdif     = RHEOL.Vdif(PhaseID(els),:);


        % Parameters for dislocation creep
        Adis     = PHYSICS.EL.Adis(els);
        Ndis     = PHYSICS.EL.Ndis(els);
        Qdis     = PHYSICS.EL.Qdis(els);
        Vdis     = PHYSICS.EL.Vdis(els);
        
        % Parameters for diffusion creep
        Adif     = PHYSICS.EL.Adif(els);
        Ndif     = PHYSICS.EL.Ndif(els);
        Qdif     = PHYSICS.EL.Qdif(els);
        Vdif     = PHYSICS.EL.Vdif(els);

        
        n        = size(Adis,2);
        C2K      = 273.15;
        R        = 8.314472; % J/mol K; univ. gas constant
        nTnodel  = 6;
        nPnodel  = 3;
        Visc_dis = zeros(nel,n);
        Visc_dif = zeros(nel,n);
        Visc_v   = zeros(nel,nip);
        
        % =========================== DEBUGGING ===========================
        debug_me = 0;
        if debug_me
            x_monitor     = 450;
            z_monitor     =  -5;
            FigNo_monitor = 666;

            % Find itegration point that is closest to monitoring coordinates
            x_ip_all = zeros(nip,nel);
            z_ip_all = zeros(nip,nel);
            for ip=1:nip
                x_ip_all(ip,:) = (reshape(GCOORD(1,EL2NOD(1:nPnodel,els)),nPnodel,nel)' * N3_ip{ip});
                z_ip_all(ip,:) = (reshape(GCOORD(2,EL2NOD(1:nPnodel,els)),nPnodel,nel)' * N3_ip{ip});
            end
            [dist,indx_ip] = sort(sqrt((x_ip_all(:)-x_monitor).^2 + (z_ip_all(:)-z_monitor).^2));
            [ip_monitor,iel_monitor] = ind2sub([nip nel],indx_ip(1));
            x_monitor = x_ip_all(ip_monitor,iel_monitor);
            z_monitor = z_ip_all(ip_monitor,iel_monitor);
            clear x_ip_all z_ip_all
        end
        % =========================== DEBUGGING ===========================
        
        for ip=1:nip
            % Evaluate pressure at integration point
            % NOTE: PRESSURE MUST BE IN UNITS OF "Pa"
            if isfield(VAR,'P_lithstat')
                P_lithstat = VAR.P_lithstat( EL2NODP(1:nPnodel,els) )';
            elseif strcmp(SETTINGS.top_surface,'free')
                P_lithstat = VAR.P( EL2NODP(1:nPnodel,els) )';
            else
                error('Lithostatic pressure must be provided.');
            end
             P_ip = P_lithstat * N3_ip{ip};
            Pd_blk = VAR.P(EL2NODP);
            Pd_blk = Pd_blk' * N3_ip{ip};
            P_ip = Pd_blk * NUMSCALE.P0; % must be in units of [Pa] !!!
            % z_ip = reshape(GCOORD(2,EL2NOD(1:nPnodel,els)),nPnodel,[])' * N3_ip{ip}; % depth of IP
            % figure(88);clf;plot(P_ip,z_ip,'k.');

            % Evaluate temperature at integration point
            if isfield(VAR,'T')
                T_ip = VAR.T( EL2NOD(1:nTnodel,els) )' * N6_ip{ip} + C2K;
            else
                T_ip = PHYSICS.T0 * ones(nel,1) + C2K;
            end
            
%             % Evaluate hydration-level at integration point
%             % Better use linear variation here (i.e. N3_ip not N6_ip) to
%             % avoid over- and undershoots
%             Hyd_ip = VAR.Hyd( EL2NOD (1:3,els) )'*N3_ip{ip};
%             if any(Hyd_ip>1) || any(Hyd_ip<0)
%                 error('Element hydration level "VAR.Hyd" is out of plausible range [0,1].');
%             end

            % 2nd strain rate invariant in units of 1/s
            if all(VAR.Er_II_totl(:)==0)
                % In first iteration of first time step only:
                Er_II    = ones(nel,nip) .* PHYSICS.ext_rate./(max(GCOORD(1,:))-min(GCOORD(1,:)));
                Er_II_SI = Er_II ./ NUMSCALE.t0;
            else
                Er_II_SI = VAR.Er_II_totl ./ NUMSCALE.t0;
            end
            
            idis   = Ndis(:,1)>0;
            idif   = Ndif(:,1)>0;
            
            if strcmp(SETTINGS.use_strainsoft,'yes')
                Eacc_II_ip        = VAR.Eacc_II(els,ip);
                Eacc_II_ip        = VAR.Eacc_visc(els,ip);
                [Adis_ss,Adif_ss] = preexp_factor_strain_softening(PHYSICS.SS,Eacc_II_ip,T_ip);
                
                % =========================================================
                % Only update the phases selected for strain softening 
                % =========================================================
                iel_noSS          = ~ismember(PhaseID(els),PHYSICS.SS.PhasesID_ss)';
                Adis_ss(iel_noSS) = 1;
                Adif_ss(iel_noSS) = 1;
            else
                Adis_ss = ones(nel,1);
                Adif_ss = ones(nel,1);
            end
            
            for i=1:n
                % Dislocation creep
                % =================
                % Factor for scaling triaxial and uniaxial experiment parameters
                % (see Gerya, 2010)
                Sc_dis        = 1 ./ (   2.^((Ndis(idis,i)-1)./   Ndis(idis,i)) ...
                                      .* 3.^((Ndis(idis,i)+1)./(2*Ndis(idis,i))) );
%                 Sc_dis        = 1/2; % factor used in Kinedyn_v1
                Visc_dis(idis,i) = Sc_dis ...
                    .* (Adis_ss(idis).*Adis(idis,i)) .^ (-1./Ndis(idis,i)) ...
                    .* Er_II_SI(idis,ip).^(1./Ndis(idis,i)-1) ...
                    .* exp( (Qdis(idis,i)+P_ip(idis).*Vdis(idis,i))./(Ndis(idis,i).*R.*T_ip(idis)) );
                
                % Diffusion creep
                % ===============
                % Factor for scaling triaxial and uniaxial experiment parameters
                % (see Gerya, 2010)
                Sc_dif           = 1/3;
%                 Sc_dif           = 1/2; % factor used in Kinedyn_v1
                Visc_dif(idif,i) = Sc_dif ...
                    .* (Adif_ss(idif).*Adif(idif,i)) .^ (-1./Ndif(idif,i)) ...
                    .* (RHEOL.Grain/RHEOL.Burger) ...
                    .* exp( (Qdif(idif,i)+P_ip(idif).*Vdif(idif,i))./(Ndif(idif,i).*R.*T_ip(idif)) );
            end
            if n==2
                error('Hydration not yet implemented.');
                % Geometric mean between wet and dry viscosities
                tmp            =    Hyd_ip(idis) .*Visc_dis(idis,1) ...
                               + (1-Hyd_ip(idis)).*Visc_dis(idis,2);
                Visc_dis       = zeros(nel,1);
                Visc_dis(idis) = tmp;
                
                tmp            =    Hyd_ip(idif) .*Visc_dif(idif,1) ...
                               + (1-Hyd_ip(idif)).*Visc_dif(idif,2);
                Visc_dif       = zeros(nel,1);
                Visc_dif(idif) = tmp;
                
% %                 % Weighted harmonic mean between wet and dry viscosities
% %                 tmp            =    X_ip(idis) ./Visc_dis(idis,1) ...
% %                                + (1-X_ip(idis))./Visc_dis(idis,2);
% %                 Visc_dis       = zeros(nel,1);
% %                 Visc_dis(idis) = 1./tmp;
% %                 
% %                 tmp            =    X_ip(idif) ./Visc_dif(idif,1) ...
% %                                + (1-X_ip(idif))./Visc_dif(idif,2);
% %                 Visc_dif       = zeros(nel,1);
% %                 Visc_dif(idif) = 1./tmp;
            end
            inone            = find(~idis & ~idif);
            if ~isempty(inone)
                disp(['PhaseID: ' num2str(unique(PhaseID(inone)))])
                error('There are elements (PhaseID above) with neither dislocation nor diffusion creep active');
            end
            iboth            = idis & idif;
            Visc_v(iboth,ip) = 1./(1./Visc_dis(iboth) + 1./Visc_dif(iboth));
            idis(iboth)      = false;
            Visc_v(idis,ip)  = Visc_dis(idis);
            idif(iboth)      = false;
            Visc_v(idif,ip)  = Visc_dif(idif);
            
            % =========================== DEBUGGING ===========================
            if debug_me && ip==ip_monitor
                figure(FigNo_monitor);
                subplot(5,1,1);
                iplot = length(get(gca,'Children'))+1;
                plot(iplot,T_ip(iel_monitor),'k.');
                if iplot==1
                    hold on; grid on;
                    xlabel('iteration');ylabel('T (C)');
                end
                subplot(5,1,2);
                plot(iplot,1e-9*P_ip(iel_monitor),'k.');
                if iplot==1
                    hold on; grid on;
                    xlabel('iteration');ylabel('P (GPa)');
                end
                subplot(5,1,3);
                plot(iplot,log10(Er_II_SI(iel_monitor)),'k.');
                if iplot==1
                    hold on; grid on;
                    xlabel('iteration');ylabel('log Er_{II} (1/s)');
                end
                subplot(5,1,4);
                plot(iplot,log10(Visc_dis(iel_monitor)),'rx');
                if iplot==1
                    hold on; grid on;
                    xlabel('iteration');ylabel('Viscosity (log Pa s)');
                end
                plot(iplot,log10(Visc_dif(iel_monitor)),'bx');
                plot(iplot,log10(Visc_v  (iel_monitor,ip_monitor)),'k.');
                if iplot==1
                    legend('dislocation','diffusion','effective');
                end
                subplot(10,1,9);
                plot(iplot,x_monitor,'k.');
                if iplot==1
                    hold on; grid on;
                    xlabel('iteration');ylabel('x-coordinate (km)');
                end
                subplot(10,1,10);
                plot(iplot,z_monitor,'r.');
                if iplot==1
                    hold on; grid on;
                    xlabel('iteration');ylabel('z-coordinate (km)');
                end
                drawnow
            end
            % =========================== DEBUGGING ===========================
        end
        
        % SCALE VISCOSITIES ACCORDING TO REFERENCE VISCOSITY
        Visc_v           = Visc_v ./ NUMSCALE.Visc0;
        
    otherwise
        error(' Unknown "case" for evaluating viscosity at integration points.');
end

% VISCOSITY CUT-OFFS (IF DEFINED)
if isfield(PHYSICS,'minVisc') && ~isempty(PHYSICS.minVisc)
    Visc_v = max(Visc_v,PHYSICS.minVisc); % lower cut-off for viscosity
end
if isfield(PHYSICS,'maxVisc') && ~isempty(PHYSICS.maxVisc)
    Visc_v = min(Visc_v,PHYSICS.maxVisc); % upper cut-off for viscosity
end

% =========================================================================
% CALCULATE "ELASTIC" VISCOSITY
% =========================================================================
if strcmp(SETTINGS.is_elastic,'yes') && dt>0
    ShearG_els = PHYSICS.ShearG(PhaseID(els));
    Visc_e     = dt*ShearG_els(:);
    Visc_e     = repmat(Visc_e(:),1,nip);
else
    Visc_e     = [];
end

% =========================================================================
% PUT ALL TOGETHER TO CALCULATE EFFECTIVE VISCOSITY (diffreent versions)
% =========================================================================
if strcmp(SETTINGS.is_elastic,'no') && strcmp(SETTINGS.is_plastic,'no')
    % ONLY VISCOUS FLOW (STOKES FLOW)
    Visc_els = Visc_v;
else
    % ELASTICITY AND/OR PLASTICITY CONSIDERED
    if nip~=size(VAR.Er_II_totl,2)
        error('Number of points with shape function values does not match number of stresses per element.');
    end
    
    switch SETTINGS.rheology_model
        % =================================================================
        % MAXWELL RHEOLOGY MODEL
        % =================================================================
        case 'maxwell' % e.g. Spiegelman et al. (2016)
                       % all creep mechanisms are additive everywhere
            if strcmp(SETTINGS.is_plastic,'yes')
                Visc_p = VAR.Tau_yld(els,:) ./ (2.*VAR.Er_II_totl(els,:));
                Visc_p(VAR.Er_II_totl(els,:)==0) = PHYSICS.maxVisc;
                % VISCOSITY CUT-OFFS (IF DEFINED)
                if isfield(PHYSICS,'minVisc') && ~isempty(PHYSICS.minVisc)
                    Visc_p = max(Visc_p,PHYSICS.minVisc); % lower cut-off for viscosity
                end
                if isfield(PHYSICS,'maxVisc') && ~isempty(PHYSICS.maxVisc)
                    Visc_p = min(Visc_p,PHYSICS.maxVisc); % upper cut-off for viscosity
                end
            else
                Visc_p = [];
            end

            % Calculate effective viscosity assuming that strain rates of all 
            % creep mechanisms are additive
            if ~isempty(Visc_p)
                if ~isempty(Visc_e)
                    Visc_els = 1 ./ (1./Visc_v + 1./Visc_e + 1./Visc_p);
                else
                    Visc_els = 1 ./ (1./Visc_v + 1./Visc_p);
                end
            else
                if ~isempty(Visc_e)
                    Visc_els = 1 ./ (1./Visc_v + 1./Visc_e);
                else
                    Visc_els = Visc_v;
                end
            end
            
        case 'moresif' % Moresi et al. (2003)
            % =================================================================
            % MORESI RHEOLOGY MODEL PLUS STRESS-FORCAST
            % =================================================================
            if strcmp(SETTINGS.is_plastic,'yes')
                if ~isempty(Visc_e)
                    Visc_els       = 1 ./ (1./Visc_v + 1./Visc_e);
                    Visc_els=max(Visc_els,PHYSICS.minVisc);
                    Visc_els=min(Visc_els,PHYSICS.maxVisc);
                    Xi             = Visc_els./Visc_e; % ==Visc_els/(ShearG*dt)
                    Tau_xx_forcast =   4/3.*Visc_els.*VAR.Er_xx(els,:) ...
                                     - 2/3.*Visc_els.*VAR.Er_zz(els,:) ...
                                     +  Xi.*VAR.Tau_xx_old(els,:);
                    Tau_zz_forcast = - 2/3.*Visc_els.*VAR.Er_xx(els,:) ...
                                     + 4/3.*Visc_els.*VAR.Er_zz(els,:) ...
                                     +  Xi.*VAR.Tau_zz_old(els,:);
                    Tau_xz_forcast =     2.*Visc_els.*VAR.Er_xz(els,:) ...
                                     +  Xi.*VAR.Tau_xz_old(els,ip);
                else
                    Visc_els       = Visc_v;
                    Visc_els=max(Visc_els,PHYSICS.minVisc);
                    Visc_els=min(Visc_els,PHYSICS.maxVisc);
                    Tau_xx_forcast =   4/3.*Visc_els.*VAR.Er_xx(els,:) ...
                                     - 2/3.*Visc_els.*VAR.Er_zz(els,:);
                    Tau_zz_forcast = - 2/3.*Visc_els.*VAR.Er_xx(els,:) ...
                                     + 4/3.*Visc_els.*VAR.Er_zz(els,:);
                    Tau_xz_forcast =     2.*Visc_els.*VAR.Er_xz(els,:); 
                end
                Tau_II_els = calc_invariant_II(Tau_xx_forcast,Tau_zz_forcast,Tau_xz_forcast);
            
                if ~isempty(Visc_e)
                    % Calculate effective strain rate; Eq.(38) in (Moresi et al., 2003)
                    % Note that Visc_e == ShearG*dt
                    Er_eff_xx = 2.*VAR.Er_xx(els,:) + 1./Visc_e.*VAR.Tau_xx_old(els,:);
                    Er_eff_zz = 2.*VAR.Er_zz(els,:) + 1./Visc_e.*VAR.Tau_zz_old(els,:);
                    Er_eff_xz = 2.*VAR.Er_xz(els,:) + 1./Visc_e.*VAR.Tau_xz_old(els,:);
                else
                    Er_eff_xx = 2.*VAR.Er_xx(els,:);
                    Er_eff_zz = 2.*VAR.Er_zz(els,:);
                    Er_eff_xz = 2.*VAR.Er_xz(els,:);
                end
                Er_II_eff = calc_invariant_II(Er_eff_xx,Er_eff_zz,Er_eff_xz);
                Visc_p    = VAR.Tau_yld(els,:) ./ Er_II_eff;
                Gamma=VAR.Tau_yld(els,:) .* (1./Visc_p - 1./Visc_els);
                YC=zeros(size(Er_II_eff));
                for ip=1:nip
                    ind_p = Tau_II_els(:,ip)>VAR.Tau_yld(els,ip) & Er_II_eff(els,ip)>0;
                    Visc_els( ind_p,ip) = Visc_p( ind_p,ip);
                    if ~isempty(Visc_e)
                        Visc_els(~ind_p,ip) = 1./(1./Visc_v(~ind_p,ip) + 1./Visc_e(~ind_p,ip));
                    else
                        Visc_els(~ind_p,ip) = Visc_v(~ind_p,ip);
                    end
                    YC(ind_p,ip)=1;
                end
            else
                if ~isempty(Visc_e)
                    Visc_els = 1 ./ (1./Visc_v + 1./Visc_e);
                else
                    Visc_els = Visc_v;
                end
                Visc_p = [];
            end
            
            Visc_els=max(Visc_els,PHYSICS.minVisc);
            Visc_els=min(Visc_els,PHYSICS.maxVisc);

        otherwise
            error('SETTINGS.rheology_model must be either "maxwell" or "harmonic"');
    end
end

% VISCOSITY CUT-OFFS (IF DEFINED)
if isfield(PHYSICS,'minVisc') && ~isempty(PHYSICS.minVisc)
    Visc_els = max(Visc_els,PHYSICS.minVisc); % lower cut-off for viscosity
end
if isfield(PHYSICS,'maxVisc') && ~isempty(PHYSICS.maxVisc)
    Visc_els = min(Visc_els,PHYSICS.maxVisc); % upper cut-off for viscosity
end

% Check if that viscosities are within valid range
if any(Visc_els<=0)
    error(' Calculated a zero or negative viscosity. STOPPING.');
end

end % END OF FUNCTION calc_visc_vep_all_ip