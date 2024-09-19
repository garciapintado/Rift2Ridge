
function [F_xx_new, F_xy_new, F_yx_new, F_yy_new, GIP_X, GIP_Y, I2, ...
          TAU_xx_new, TAU_xy_new, TAU_yy_new, Mu_all_new, E2_all_new] = ...
    remesh_F_TAU_Mu_E2(GCOa, E2Na, GCOb, E2Nb, ...
        F_xx, F_xy, F_yx, F_yy, TAU_xx, TAU_xy, TAU_yy, Mu_all, E2_all, I2, ...
        extrap_scheme)
  %
  % GCOa :: updated coordinates
  % E2Na :: updated element to node maping
  % GCOb :: background coordinates
  % E2Nb :: background element to node mapping
  % F_xx, F_xy, F_yx, F_yy, TAU_xx, TAU_xy, TAU_yy, Mu_all, I2 :: variables
  %   (I2 is a struct) defined at integration points in the background mesh
  % extrap_scheme :: string, in {"ip2nodb2ip_N3","ip2nodb2ip_N7"}
  %
  % First version: Miguel Martinez
  %
  % Last modification: Javier GP, MARUM, 2020 
  %               - restructured by iterating through
  %                 variables: introduces a few loops albeit simpler, and memory
  %                 safer.
  %               - nonlinear scheme introduced into 'ipsb2nodesb2ipsa'
  %                  interpolation cheme
  
  % Switch for the interpolation scheme
  
  I2_p = I2.p;
  I2_c = I2.c;

  eids = "645";                                                               % rift2ridge convention
  nela = size(E2Na,2);
  nelb = size(E2Nb,2);                                                        % background: number of elements
  vnames = ["F_xx" "F_xy" "F_yx" "F_yy" "TAU_xx" "TAU_xy" "TAU_yy" "Mu_all" "E2_all" "I2_p" "I2_c"];
  nv = length(vnames);
  nip = size(F_xx,2);
  
  % Calculates physical coordinates of the integration points for the new mesh
  [GIP_X, GIP_Y] = ip_coord(GCOa, E2Na, nela, nip);                        % [nela,nip]

  if contains(extrap_scheme,"nodb2ip")                                     % background IPs -> background nodes -> updated IPs
      if contains(extrap_scheme,"N3")
          linearint = true;
      elseif contains(extrap_scheme,"N7")
          linearint = false;
      else
          error("remesh_F_TAU_Mu_E2:: unknown extrap_scheme")
      end
      % Old element indexes in which the new ips are contained
      Tris_F = tsearch2(GCOb, uint32(E2Nb(1:3,:)),[GIP_X(:)';GIP_Y(:)']);  % [1, nela*nip]
      Tris_F = reshape(Tris_F,size(GIP_X));                                % [nela, nip]
      Ind = find(Tris_F==0);
        
      % Check all elements were found, and if not find the closest element center
      if ~isempty(Ind)
          for i=1:length(Ind)
              [~, Tris_F(Ind(i))] = min(sqrt( ...
                  (GCOb(1,E2Nb(7,:)) - GIP_X(Ind(i))).^2 + ...
                  (GCOb(2,E2Nb(7,:)) - GIP_Y(Ind(i))).^2));
          end
      end
        
      if(any(isnan(Tris_F)))
          error('remeshing failed in finding closest element');
      end
        
      % Calculates the local coordinates of the new ip
      xp = GIP_X(:)';                                                    % [1, nela*nip] new ip physical coordinates
      yp = GIP_Y(:)';                                                    % [1, nela*nip] new ip physical coordinates
        
      x = reshape(GCOb(1,E2Nb(1:3,Tris_F)),3, numel(xp));                % [3, nela*nip] reshaped old vertices surrounding each new ip
      y = reshape(GCOb(2,E2Nb(1:3,Tris_F)),3, numel(yp));                % [3, nela*nip] reshaped old vertices surrounding each new ip
        
      xi = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:) ...
          .*yp-xp.*y(3,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:) ...
          +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));                % [1, nela*nip] canonical x [$\xi$]  coord for new ip within background triangles
      yi = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp ...
          .*y(1,:)-xp.*y(2,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)...
          +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));                % [1, nela*nip] canonical y [$\eta$] coord for new ip within background triangles
        
      % canonical shape function evaluation at new integration points
      if linearint
          N3 = shp_triangle([xi' yi'], 3);                                   % [3, nela*nip] linear shape functions
      else
          N7 = shp_triangle([xi' yi'], 7);                                   % [7, nela*nip] quadratic+bubble node shape functions
      end
      
      % linear interpolation of variables from:
      % a) old interpolation points into old triangle vertices
      % b) old triangle vertices into new ip coordinates
       
      for iv = 1:nv
           if linearint
               Vn  = ipval2nodval(E2Nb(1:3,:), eval(vnames(iv)), eids);   % [3,nelb]   this variable at old vertices
               Vip = reshape(sum(N3 .* Vn(:,Tris_F)),nela,nip);            % [nela,nip] this variable at new integration points
           else
               Vn = ipval2nodval(E2Nb, eval(vnames(iv)), eids, true);     % [nnod7b,1] as continuous by simple averaging
               Vn = Vn(E2Nb);                                              % [7,nelb] 
               Vip = reshape(sum(N7 .* Vn(:,Tris_F)),nela,nip);            % [nela,nip] this variable at new integration points
           end
           eval(vnames(iv) + "_new = Vip;")
      end
  else
      error("remesh_F_TAU_Mu_E2: unknown extrap_scheme")
  end % if
  
  I2 = [];      % I2 at new integration points
  I2.f = F2I2f(F_xx_new, F_yy_new, F_xy_new, F_yx_new);
  I2.p = I2_p_new;
  I2.c = I2_c_new;  

end % function
