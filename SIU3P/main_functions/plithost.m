function PLITHOST = plithost(GCOORD,ELEM2NODE,Point_id,Phases,...
    rho,G,xip,zip)
% PLITHOS = PLITHOS(GCOORD,ELEM2NODE,POINT_ID,PHASES,RHO,G) calculates the
% lithostatic pressure at the integration points using the mesh defined by
% GCOORD, ELEM2NODE, POINT_ID, PHASES and the density RHO and gravity
% vector G.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 03-09-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------
  error("JGP: I believe there is an error in this function. Check before use...")
  t = tic;

  nip = size(xip,2);

  % Vectorize GIPs
  Gip_x = xip(:);
  Gip_y = zip(:);

  % Find interfaces id
  Interf_id = 3:3:max(Point_id);

  % Initialize vertical distance to interface matrix
  DIST2INT = zeros(length(Interf_id),size(Gip_x,1));

  % calculate vertical distances from Gip to interfaces on top
  for n = 1:length(Interf_id)
      % Interface evaluated
      INTn = GCOORD(:,Point_id==Interf_id(n));
      % Sort interface
      [INTn(1,:),indx] = sort(INTn(1,:));
      INTn(2,:) = INTn(2,indx);
      % Y coordinate of the projection of the global coordinates into the
      % interface evaluated
      INTy = interp1(INTn(1,:),INTn(2,:),Gip_x); % 
   
      DIST2INT(n,:) = INTy'-Gip_y';                                          % vertical distances
  end

  % Remove negative distances (when points are above interfaces)
  DIST2INT(DIST2INT<0) = 0;

  % Build density matrix
  RHO = rho(repmat(unique(Phases)',1,size(Gip_x,1)));
  RHO(1:end-1,:) = -diff(RHO);

  nel = size(ELEM2NODE,2);
  PLITHOST = reshape(sum(DIST2INT.*RHO.*(-G(2))), nel, nip);                 % lithostatic pressure for each integration point
end
 
%disp(toc(t))

% % Plot (uncomment)
% plp = sum(PLITHOST,2)/6;
% patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000,'facevertexcdata', ...
%     plp(:)/1e6,'FaceColor','flat')
% shading flat
% colorbar
% title('Pressures [MPa]')
% xlabel('Distance [Km]')
% ylabel('Depth [Km]')