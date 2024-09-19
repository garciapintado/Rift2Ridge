function L1 = clip2GEOM(L1, L2, clipdis, rmdupl)
    % clip polyline L1 to polyline L2 considering a minimum vertical
    % distance. Note clipping does not project L1 into L2; simply replaces
    % coordinates. Duplicates in the updated L1 are removed.
    % 
    % L{1} :: [2,n1] polyline, mutable
    % L{2} :: [2,n2] polyline, immutable
    % clipdis:: threshold layer thicknes to trigger clipping
    % returns L1a and corresponding indices in L1 input vector, referring to
    % surviving indices after duplicate removal
    % 
    % Sharp direction change check: within-polyline angles are evaluated
    % and a cost function evaluating direction change and Euclidean
    % distance is evaluated to choose destination clipping node
    % in L{2} where the L{1} node in to be clipped into 
    
    % Author: Javier García-Pintado, MARUM, 2020
    if nargin < 4
        rmdupl = true;                                                     % duplicate removal
    end
    
    function cost = costf(dalpha, dd)
       cost = dalpha ./ 180. + dd ./ clipdis;   
    end
    
    L1dz = interp1(L2(1,:),L2(2,:),L1(1,:)) - L1(2,:);                     % L1-L2 thickness at L1 nodes
    L1boo = L1dz <= clipdis;                                               % LOGICAL [n1]    
    L1boo(([0 diff(L1boo)] == 1) & ([diff(L1boo) 0] == -1)) = false;       % do not produce widow clipped nodes
    
    [dd2,I2] = pdist2(L2',L1','euclidean','Smallest',2);                   % [2,n1] two nearest neighbours as candidate destinations for each L{1} node
    L1boo = L1boo & dd2(1,:) > 0.;                                         % remove already clipped nodes from moving set
    
    if ~any(L1boo)
        ids = 1:size(L1,2);
        return;
    end
    
    % stage changes
    L1ids = find(L1boo);                                                   % point indices in L1 within clipdis to L2 
    kids  = repelem(1,length(L1ids));                                      % Euclidean nearest index within search (1 or 2 for 1st/2nd nearest one)
    js = [-1 1];                                                           % -1:junction before node, 1:junction after node
    for ii=1:length(L1ids)                                                 % hold on; plot(L1(1,i)/1000,L1(2,i)/1000,'+','markersize',20)                                       
       done = false;                                                       % hold on; plot(L2(1,I2(1,i))/1000,L2(2,I2(1,i))/1000,'o','markersize',25)
       i = L1ids(ii);                                                      % index within L1
       
       if all(ismember([i-1,i+1],L1ids))                                   % surrounded by collapsing nodes: just collapse into nearest neighbour
           continue;                                                      
       end
       for j=js                                                            % 1st junction criterion
           if (I2(1,i) == I2(1,i+j)) && (dd2(1,i+j) == 0)                    % 1st neighbour is junction => choose this target
               kids(ii) = 1;
               done = true;
               break;
           end
       end
       if done
           continue;
       end
       for j=js                                                            % 2nd junction criterion
           if (I2(2,i) == I2(1,i+j)) && (dd2(1,i+j) == 0)                  % 2nd neighbour is junction (before or after this node)
               xy1 = [L1(:,[i i-j]) L2(:,I2(1,i))];                      
               dalpha1 = getAnglesPolyline(xy1);                            % 1st option angle variation
               xy2 = [L1(:,[i i-j]) L2(:,I2(2,i))];                            
               dalpha2 = getAnglesPolyline(xy2);                            % 2nd option angle variation                             
               if dalpha2 < dalpha1                                        % and angle variation more conservative 
                   kids(ii) = 2;
                   done = true;
                   break;
               end    
           end  
       end % if 2nd junction criterion
       if done
           continue;
       end
       if ismember(i-1,L1ids)                                              % angle + distance criterion  
           j = i + 1;                                                      % io: index in L1 of inmmutable node for angle criterion 
       else
           j = i - 1;
       end   
       xy1 = [L1(:,[i j]) L2(:,I2(1,i))];                      
       dalpha1 = getAnglesPolyline(xy1);                                   % 1st option angle variation
       xy2 = [L1(:,[i j]) L2(:,I2(2,i))];                            
       dalpha2 = getAnglesPolyline(xy2);                                   % 2nd option angle variation      
       cost = costf([dalpha1; dalpha2], dd2(:,i));
       [~,kids(ii)] = min(cost);    
    end   
                                                                           % figure(); plot(L1(1,:),L1(2,:),'.-','color','blue') 
                                                                           % hold on; plot(L2(1,:),L2(2,:),'x-','color',[.4 .4 .4])
                                                                           % hold on; plot(L1(1,L1boo),L1(2,L1boo),'o','color','green') 
    L2ids = I2(sub2ind(size(I2),kids,L1ids));                              % indices in L2 to replace L1 ones
    
    % commit changes
    L1(:,L1ids) = L2(:,L2ids);
    ids = 1:size(L1,2);
    if rmdupl
       %[~,ids,~] =  unique(L1','stable','rows');
       [~,ids,~] =  unique(L1(1,:),'stable');
       L1 = L1(:,ids);
    end
    L = addMissingNodes({L1,L2});                                          % remove widow non-clipped nodes
    L1 = L{1};
end % function
