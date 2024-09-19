function V_nod = ipval2nodval(EL2NOD, V_ip, eids, continuous, quadratic)
% Usage: V_nod = ipval2nodval(EL2NOD, V_ip, eids)
% 
% Purpose: Calculates nodal values of one (or more) variables defined at
%          integration points. Gives a discontinuous (per element)
%          solution.
%
% Input:
%   EL2NOD : [matrix] : finite element connectivity matrix (nnodel x nel)
%   V_ip   : [cell]   : cell array with variable values at integration
%                       points; (nvar x 1), each cell (nel x nip)
%
% Output:
%   V_nod  : [cell]   : - Default: cell array of discontinuous variable fields at
%                       nodes; (nvar x 1), each cell (nnodel x nel)
%                       - Continuous option: [nnod,1] 

% Notes: this function map values at integration point to values at triangle nodes by
%        simple linear regression over the corners of each element. Higher order
%        variation from the integration point is neglected
%
%        TODO: this is a simplistic approach to recovery of nodal values from ip values.  
%              A better procedure, e.g. as the recovery by
%              equilibration of patches (Boroomand and Zienkiewicz, 1997), is due.
%
% Author:       J. García-Pintado: 
%               - augmented to use standard 'triangle' node convention
%               - include 'continuous' option as returned value by simple averaging
%

if nargin < 4
    continuous = false;
end

if nargin < 5
    quadratic = true; % recalculate the 'bubble' 7th node by quadratic shape functions from the continuous values, rather than from linear interpolation 
end

if ~iscell(V_ip)
    V_ip = {V_ip};
    format_V_nod = 'matrix';
else
    format_V_nod = 'cell';
end

[nnodel, nel] = size(EL2NOD);
nip  = size(V_ip{1},2);
nvar = length(V_ip);

if ~ismember(nnodel,[3 6 7])
    error(' Cannot identify element type.');
end   
    
nvert = 3;
IP_X  = ip_triangle(nip);
NL    = sf_dsf_tri367(IP_X',nvert,'matrix');                               % [nvert,nip]

V_nod = cell(1,nvar);

if eids == "456"                     % "kinedyn" edge nodes
    eid = [4 5 6];                                                     
elseif eids == "645"                 % "rift2ridge" edge nodes
    eid = [6 4 5];
else
   error('edge node sorting not identified') 
end

for k=1:nvar                                                               % each variable
    Vi_ip  = V_ip{k};
    Vi_nod = zeros(nnodel,nel);
    Vi_nod(1:nvert,:) = NL' \ Vi_ip';
    if nvert < nnodel
        Vi_nod(eid(1),:) = 0.5.*(Vi_nod(1,:) + Vi_nod(2,:));                    % Interpolate edge node values
        Vi_nod(eid(2),:) = 0.5.*(Vi_nod(2,:) + Vi_nod(3,:));
        Vi_nod(eid(3),:) = 0.5.*(Vi_nod(1,:) + Vi_nod(3,:));
        if nnodel==7
            Vi_nod(7,:) = (1/3).*sum(Vi_nod(1:3,:),1);
        end  
    end
    if continuous
        Vi_nod = nelval2nodval(EL2NOD, Vi_nod); % [nnod,1]
    end
    if quadratic && nnodel==7
        nnod6 = max(EL2NOD(1:6,:),[],'all');
        if continuous
            Vi_nod = addCentroidVal(EL2NOD(1:6,:), Vi_nod(1:nnod6)')';
        else
            Vi_nod_con = nelval2nodval(EL2NOD, Vi_nod);
            Vi_nod_con = addCentroidVal(EL2NOD(1:6,:), Vi_nod_con(1:nnod6)');
            Vi_nod(7,:) = Vi_nod_con((nnod6+1):end);
        end
    end

    mnmx = minmax(Vi_ip(:)');
    Vi_nod(Vi_nod < mnmx(1)) = mnmx(1);
    Vi_nod(Vi_nod > mnmx(2)) = mnmx(2);
    if format_V_nod == "matrix"
        V_nod = Vi_nod;
    else
        V_nod{k} = Vi_nod;
    end
end

end % END OF FUNCTION ipval2nodval
