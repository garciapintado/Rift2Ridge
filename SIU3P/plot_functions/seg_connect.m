function [cont_line] = seg_connect(coord,e2n)
% [CONT_LINE] = SEG_COONECT(COORD,E2N) makes a CONT_LINE 2xNNODES matrix
% with x and y coordinates by connecting segments defined by a connectivity
% matrix E2N which coordinates are COORD.

% Author: Miguel Andres-Martinez, University of Bremen
%         andresma@uni-bremen.de

% Find the two edges of the line
[node_edge,el_edge] = find(reshape(sum(e2n(:)==e2n(:)'),2,size(e2n,2))==1);
prev_el = e2n(el_edge(1));
prev_nod = node_edge(1);
cont_line = coord(:,e2n(:,prev_el));

% Loop through segments
for n = 1:size(e2n,2)-1
    next_nod = [1; 1]; 
    next_nod(prev_nod) = 0;
    next_nod = next_nod==1;
    
    [ii,jj] = find(ismember(e2n,e2n(next_nod,prev_el)));
    el_col = jj~=prev_el;
    if sum(el_col)>1
        el_col = find(el_col);
        el_col = el_col(2);
    end
    prev_el = jj(el_col);
    prev_nod = ii(el_col);
    
    cont_line = [cont_line coord(:,e2n(:,prev_el))];
end


