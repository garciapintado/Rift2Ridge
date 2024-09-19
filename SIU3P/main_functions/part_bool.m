function [v_segment] = part_bool(v_bool)

% [V_SEGMENT] = PART_BOOL(V_BOOL) finds in the vector V_BOOL, with values
% of 0 or 1, where the 1s are contiguous, and assigns to the different
% groups of 1s, separated by 0s, a different index. The output V_SEGMENT is
% a vector of the same size of V_BOOL, with 0s where V_BOOL is 0, and with
% the indexes of the 1-groups where V_BOOL is 1.
%
% Example:
%
% v_bool = [1 1 0 0 1 1 1 0 1];
% v_segment = part_bool(v_bool);
% 
% Then v_segment will be:
% v_segment = [1 1 0 0 2 2 2 0 3]

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 03-04-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

v_segment = zeros(size(v_bool));
c = 1;
u = 0;
for j = 1:length(v_bool)
    if v_bool(j) ~= 0
        v_segment(j) = c;
        u = 1;
    elseif u == 1
        c = c+1;
        u = 0;
    end
end