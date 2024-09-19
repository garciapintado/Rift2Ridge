function LINE = cut(LINE,LINEC,rm_above)

% plot(LINE(1,:),LINE(2,:),'.-k')
% hold on
% plot(LINEC(1,:),LINEC(2,:),'.-r')

[int_points,~] = intersect_lines(LINE,LINEC);

L12 = interp1(LINE(1,:),LINE(2,:),LINEC(1,:));
L21 = interp1(LINEC(1,:),LINEC(2,:),LINE(1,:));

L12I = unique([LINE(1,:),LINEC(1,:),int_points(1,:)]);
L12I = [L12I; repmat(Inf*ones(1,size(L12I,2)),5,1)];
[i,ii] = ismember(L12I(1,:),LINE(1,:));
L12I(2,i) = LINE(2,ii(ii~=0));
[i,ii] = ismember(L12I(1,:),LINEC(1,:));
L12I(3,i) = LINEC(2,ii(ii~=0));
[i,ii] = ismember(L12I(1,:),LINEC(1,:));
L12I(4,i) = L12(ii(ii~=0));
[i,ii] = ismember(L12I(1,:),LINE(1,:));
L12I(5,i) = L21(ii(ii~=0));
[i,ii] = ismember(L12I(1,:),int_points(1,:));
L12I(6,i) = int_points(2,ii(ii~=0));

LINE = [L12I(1,:); min(L12I(2:6,:))];

% LINE = ismember(LINE,L12I)
% plot(LINE(1,:),LINE(2,:),'g')


