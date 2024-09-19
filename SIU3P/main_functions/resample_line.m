function L = resample_line(L0,s,type)
% L = RESAMPLE_LINE(L0,S,TYPE) resamples the line given by L0 with an S 
% spacing. TYPE could be 'linear' or 'spline' which stands for the type of
% interpolation.

d = [0 cumsum(sqrt((L0(1,2:end)-L0(1,1:end-1)).^2 + ...
    (L0(2,2:end)-L0(2,1:end-1)).^2))];

d_start = d(1);
d_end = d(end);

np = round(abs(d_start-d_end)/s);

ds = linspace(d_start,d_end,np);

L = [interp1(d,L0(1,:),ds,type); interp1(d,L0(2,:),ds,type)];

% plot(L0(1,:),L0(2,:),'.-')
% hold on
% plot(L(1,:),L(2,:),'.-')
% hold off