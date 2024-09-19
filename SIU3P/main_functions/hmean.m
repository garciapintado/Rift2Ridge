function [hm] = hmean(V)
hm = size(V,1)./sum(1./V);
