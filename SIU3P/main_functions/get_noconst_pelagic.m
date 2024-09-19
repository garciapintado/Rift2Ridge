function nconstpe = get_noconst_pelagic(belowSea,TopoXY,Pe_mar,Pe_hemi,ccdthk,dh,sealevel)
  % Non pelagic rate as a function to ccdlevel
  ccdlevel = sealevel - ccdthk; % ccdlevel 
  
  acoeff = (Pe_hemi - Pe_mar) / dh;
  bcoeff =  Pe_mar - acoeff*ccdlevel;
  
  Hw       = sealevel - TopoXY(2,belowSea);
  
  nconstpe = acoeff.*Hw + bcoeff;
  
  nconstpe(Hw <= ccdlevel) = Pe_mar;
  nconstpe(Hw >= ccdlevel) = Pe_hemi;
  
  plot(Hw,nconstpe);
  
end
