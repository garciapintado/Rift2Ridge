for i=1:size(SED.E2N,2)
    is = (i-1)*3+1; ie = (i-1)*3+3;
    EL2N(i,:) = is:ie; 
end

patch('faces',SED.E2N(1:3,SED.Phases==max(SED.Phases))','vertices',SED.GCO'/1000,'facevertexcdata',Temp','FaceColor','flat')
caxis([0 max(max(Temp(SED.E2N(:,SED.Phases==max(SED.Phases)))))])
shading interp
axis tight
colorbar