hold on

Mu_all_m = Mu_all(Phases==4,:);
ELEM2NODE_m = ELEM2NODE(:,Phases==4);

while 1
    visco_min = min(min(Mu_all_m));
    plot(GCOORD(1,ELEM2NODE_m(Mu_all_m==visco_min))/km,...
        GCOORD(2,ELEM2NODE_m(Mu_all_m==visco_min))/km,'.k')
    ELEM2NODE_m(Mu_all_m==visco_min)=[];
    Mu_all_m(Mu_all_m==visco_min)=[];
    pause(0.5)
end