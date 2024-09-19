function [Tei] = Te_Lowry(MESH,VAR,PHYSICS,INPUT_COL)

nnc=INPUT_COL.nnc;

st=0.0;

Beta_N = 0.75;
Beta_T = 3;

lambda=0.42;

pr=0.25;

ext=zeros(length(INPUT_COL.rev_depth_profile),1);
cmp=zeros(length(INPUT_COL.rev_depth_profile),1);
s0=zeros(length(INPUT_COL.rev_depth_profile),1);
sxx=zeros(length(INPUT_COL.rev_depth_profile),1);

amom=0.0;


rate_c=PHYSICS.ext_rate*PHYSICS.km2m/PHYSICS.Myr2s/(MESH.xmax-MESH.xmin);
%rate_c=INPUT_COL.rate;
rate_m=rate_c;

%upper crust
ac=PHYSICS.RHEOL.Adis(3);
nc=PHYSICS.RHEOL.Ndis(3);
qc=PHYSICS.RHEOL.Qdis(3);

%lower crust
acl=PHYSICS.RHEOL.Adis(2);
ncl=PHYSICS.RHEOL.Ndis(2);
qcl=PHYSICS.RHEOL.Qdis(2);

%mantle
am=PHYSICS.RHEOL.Adis(1);
nm=PHYSICS.RHEOL.Ndis(1);
qm=PHYSICS.RHEOL.Qdis(1);

tc_up=INPUT_COL.tc_up;
tc_low=INPUT_COL.tc_low;
tm=INPUT_COL.tm;

st0=0;
r=8.31;

a1=0;


mu=tan(PHYSICS.PhiFric0(4)/180*3.14);
rp=1.0/((sqrt(1.0+mu*mu)-mu)^2);

ec=PHYSICS.ShearG(3)*(2*(1+0.25));
em=PHYSICS.ShearG(1)*(2*(1+0.25));

amom=0;

for j=1:nnc

    dep=-INPUT_COL.rev_depth_profile(j)*1000;
    dev=(((rate_c/ac)^(1.0/nc)))*exp(qc/(nc*r*INPUT_COL.rev_T_profile(j)));
    
%    dev=dev*1.0e6;
    
    cmp(j)=(rp-1.0)*INPUT_COL.rev_Litho(j)*(1.0-lambda);         %!compress brittle strength
    ext(j)=((rp-1.0)/rp)*INPUT_COL.rev_Litho(j)*(1.0-lambda);    %!ext brittle strength
    cmp(j)=-min(cmp(j),dev);                      %!compressional YSE
    ext(j)=min(ext(j),dev);                       %!extensional YSE

    if(dev>2.0e5) 
        nnn=j;
    end
    
    if(st0>0.0)
        s0(j)=min(st0,ext(j));
    else
        s0(j)=max(st0,cmp(j));
    end
end

ndim=nnn;
ndim_ct=nnn;
sgrad=ec*INPUT_COL.curv;

tmom=1.0e20;
zn=0.0;
space=INPUT_COL.space;

%%call integ(ndim,st0,space,zn,sgrad,amom,a1)      

[amom,a1,ext,cmp,s0,sxx]=integ(ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

dmom=amom-2.0*a1;
fxp=zn;
fyp=dmom;
zn=space*double(ndim);

%%call integ(ndim,st0,space,zn,sgrad,amom,a1);
[amom,a1,ext,cmp,s0,sxx]=integ(ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

dmom=amom-2.00*a1;
fxm=zn;
fym=dmom;

while(fym~=0.0)
    fxc=(fyp*fxm-fxp*fym)/(fyp-fym);
    zn=fxc;
    %%call integ(ndim,st0,space,zn,sgrad,amom,a1)
    [amom,a1,ext,cmp,s0,sxx]=integ(ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);
    dmom=amom-2.0*a1;
    fyc=dmom;
    
    if(dmom>0.0)
        fxp=fxc;
        fyp=fyc;
    else
        fxm=fxc;
        fym=fyc;
    end
    
    if(abs(dmom)<=1.0e28)
        break;
    end
end

%call integ(ndim,st0,space,zn,sgrad,amom,a1)
[amom,a1,ext,cmp,s0,sxx]=integ(ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);
if (amom<0)
    amom=amom*(-1);
end

tec=(12.0*amom*(1.0-pr*pr)/(abs(INPUT_COL.curv)*ec))^(1.0/3.0);
znc=zn;


% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c 	Calculate lower crustal component of Te from the geotherm
% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
nncl=INPUT_COL.nncl;

mu=tan(PHYSICS.PhiFric0(3)/180*3.14);
rp=1.0/((sqrt(1.0+mu*mu)-mu)^2);

ec=PHYSICS.ShearG(2)*(2*(1+0.25));      
amom=0.0;

for j=nnc:nncl
      
    dep=-INPUT_COL.rev_depth_profile(j)*1000;
    dev=(((rate_c/acl)^(1.0/ncl)))*exp(qcl/(ncl*r*INPUT_COL.rev_T_profile(j)));
    
%    dev=dev*1.0e6;
    
    cmp(j)=(rp-1.0)*INPUT_COL.rev_Litho(j)*(1.0-lambda);         %!compress brittle strength
    ext(j)=((rp-1.0)/rp)*INPUT_COL.rev_Litho(j)*(1.0-lambda);    %!ext brittle strength
    cmp(j)=-min(cmp(j),dev);                      %!compressional YSE
    ext(j)=min(ext(j),dev);                       %!extensional YSE
 
    if(dev>2.0e5)
        nnn=j;
    end
    if(st0>0)
        s0(j)=min(st0,ext(j));
    else
        s0(j)=max(st0,cmp(j));
    end
end

ndim=nnn;
ndim_ctl=nnn;
sgrad=ec*INPUT_COL.curv;     
    
tmom=1.0e20;

zn=tc_up;

%call integm(nnc,ndim,st0,space,zn,sgrad,amom,a1)      
[amom,a1,ext,cmp,s0,sxx]=integm(nnc,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

dmom=amom-2.0*a1;
fxp=zn;
fyp=dmom;
zn=space*double(ndim);

%call integm(nnc,ndim,st0,space,zn,sgrad,amom,a1)      
[amom,a1,ext,cmp,s0,sxx]=integm(nnc,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);
dmom=amom-2.0*a1;
fxm=zn;
fym=dmom;

while(fym~=0.0)
    fxc=(fyp*fxm-fxp*fym)/(fyp-fym);
    zn=fxc;
    %call integm(nnc,ndim,st0,space,zn,sgrad,amom,a1) 
    [amom,a1,ext,cmp,s0,sxx]=integm(nnc,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);
    dmom=amom-2.0*a1;
    fyc=dmom;
    
    if(dmom>0)
        fxp=fxc;
        fyp=fyc;
    else
        fxm=fxc;
        fym=fyc;
    end
    
    if(abs(dmom)<=1.0e28) 
        break;
    end
end

%call integm(nnc,ndim,st0,space,zn,sgrad,amom,a1)      
[amom,a1,ext,cmp,s0,sxx]=integm(nnc,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

if (amom<0)
    amom=amom*(-1);
end

tecl=(12.0*amom*(1.0-pr*pr)/(abs(INPUT_COL.curv)*ec))^(1.0/3.0);
zncl=zn;


% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c******* Calculate mantle component of Te from the geotherm *******
% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

nnm=int16(tm/space)+1;
nnm=length(INPUT_COL.rev_T_profile);
amom=0;
for j=nncl:nnm
    dev=(((rate_m/am)^(1.0/nm)))*exp(qm/(nm*r*INPUT_COL.rev_T_profile(j)));
	if(dev>10e6)
        l=j;
    end
%    dev=dev*1.0e6;
end
%%%%%%%%%%%%%%%%suspicious!!!!!!!!!!!!!!!!! nnc--> nncl
for j=nncl:l
    dep=-space*double(j-1)*1000; %          !depth is in meters
      
%c		Calculate yield strength envelope ************************

	dev=(((rate_m/am)^(1.0/nm)))*exp(qm/(nm*r*INPUT_COL.rev_T_profile(j)));
%	dev=dev*1.e6;

    cmp(j)=(rp-1.0)*INPUT_COL.rev_Litho(j)*(1.0-lambda);         %!compress brittle strength
    ext(j)=((rp-1.0)/rp)*INPUT_COL.rev_Litho(j)*(1.0-lambda);    %!ext brittle strength
        
	cmp(j)=-min(cmp(j),dev);     
    ext(j)=min(ext(j),dev);      

	if(dev>2.0e5) 
        nnn=j;
    end
    
    if(st0>0.0)
         s0(j)=min(st0,ext(j));     
    else
        s0(j)=max(st0,cmp(j));
    end
end

ndim=nnn;
ndim_mt=nnn;
sgrad=ec*INPUT_COL.curv; 
tmom=1e20;
%cA6   zn=tc                       !change 
zn=tc_up+tc_low;   %                    !change 
      
%call integm(nncl,ndim,st0,space,zn,sgrad,amom,a1)      

[amom,a1,ext,cmp,s0,sxx]=integm(nncl,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

 
dmom=amom-2.0*a1;
fxp=zn;
fyp=dmom;
zn=space*double(ndim);

%call integm(nncl,ndim,st0,space,zn,sgrad,amom,a1)
 [amom,a1,ext,cmp,s0,sxx]=integm(nncl,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

 
dmom=amom-2.0*a1;
fxm=zn;
fym=dmom;

	while(fym~=0.0)
        fxc=(fyp*fxm-fxp*fym)/(fyp-fym);
        zn=fxc;

 %       	call integm(nncl,ndim,st0,space,zn,sgrad,amom,a1)
        [amom,a1,ext,cmp,s0,sxx]=integm(nncl,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

        dmom=amom-2.0*a1;
	
		fyc=dmom;
        if(dmom>0.0)
            fxp=fxc;
            fyp=fyc;
        else
            fxm=fxc;
            fym=fyc;
        end
        	

		if(abs(dmom)<=1.0e28)
            break;
        end
    end

%	call integm(nncl,ndim,st0,space,zn,sgrad,amom,a1)
    [amom,a1,ext,cmp,s0,sxx]=integm(nncl,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx);

if(amom<0)
    amom=amom*(-1);
end
   
tem=(12.0*amom*(1d0-pr*pr)/(abs(INPUT_COL.curv)*em))^(1.0/3.0);
znm=zn;

% C############################################################
% C###### WRITE OUTPUT ########################################
%       
%       open(10, file='yse.lst')

for j=1:nnm
    
	dep=space*double(j-1); %           !depth is in meters
    if(j<nnc) 
        dev=(((rate_c/ac)^(1.0/nc)))*exp(qc/(nc*r*INPUT_COL.rev_T_profile(j)));
        

% C 			If A is iven in Mpa then 

%        dev=dev*1.0e6;
        cmp(j)=(rp-1.0)*INPUT_COL.rev_Litho(j)*(1.0-lambda);         %!compress brittle strength
        ext(j)=((rp-1.0)/rp)*INPUT_COL.rev_Litho(j)*(1.0-lambda);    %!ext brittle strength
        cmp(j)=-min(cmp(j),dev);    % !compressional YSE
        ext(j)=min(ext(j),dev);     % !extensional YSE
	 	sgrad=ec*INPUT_COL.curv; 
	 	sxx(j)=-(sgrad*(dep-znc)+st0);
    elseif((j>nnc) && (j<nncl))
        dev=(((rate_c/acl)^(1.0/ncl)))*exp(qcl/(ncl*r*INPUT_COL.rev_T_profile(j)));

% C 			if A is iven in Mpa then 

  %       dev=dev*1.0e6;

         cmp(j)=(rp-1.0)*INPUT_COL.rev_Litho(j)*(1.0-lambda);         %!compress brittle strength
         ext(j)=((rp-1.0)/rp)*INPUT_COL.rev_Litho(j)*(1.0-lambda);    %!ext brittle strength
         cmp(j)=-min(cmp(j),dev); %     !compressional YSE
         ext(j)=min(ext(j),dev); %      !extensional YSE
	 		
         sgrad=ec*INPUT_COL.curv; 
% c	 		write(*,*) curv
	 	 sxx(j)=-(sgrad*(dep-znc)+st0); 
    else             
	 	dev=(((rate_m/am)^(1.0/nm)))*exp(qm/(nm*r*INPUT_COL.rev_T_profile(j)));

% C			If A is iven in Mpa then          

%	 	dev=dev*1.0e6;
        cmp(j)=(rp-1.0)*INPUT_COL.rev_Litho(j)*(1.0-lambda);         %!compress brittle strength
        ext(j)=((rp-1.0)/rp)*INPUT_COL.rev_Litho(j)*(1.0-lambda);    %!ext brittle strength

        cmp(j)=-min(cmp(j),dev);     
        ext(j)=min(ext(j),dev);
	 	sgrad=ec*INPUT_COL.curv;  
	 	sxx(j)=-(sgrad*(dep-znm)+st0);
    end
		

% c	write(10, '(i4,f8.3,4d13.4)')
% c	1	box,time,dep/1000,cmp(j)/1d6,ext(j)/1d6,sxx(j)/1d6
% 
% c-tiago  	ACCOUNT FOR COUPLING-DECOUPLING EFFECTS       	

	dev_up=(((rate_c/ac)^(1d0/nc)))*exp(qc/(nc*r*INPUT_COL.rev_T_profile(nnc)));
	dev_low=(((rate_c/acl)^(1d0/ncl)))*exp(qcl/(ncl*r*INPUT_COL.rev_T_profile(nncl)));

% 
% C 		if A is given in Mpa then 
        

		
		if((dev_up>20)&&(dev_low>20))
            te_sp=tec+tecl+tem;
        elseif((dev_up>20)&&(dev_low<20))
            te_sp=(((tec+tecl)^3)+tem^3)^(1.0/3.0);
        elseif((dev_up<20)&&(dev_low>20))
            te_sp=(tec^3+((tecl+tem)^3))^(1.0/3.0);
		else 
	 		te_sp=(tec^3+tecl^3+tem^3)^(1.0/3.0);
        end
        
        Tei=te_sp;
end
	
% c      if (box==65) then
% c		write (*,'(1x,a,4f12.3)') "End of TE_YSE for box 65", tec,
% c	1			tecl, tem, te_sp
% c	endif

end




%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

function [amom,a1,ext,cmp,s0,sxx]=integ(ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx)


%	===========================================
%	Calculates the bending moment for the CRUST
%	===========================================

% 	implicit none
% 
% 	integer ndim
% 	double precision st0,space,zn,sgrad,amom,a1
% 
% 	double precision ext(10001),cmp(10001),s0(10001),sxx(10001)
%       common /teyse/ext,cmp,s0,sxx
% 
% 	integer i,j,ii,i1,i2
% 	double precision temp(156),z,z1,z2,zn1,zn2,elmt

% c	START
% c	=====



amom=0.0;
a1=0.0;
j=int16(zn/space)+1;                                                       
      sxx(1)=0.0;
       
      for i=2:ndim     %!where ductile stress is > 0.2 MPa) 
         	z=double(i-1)*space;
         	sxx(i)=sgrad*(z-zn)+st0;
         	if(sxx(i)>0.0)
          		sxx(i)=min(sxx(i),ext(i));
         	else
          		sxx(i)=max(sxx(i),cmp(i));
            end
      end
%c10	continue
      if(s0(j)==st0)
        	for i=2:ndim
         		z1=double(i-2)*space;
         		z2=double(i-1)*space;
         		elmt=0.5*space*((z1-zn)*(sxx(i-1)-s0(i-1))+(z2-zn)*(sxx(i)-s0(i)));
         		amom=amom+elmt; %      !amom stands for area moment
         		if(sxx(i)<s0(i))
                    a1=a1+elmt;
                end
            end
      
		 
%c20     	continue
        
      else

%c		WHEN THERE IS HORIZNTAL STRESSES	
        
		if(sgrad>0.0 && st0>0.0)
            for  i=2:(ndim-1)
         			if(sxx(i) < s0(i))
                        i1=i;
                    end
         			ii=ndim-i+1;
         			if(s0(ii) < ext(ii))
                        i2=ii;
                    end
            end
			 
%c30     		continue
            zn1=space*double(i1);
        	zn2=space*double(i2-1);
        	for i=2:ndim
         		z1=double(i-2)*space;
         		z2=double(i-1)*space;
         		z=0.5*(z1+z2);
         		if(z<zn1)
                    elmt=0.5*space*((z1-zn1)*(sxx(i-1)-s0(i-1))+(z2-zn1)*(sxx(i)-s0(i)));
          			amom=amom+elmt;
          			a1=a1+elmt;
         		else
          			elmt=0.5*space*((z1-zn2)*(sxx(i-1)-s0(i-1))+(z2-zn2)*(sxx(i)-s0(i)));
          			amom=amom+elmt;
                    
                end
            end
%c40     		continue

        end
      end

       
    end			%!End INTEG

%c	ccccccccccccccccccccccccccccccccccccccc
%c	ccccccccccccccccccccccccccccccccccccccc
                           
function [amom,a1,ext,cmp,s0,sxx]=integm(nnc,ndim,st0,space,zn,sgrad,amom,a1,ext,cmp,s0,sxx)

% c	===========================================
% c	Calculates the bending moment for the CRUST
% c	===========================================
% 
% 	implicit none
% 
% 	integer nnc,ndim
% 	double precision st0,space,zn,sgrad,amom,a1
% 
% 	double precision ext(10001),cmp(10001),s0(10001),sxx(10001)
%       common /teyse/ext,cmp,s0,sxx
% 
% 	integer i,j,ii,i1,i2
% 	double precision temp(156),z,z1,z2,zn1,zn2,elmt
% 
% 
% c	START
% c	=====

      amom=0.0;
      a1=0.0;
      j=int16(zn/space)+1;                                                       
      sxx(1)=0.0;
      for i=(nnc+1):ndim     %!where ductile stress is > 0.2 MPa) 
         	z=double(i-1)*space;
         	sxx(i)=sgrad*(z-zn)+st0;
         	if(sxx(i)>0.0)
                sxx(i)=min(sxx(i),ext(i));
            else
          		sxx(i)=max(sxx(i),cmp(i));
            end
      end 
%c10   continue
      if(s0(j)==st0)
          for i=(nnc+1):ndim
         		z1=double(i-2)*space;
         		z2=double(i-1)*space;
         		elmt=0.5*space*((z1-zn)*(sxx(i-1)-s0(i-1))+(z2-zn)*(sxx(i)-s0(i)));
         		amom=amom+elmt;    %  !amom stands for area moment
                if(sxx(i)<s0(i))
                    a1=a1+elmt;
                end
          end
       
%c20     	continue
     else

%c		WHEN THERE IS HORIZONTAL STRESSES	
       
    if(sgrad>0.0 && st0>0.0)
        for i=(nnc+1):(ndim-1)
            if(sxx(i)<s0(i))
                i1=i;
            end
            ii=ndim-i+1;
         	if(s0(ii)<ext(ii))
                i2=ii;
            end
        end 
%c30     		continue
        zn1=space*double(i1);
        zn2=space*double(i2-1);
        for i=(nnc+1):ndim
            z1=double(i-2)*space;
            z2=double(i-1)*space;
         	z=0.5*(z1+z2);
         	if(z<zn1)
                elmt=0.5*space*((z1-zn1)*(sxx(i-1)-s0(i-1))+(z2-zn1)*(sxx(i)-s0(i)));
          		amom=amom+elmt;
          		a1=a1+elmt;
            else
                elmt=0.5*space*((z1-zn2)*(sxx(i-1)-s0(i-1))+(z2-zn2)*(sxx(i)-s0(i)));
          		amom=amom+elmt;
            end
        end
%c40     		continue

    end
	
      end
end %!End INTEGM
