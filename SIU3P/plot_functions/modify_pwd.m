function [ISOCpwd_new] = modify_pwd(ISOCHRONS,Topography)

%

[ISOCx,ISOCy_ini,ISOCyy,ISOCpwd,ero_bool]=find_ero(ISOCHRONS,Topography);

ISOCpwd_new = ISOCpwd;

ero_bool = (ISOCy_ini~=ISOCyy);

NPoints = size(ISOCpwd(1,:),2);

isoc_i = unique(ISOCHRONS(3,:));

for k = 2:size(isoc_i,2)
    
    for j = 1:NPoints
       
        if ero_bool(k,j) == 1
            
            py     = ISOCyy(k,j);
            py_pre = ISOCyy(k-1,j);
            py_ero = ISOCy_ini(k,j);
            
            pwd_pre = ISOCpwd(k-1,j);
            pwd_ini = ISOCpwd(k,j);
          
            if abs(py_pre-py_ero)> 0 
                
                pwd     = interp1([py_pre py_ero],[pwd_pre pwd_ini],py);
                
                ISOCpwd_new(k,j) = pwd; 
            end
            
        end
        
    end
    
end