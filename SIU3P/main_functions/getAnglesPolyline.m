function angles = getAnglesPolyline(xy)
   % +++ purpose +++
   % get (minimum) angles [in (0,180)] between consecutive segments along a polyline
   
   diffx = diff(xy(1,:));
   diffy = diff(xy(2,:));
   alphas = atand(diffy./diffx);
   plus180 = (diffx < 0 & diffy < 0) | (diffx < 0 & diffy > 0);
   alphas(plus180) = alphas(plus180) + 180;
   plus360 = alphas < 0;
   alphas(plus360) = alphas(plus360) + 360;                                % text(xy(1,1:end-1), xy(2,1:end-1),string(alphas))
   angles = abs(180 - abs(diff(alphas)));                                       

   % test
   if 1 > 2 % not run
       xy = rand(2,10);
       figure(); plot(xy(1,:),xy(2,:),'.-','markersize',10)
       angles = getAnglesPolyline(xy);
       text(xy(1,2:end-1), xy(2,2:end-1)+0.01,string(angles),'color','red') 
   end
end % function
