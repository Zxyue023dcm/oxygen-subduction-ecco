function Wek = cal_Wek(oceanlon,oceanlat,taue,taun)
%CAL_WEK Calculate Ekman pumping velocity from wind stress.

[lon_l,lat_l]=distcoordinate(oceanlon(:,1),oceanlat(1,:)');
f=sw_f(oceanlat);
rou0=1026; % kg m^-3
Tx=taue./(f*rou0);
Ty=taun./(f*rou0);
Wek=nan(size(taue));

for t=1:size(taue,3)
    [~,gXy,gYx]=curldffs_region(Tx(:,:,t),lon_l,Ty(:,:,t),lat_l);
    Wek(:,:,t)=gYx-gXy; % m s^-1
end
end

