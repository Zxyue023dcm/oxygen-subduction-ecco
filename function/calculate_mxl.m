function Hml = calculate_mxl(theta,salt,dep,lon,lat)
%CALCULATE_MXL Calculate mixed-layer depth using a sigma0 threshold of 0.03.
%   THETA and SALT are ordered as longitude-by-latitude-by-depth-by-time.

mxl=nan(size(theta,1),size(theta,2),size(theta,4));
theta(theta==0)=nan;
salt(salt==0)=nan;

% Calculate potential density anomaly.
depth=repmat(reshape(dep,[1 1 size(theta,3) 1]), ...
    [size(theta,1) size(theta,2) 1 size(theta,4)]);
LAT=repmat(lat,[1 1 size(theta,3) size(theta,4)]);
LON=repmat(lon,[1 1 size(theta,3) size(theta,4)]);
p=gsw_p_from_z(depth,LAT);
SA=gsw_SA_from_SP(salt,p,LON,LAT);
CT=gsw_CT_from_pt(SA,theta);
sigma=gsw_sigma0(SA,CT);
clear p SA CT theta salt depth LON LAT

% Locate the 0.03 kg m^-3 density increase at each grid point.
for m=1:size(sigma,4)
    disp(m);tic
    for i=1:size(sigma,1)
        for j=1:size(sigma,2)
            SIG=squeeze(sigma(i,j,:,m));
            if isnan(SIG(1))
                continue
            end

            valid=isfinite(SIG);
            SIG=SIG(valid);
            dep1=dep(valid);
            if numel(SIG)<2
                continue
            end

            target=SIG(1)+0.03;
            if issorted(SIG,'ascend')
                mxl(i,j,m)=interp1(SIG,-dep1,target,'linear');
            elseif SIG(2)>=SIG(1)
                idx=find_peak_plateau(SIG);
                if ~isempty(idx) && SIG(idx(1))>target
                    mxl(i,j,m)=interp1(SIG(1:idx(1)), ...
                        -dep1(1:idx(1)),target,'linear');
                else
                    idx=find_valley_plateau(SIG);
                    if ~isempty(idx)
                        mxl(i,j,m)=interp1(SIG(idx(end):end), ...
                            -dep1(idx(end):end),target,'linear');
                    end
                end
            else
                idx=find_valley_plateau(SIG);
                if ~isempty(idx)
                    mxl(i,j,m)=interp1(SIG(idx(end):end), ...
                        -dep1(idx(end):end),target,'linear');
                end
            end
        end
    end
    toc
end

Hml=mxl;
end
