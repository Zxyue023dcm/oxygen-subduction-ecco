function idx = find_peak_plateau(x)
%FIND_PEAK_PLATEAU Return all indices belonging to local peak plateaus.
%   A plateau is returned when its adjacent values are both lower. End
%   points are not classified as peaks.

x=x(:);
n=numel(x);
if n<3
    idx=[];
    return
end

change_pos=find(diff(x)~=0);
starts=[1;change_pos+1];
ends=[change_pos;n];
isPeak=false(n,1);

for k=1:numel(starts)
    left=starts(k);
    right=ends(k);
    if left>1 && right<n && ...
            x(left-1)<x(left) && x(right+1)<x(right)
        isPeak(left:right)=true;
    end
end

idx=find(isPeak);
end
