function mean_val = weighted_average_depth(depth,value)
%WEIGHTED_AVERAGE_DEPTH Calculate a layer-thickness-weighted mean.
%   DEPTH contains cell-center depths and may be unevenly spaced. Layer
%   boundaries are the midpoints between adjacent depths.

if ~isvector(depth) || ~isvector(value)
    error('weighted_average_depth:InvalidInput', ...
        'Depth and value must both be vectors.');
end

depth=depth(:);
value=value(:);
if numel(depth)~=numel(value)
    error('weighted_average_depth:SizeMismatch', ...
        'Depth and value must have the same length.');
end

n=numel(depth);
if n==1
    mean_val=value;
    return
end

[depth_sorted,idx]=sort(depth);
value_sorted=value(idx);
thickness=zeros(n,1);

thickness(1)=(depth_sorted(2)-depth_sorted(1))/2;
for i=2:n-1
    upper=(depth_sorted(i-1)+depth_sorted(i))/2;
    lower=(depth_sorted(i)+depth_sorted(i+1))/2;
    thickness(i)=lower-upper;
end
thickness(n)=(depth_sorted(n)-depth_sorted(n-1))/2;

mean_val=sum(value_sorted.*thickness)/sum(thickness);
end
