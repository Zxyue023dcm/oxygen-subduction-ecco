function S = positiveCentroidBands(Q, lat, lon, pacificLon, power, negFactor, oceanMask, latBands, bandNames, qInfluence)
% positiveCentroidBands
%
% Signed-contrast centroid:
%   alpha = negFactor
%   alpha = 0   -> original positive-only centroid
%   alpha > 0   -> negative values contribute oppositely
%
% Q          : [nlat, nlon, nt]
% lat        : latitude vector, degrees
% lon        : longitude vector, degrees
% pacificLon : [lon1 lon2], e.g. [120 290] for 120E-70W
% power      : exponent for |Q| weighting, usually 1 or 2
% negFactor  : alpha, negative-value contribution factor, usually 0~1
% oceanMask  : optional [nlat, nlon] logical mask
% latBands   : [nBand,2], each row = [latMin latMax]
% bandNames  : optional cellstr
% qInfluence : 0~1, degree of Q-amplitude participation in weights
%              1 -> original amplitude-weighted scheme
%              0 -> ignore Q magnitude and use support-only area weights
%
% Main outputs:
% S.lon360, S.lon180, S.latSph
%   signed-contrast centroid on the sphere
%
% Diagnostics:
% S.latMeanSigned
% S.posAbsLatMean, S.negAbsLatMean
% S.netWeight, S.magWeight, S.coherence

if nargin < 4 || isempty(pacificLon)
    pacificLon = [120 290];
end

if nargin < 5 || isempty(power)
    power = 1;
end

if nargin < 6 || isempty(negFactor)
    negFactor = 1;
end

if nargin < 7 || isempty(oceanMask)
    oceanMask = [];
end

if nargin < 8 || isempty(latBands)
    latBands = [-90 90];
end

if nargin < 10 || isempty(qInfluence)
    qInfluence = 1;
end

if ndims(Q) == 2
    Q = reshape(Q, size(Q,1), size(Q,2), 1);
end

[nlat, nlon, nt] = size(Q);

lat = lat(:);
lon = lon(:);

if length(lat) ~= nlat
    error('Length of lat does not match size(Q,1).');
end

if length(lon) ~= nlon
    error('Length of lon does not match size(Q,2).');
end

if size(latBands,2) ~= 2
    error('latBands must be [nBand,2], each row = [latMin latMax].');
end

if power < 1
    error('power must be >= 1.');
end

if negFactor < 0
    error('negFactor must be >= 0.');
end

if qInfluence < 0 || qInfluence > 1
    error('qInfluence must be within [0, 1].');
end

nBand = size(latBands,1);

for k = 1:nBand
    if latBands(k,1) > latBands(k,2)
        error('Each row of latBands must satisfy latMin <= latMax.');
    end
    if latBands(k,1) < -90 || latBands(k,2) > 90
        error('latBands must stay within [-90, 90].');
    end
end

if nargin < 9 || isempty(bandNames)
    bandNames = cell(nBand,1);
    for k = 1:nBand
        bandNames{k} = sprintf('Band_%d_[%g,%g]', k, latBands(k,1), latBands(k,2));
    end
end

if numel(bandNames) ~= nBand
    error('Length of bandNames must equal number of rows in latBands.');
end

[lat, ilat] = sort(lat, 'ascend');

lon360 = mod(lon, 360);
[lon360, ilon] = sort(lon360, 'ascend');

Q = Q(ilat, ilon, :);

if ~isempty(oceanMask)
    oceanMask = oceanMask(ilat, ilon);
else
    oceanMask = true(nlat, nlon);
end

[Lon, Lat] = meshgrid(lon360, lat);

latEdge = centerToEdge(lat);
latEdge(1)   = max(latEdge(1), -90);
latEdge(end) = min(latEdge(end), 90);

lonEdge = centerToEdge(lon360);

dSinLat = sind(latEdge(2:end)) - sind(latEdge(1:end-1));
dLon    = deg2rad(diff(lonEdge));
cellArea = dSinLat(:) * dLon(:).';

lonRange = mod(pacificLon, 360);

if lonRange(1) <= lonRange(2)
    pacificLonMask = (Lon >= lonRange(1)) & (Lon <= lonRange(2));
else
    pacificLonMask = (Lon >= lonRange(1)) | (Lon <= lonRange(2));
end

baseMask = pacificLonMask & oceanMask;

bandMask = cell(nBand,1);
for k = 1:nBand
    bandMask{k} = (Lat >= latBands(k,1)) & (Lat <= latBands(k,2));
end

S.lon360         = nan(nt, nBand);
S.lon180         = nan(nt, nBand);
S.latSph         = nan(nt, nBand);
S.latMeanSigned  = nan(nt, nBand);

S.posAbsLatMean  = nan(nt, nBand);
S.negAbsLatMean  = nan(nt, nBand);
S.posLatMean     = nan(nt, nBand);
S.negLatMean     = nan(nt, nBand);

S.netWeight      = nan(nt, nBand);
S.magWeight      = nan(nt, nBand);
S.coherence      = nan(nt, nBand);

S.bandNames      = bandNames(:);
S.latBands       = latBands;
S.power          = power;
S.negFactor      = negFactor;
S.qInfluence     = qInfluence;

for it = 1:nt
    q = Q(:,:,it);

    qPos = max(q, 0).^power;
    qNeg = max(-q, 0).^power;

    bad = ~isfinite(q);
    qPos(bad) = 0;
    qNeg(bad) = 0;

    % Blend between amplitude weighting and support-only weighting.
    % qInfluence = 1 keeps the original Q-based weights.
    % qInfluence = 0 uses only the presence of positive/negative values.
    qPosEff = qInfluence * qPos + (1 - qInfluence) * double(qPos > 0);
    qNegEff = qInfluence * qNeg + (1 - qInfluence) * double(qNeg > 0);

    for k = 1:nBand
        mask = baseMask & bandMask{k};

        wPos = qPosEff .* cellArea .* double(mask);
        wNeg = negFactor * qNegEff .* cellArea .* double(mask);

        wSigned = wPos - wNeg;
        wMag    = wPos + wNeg;

        posTotal = sum(wPos(:));
        negTotal = sum(wNeg(:));
        magTotal = sum(wMag(:));
        netTotal = sum(wSigned(:));

        S.netWeight(it,k) = netTotal;
        S.magWeight(it,k) = magTotal;

        if posTotal > 0
            S.posLatMean(it,k)    = sum(wPos(:) .* Lat(:)) / posTotal;
            S.posAbsLatMean(it,k) = sum(wPos(:) .* abs(Lat(:))) / posTotal;
        end

        if negTotal > 0
            S.negLatMean(it,k)    = sum(wNeg(:) .* Lat(:)) / negTotal;
            S.negAbsLatMean(it,k) = sum(wNeg(:) .* abs(Lat(:))) / negTotal;
        end

        if magTotal <= 0
            continue;
        end

        x = sum(wSigned(:) .* cosd(Lat(:)) .* cosd(Lon(:)));
        y = sum(wSigned(:) .* cosd(Lat(:)) .* sind(Lon(:)));
        z = sum(wSigned(:) .* sind(Lat(:)));

        vecNorm = hypot(hypot(x, y), z);
        S.coherence(it,k) = vecNorm / magTotal;

        if vecNorm <= 0
            continue;
        end

        lonC = atan2d(y, x);
        latC = atan2d(z, hypot(x, y));

        S.lon360(it,k) = mod(lonC, 360);
        S.lon180(it,k) = mod(lonC + 180, 360) - 180;
        S.latSph(it,k) = latC;

        S.latMeanSigned(it,k) = sum(wSigned(:) .* Lat(:)) / magTotal;
    end
end

end


function edge = centerToEdge(center)

center = center(:);

if length(center) == 1
    edge = [center - 0.5; center + 0.5];
    return;
end

dc = diff(center);

edge = [
    center(1) - dc(1)/2
    center(1:end-1) + dc/2
    center(end) + dc(end)/2
];

end

