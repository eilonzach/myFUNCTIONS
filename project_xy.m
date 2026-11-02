function [result1, result2] = project_xy_v2(par, input1, input2, direction)
%PROJECT_XY_V2 Convert between geographic coordinates and local x-y.
%
%   [x, y] = project_xy_v2(par, lat, lon, 'forward')
%   [lat, lon] = project_xy_v2(par, x, y, 'inverse')
%
%   The projection method depends on par.map_proj:
%     - 'naive'   : simple local approximation in km
%     - 'flat'    : flat-earth style approximation in km
%     - map object : MATLAB projection forward/inverse functions

narginchk(3, 4)

% Default to forward projection if direction is omitted.
if nargin == 3
    direction = 'forward';
elseif ~(ischar(direction) || isstring(direction)) || ...
        (~strcmpi(direction, 'forward') && ~strcmpi(direction, 'inverse'))
    error('incorrect 4th argument');
end

switch lower(char(direction))

    case 'forward'
        % Convert latitude/longitude to local x-y coordinates.
        lat = input1;
        lon = input2;

        if isnaive(par)
            % Naive local tangent-plane approximation.
            y = deg2km(lat - par.origin(1));
            x = deg2km(lon - par.origin(2)) * cos(deg2rad(par.origin(1)));

        elseif isflat(par)
            % Flat-earth approximation in km.
            y = deg2km(lat - par.origin(1));
            x = deg2km(lon - par.origin(2));

        elseif ismaptoolbox(par)
            % Project geographic coordinates to planar coordinates.
            % projfwd returns meters, so convert to km for consistency.
            [xMeters, yMeters] = projfwd(par.map_proj, lat, lon);
            x = xMeters / 1000;
            y = yMeters / 1000;

        else
            error('invalid choice for map_proj');
        end

        result1 = x;
        result2 = y;

    case 'inverse'
        % Convert local x-y coordinates back to latitude/longitude.
        x = input1;
        y = input2;

        if isnaive(par)
            error('haven''t implemented naive inverse projection');

        elseif isflat(par)
            % Reverse the flat-earth approximation.
            lat = km2deg(y) + par.origin(1);
            lon = km2deg(x) + par.origin(2);

        elseif ismaptoolbox(par)
            % projinv expects meters, so convert from km first.
            xMeters = x * 1000;
            yMeters = y * 1000;
            [lat, lon] = projinv(par.map_proj, xMeters, yMeters);

        else
            error('invalid choice for map_proj');
        end

        result1 = lat;
        result2 = lon;

    otherwise
        error('direction must be forwards or backwards');
end

end

function truefalse = isflat(par)
%ISFLAT True when par.map_proj requests the flat-earth approximation.
truefalse = ~isstruct(par.map_proj) && strcmpi('flat', par.map_proj);
end

function truefalse = isnaive(par)
%ISNAIVE True when par.map_proj requests the naive local approximation.
truefalse = ~isstruct(par.map_proj) && strcmpi('naive', par.map_proj);
end

function truefalse = ismaptoolbox(par)
%ISMAPTOOLBOX True when par.map_proj is a MATLAB map/projection object.
truefalse = false;

if isstruct(par.map_proj) && isfield(par.map_proj, 'mapprojection')
    truefalse = true;
elseif exist('projcrs', 'class') == 8 && isa(par.map_proj, 'projcrs')
    truefalse = true;
end
end
