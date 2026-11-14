function [L, xpk] = cell_spacing_from_axis(x, p, pa, x0, Lguess)
%CELL_SPACING_FROM_AXIS  Shock-cell length from the pressure along the jet axis.
%
%   [L, xpk] = cell_spacing_from_axis(x, p, pa, x0, Lguess)
%
%   Pressure on the axis of an under-expanded jet rises and falls once per
%   shock cell. The peaks are the local maxima downstream of x0 (the exit) that
%   exceed ambient pressure, kept at least 0.3 Lguess apart; L is the mean
%   distance between successive peaks (NaN if fewer than two).

    x = x(:); p = p(:);
    s = x > x0 & isfinite(p);
    x = x(s); p = p(s);
    xpk = [];
    for k = 2:numel(p) - 1
        if p(k) >= p(k-1) && p(k) > p(k+1) && p(k) > pa
            if isempty(xpk) || x(k) - xpk(end) > 0.3*Lguess
                xpk(end+1) = x(k); %#ok<AGROW>
            end
        end
    end
    L = NaN;
    if numel(xpk) >= 2, L = mean(diff(xpk)); end
end
