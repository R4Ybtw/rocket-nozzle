function T = xyplot_table(sheets)
%XYPLOT_TABLE  Columns of a Flow Simulation XY-plot export, whatever its layout.
%
%   T = xyplot_table(sheets)     sheets: cell array of worksheets, each a cell
%                                array as readcell returns it
%
%   In every sheet the header row is the first row holding an abscissa cell
%   ("Model X", "Model Y", ...) and at least one parameter. Parameters found:
%     p    "Pressure" (static; not total or dynamic)       Pa
%     rho  "Density"                                        kg/m^3
%     u    "Velocity (X)" or "Velocity X"                   m/s
%     M    "Mach Number"                                    -
%   T.s is the abscissa (m) of the first parameter found; every other
%   parameter is interpolated onto it, so one table results whether the export
%   put all parameters on one sheet or one per sheet. Units in brackets are
%   honoured for lengths (m, mm) and pressures (Pa, kPa, MPa).

    keys = {'p', '^\s*(static\s+)?pressure(\s|\[|\(|$)'; 'rho', 'density'; 'u', 'velocity\s*\(?x\)?'; 'M', 'mach'};
    T = struct('s', []);
    found = struct();
    for w = 1:numel(sheets)
        C = sheets{w};
        [nr, nc] = size(C);
        for r = 1:nr
            ia = 0;
            for c = 1:nc
                if istext(C{r,c}) && ~isempty(regexpi(C{r,c}, 'model\s*[xyz]', 'once')), ia = c; break; end
            end
            if ~ia, continue; end
            any_param = false;
            for c = 1:nc
                if c == ia || ~istext(C{r,c}), continue; end
                for k = 1:size(keys, 1)
                    if isfield(found, keys{k,1}), continue; end
                    if ~isempty(regexpi(C{r,c}, keys{k,2}, 'once')) && isempty(regexpi(C{r,c}, 'total|dynamic', 'once'))
                        [s, v] = column_pair(C, r, ia, c);
                        if numel(s) >= 2
                            found.(keys{k,1}) = struct('s', s*length_scale(C{r,ia}), 'v', v*value_scale(C{r,c}));
                            any_param = true;
                        end
                    end
                end
            end
            if any_param, break; end
        end
    end
    f = fieldnames(found);
    if isempty(f), return; end
    T.s = found.(f{1}).s;
    for k = 1:numel(f)
        F = found.(f{k});
        [ss, i] = unique(F.s);
        T.(f{k}) = interp1(ss, F.v(i), T.s, 'linear', NaN);
    end
end

function tf = istext(v)
    tf = ischar(v) && ~isempty(v);
end

function [s, v] = column_pair(C, r, ia, c)
    s = []; v = [];
    num = @(x) isnumeric(x) && isscalar(x) && isfinite(x);
    for q = r+1:size(C, 1)
        if num(C{q,ia}) && num(C{q,c})
            s(end+1, 1) = C{q,ia}; v(end+1, 1) = C{q,c}; %#ok<AGROW>
        end
    end
end

function k = length_scale(h)
    k = 1; if ~isempty(regexpi(h, '\[\s*mm\s*\]', 'once')), k = 1e-3; end
end

function k = value_scale(h)
    k = 1;
    if ~isempty(regexpi(h, '\[\s*kPa\s*\]', 'once')), k = 1e3; end
    if ~isempty(regexpi(h, '\[\s*MPa\s*\]', 'once')), k = 1e6; end
end
