function R = read_runs(file)
%READ_RUNS  data/results/runs.csv: case,mesh,mdot_inlet,span_m,converged
%   mdot_inlet  mass flow through the inlet lid as Flow Simulation reports it (kg/s)
%   span_m      Z max - Z min of the 2D domain (m); blank rows inherit the nearest given
%   converged   1 or 0 (blank = 1). Rows without a mass flow are ignored.
    R = struct('case', {{}}, 'mesh', {{}}, 'mdot', [], 'span', [], 'converged', [], 'n', 0);
    if ~exist(file, 'file'), return; end
    L = regexp(fileread(file), '\r?\n', 'split');
    for k = 2:numel(L)
        q = regexp(strtrim(L{k}), ',', 'split');
        if numel(q) < 3, continue; end
        m = str2double(q{3}); if isnan(m), continue; end
        q(end+1:5) = {''};
        R.case{end+1, 1} = strtrim(q{1}); R.mesh{end+1, 1} = strtrim(q{2});
        R.mdot(end+1, 1) = abs(m); R.span(end+1, 1) = str2double(q{4});
        c = str2double(q{5}); if isnan(c), c = 1; end
        R.converged(end+1, 1) = logical(c);
    end
    R.n = numel(R.mdot);
    given = find(isfinite(R.span));
    for k = 1:R.n
        if isnan(R.span(k)) && ~isempty(given), [~, j] = min(abs(given - k)); R.span(k) = R.span(given(j)); end
    end
end
