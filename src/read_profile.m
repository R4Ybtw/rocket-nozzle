function P = read_profile(file)
%READ_PROFILE  A converted XY-plot profile: data/results/<case>_<mesh>_<line>.csv
%   P.s, P.p, P.rho, P.u, P.M (column vectors, sorted by s); P.n = 0 if absent.
    P = struct('s', [], 'p', [], 'rho', [], 'u', [], 'M', [], 'n', 0);
    if ~exist(file, 'file'), return; end
    L = regexp(fileread(file), '\r?\n', 'split'); D = [];
    for k = 2:numel(L)
        v = str2double(regexp(L{k}, ',', 'split'));
        if numel(v) >= 2 && ~isnan(v(1)), D(end+1, 1:numel(v)) = v; end %#ok<AGROW>
    end
    if isempty(D), return; end
    D(:, end+1:5) = NaN;
    [~, i] = sort(D(:,1)); D = D(i,:);
    P.s = D(:,1); P.p = D(:,2); P.rho = D(:,3); P.u = D(:,4); P.M = D(:,5); P.n = size(D, 1);
end
