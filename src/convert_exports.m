function n = convert_exports(root)
%CONVERT_EXPORTS  Turn the XY-plot exports into plain CSV profiles (MATLAB only).
%
%   n = convert_exports(root)
%
%   Reads data/flowsim_exports/<case>_<mesh>_<line>.xlsx (or .csv), where
%     case  design | over | under      mesh  coarse | medium | fine
%     line  exit  (across the exit plane: Pressure, Density, Velocity (X) against Model Y)
%           wall  (along the inner wall: Pressure against Model X)
%           axis  (along the axis behind the exit: Pressure against Model X)
%   and writes data/results/<same name>.csv with columns s,p,rho,u,M
%   (missing parameters as NaN). Returns the number converted.

    n = 0;
    if exist('OCTAVE_VERSION', 'builtin'), return; end
    d = [dir(fullfile(root, 'data', 'flowsim_exports', '*.xlsx')); dir(fullfile(root, 'data', 'flowsim_exports', '*.csv'))];
    for k = 1:numel(d)
        tok = regexp(d(k).name, '^(design|over|under)_(coarse|medium|fine)_(exit|wall|axis)\.(xlsx|csv)$', 'tokens', 'once');
        if isempty(tok), fprintf('  skipped %s\n', d(k).name); continue; end
        f = fullfile(d(k).folder, d(k).name);
        if strcmpi(tok{4}, 'csv'), sh = {''}; else, sh = cellstr(sheetnames(f)); end
        sheets = {};
        for w = 1:numel(sh)
            try
                if isempty(sh{w}), C = readcell(f); else, C = readcell(f, 'Sheet', sh{w}); end
            catch
                continue
            end
            C(cellfun(@(v) isa(v, 'missing'), C)) = {[]};
            sheets{end+1} = cellfun(@(v) convertStringsToChars(v), C, 'UniformOutput', false); %#ok<AGROW>
        end
        T = xyplot_table(sheets);
        if isempty(T.s), fprintf('  %s: no Model X/Y column with a parameter beside it\n', d(k).name); continue; end
        cols = {'p', 'rho', 'u', 'M'}; D = T.s;
        for c = 1:4
            if isfield(T, cols{c}), D(:, c+1) = T.(cols{c}); else, D(:, c+1) = NaN; end
        end
        out = fullfile(root, 'data', 'results', sprintf('%s_%s_%s.csv', tok{1:3}));
        fid = fopen(out, 'w'); fprintf(fid, 's,p,rho,u,M\n'); fprintf(fid, '%.8g,%.8g,%.8g,%.8g,%.8g\n', D'); fclose(fid);
        n = n + 1;
    end
end
