function out = rpt(out, fmt, varargin)
%RPT  Print a line of the analysis report and keep it for summary.md.
    s = sprintf(fmt, varargin{:});
    fprintf('%s\n', s);
    out{end+1} = s;
end
