function G = mesh_study(f_coarse, f_medium, f_fine, r)
%MESH_STUDY  Grid convergence from three meshes (Richardson extrapolation + GCI).
%
%   G = mesh_study(f_coarse, f_medium, f_fine)        refinement ratio r = 2
%   G = mesh_study(f_coarse, f_medium, f_fine, r)
%
%   The three values are one result (e.g. CL at 10 deg) on meshes whose cell
%   size halves at each step. Method: Roache's Grid Convergence Index, as set
%   out by Celik et al., J. Fluids Eng. 130 (2008) 078001.
%
%   G.type     'monotonic'   differences shrink in the same direction -> extrapolate
%              'oscillatory' differences change sign -> uncertainty = half the range
%              'diverging'   differences grow -> no mesh-independent value yet
%              'identical'   all three values equal -> uncertainty zero
%   G.p        observed order of accuracy (about 2 expected for these solvers)
%   G.f_ext    extrapolated zero-cell-size value (monotonic only)
%   G.U        absolute numerical uncertainty of each mesh [coarse medium fine]
%   G.Fs       safety factor used: 1.25 normally; 3 (with p = 2) if the observed
%              order is outside 0.5..4, where the extrapolation is not trusted

    if nargin < 4, r = 2; end
    f = [f_coarse f_medium f_fine];
    e21 = f_fine - f_medium;  e32 = f_medium - f_coarse;
    tol = 1e-9*max(abs(f));                       % differences below round-off count as zero
    if abs(e21) <= tol, e21 = 0; end
    if abs(e32) <= tol, e32 = 0; end
    G.f = f; G.r = r; G.p = NaN; G.f_ext = NaN; G.Fs = NaN; G.U = NaN(1, 3);
    if e21 == 0 && e32 == 0
        G.type = 'identical'; G.f_ext = f_fine; G.Fs = 1.25; G.U = zeros(1, 3);
        return
    end
    R = e21/e32;
    if R < 0
        G.type = 'oscillatory';
        G.U(:) = 0.5*(max(f) - min(f));
    elseif R >= 1
        G.type = 'diverging';
    else
        G.type = 'monotonic';
        p = log(1/R)/log(r);
        Fs = 1.25; p_used = p;
        if p < 0.5 || p > 4, Fs = 3; p_used = 2; end
        G.p = p; G.Fs = Fs;
        G.f_ext = f_fine + e21/(r^p_used - 1);
        G.U = Fs*abs(f - G.f_ext);
    end
end
