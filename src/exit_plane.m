function E = exit_plane(P, pa, pc, hs, g)
%EXIT_PLANE  Mass flow and thrust from the CFD profile across the exit plane.
%
%   E = exit_plane(P, pa, pc, hs, g)
%
%   P   profile across the exit of the upper half (s = y from the axis), with
%       pressure, density and axial velocity
%   Per unit depth, for the WHOLE nozzle (twice the half):
%     mdot = 2 int rho u dy
%     F    = 2 int (rho u^2 + p - pa) dy        (momentum + pressure thrust)
%     C_F  = F / (pc 2 hs)
%   E.M_mean is the mass-flow-weighted exit Mach number, from u / sqrt(g p/rho).

    y = P.s; p = P.p; r = P.rho; u = P.u;
    ok = isfinite(y) & isfinite(p) & isfinite(r) & isfinite(u);
    y = y(ok); p = p(ok); r = r(ok); u = u(ok);
    E.mdot = 2*trapz(y, r.*u);
    E.F = 2*trapz(y, r.*u.^2 + p - pa);
    E.CF = E.F/(pc*2*hs);
    M = u./sqrt(g*p./r);
    E.M_mean = trapz(y, r.*u.*M)/trapz(y, r.*u);
end
