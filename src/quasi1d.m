function Q = quasi1d(ratio_e, pb_p0, g)
%QUASI1D  Exit state of a choked convergent-divergent nozzle, quasi-1D theory.
%
%   Q = quasi1d(ratio_e, pb_p0, g)
%
%   ratio_e   exit-to-throat area ratio (height ratio for a planar nozzle)
%   pb_p0     back (ambient) pressure over chamber stagnation pressure
%
%   Regimes, by falling back pressure:
%     'subsonic'      not choked at the throat (outside this project's cases)
%     'shock inside'  normal shock in the divergent part, subsonic exit at pb
%     'overexpanded'  supersonic exit, pe < pb: oblique shocks outside (inviscid)
%     'design'        pe = pb (within 1%)
%     'underexpanded' pe > pb: expansion fans outside
%   Q.regime, Q.Me, Q.pe_p0, Q.A_shock (shock area / throat area, NaN if none),
%   Q.pb_limits = [pb/p0 at choking, with a shock at the exit, at design].

    Mdes = gasdyn('M_from_A', ratio_e, g, 'sup');
    Msub = gasdyn('M_from_A', ratio_e, g, 'sub');
    p_des = 1/gasdyn('p0_p', Mdes, g);
    p_sub = 1/gasdyn('p0_p', Msub, g);
    S = gasdyn('normal_shock', Mdes, g);
    p_exit_shock = p_des*S.p2_p1;
    Q.pb_limits = [p_sub p_exit_shock p_des];
    Q.A_shock = NaN;
    if pb_p0 >= p_sub
        Q.regime = 'subsonic'; Q.Me = gasdyn('M_from_p', 1/pb_p0, g); Q.pe_p0 = pb_p0;
    elseif pb_p0 > p_exit_shock
        % shock at area ratio As: total pressure falls by p02/p01, the sonic
        % area grows to A*/(p02/p01), and the subsonic exit must reach pb
        f = @(As) shock_exit_p(As, ratio_e, g) - pb_p0;
        As = bisect(f, 1 + 1e-9, ratio_e);
        Q.regime = 'shock inside'; Q.A_shock = As;
        [Q.pe_p0, Q.Me] = shock_exit_p(As, ratio_e, g);
    else
        Q.Me = Mdes; Q.pe_p0 = p_des;
        if abs(pb_p0/p_des - 1) < 0.01, Q.regime = 'design';                  % within 1%
        elseif pb_p0 > p_des, Q.regime = 'overexpanded';
        else, Q.regime = 'underexpanded'; end
    end
end

function [pe, Me] = shock_exit_p(As, ratio_e, g)
    M1 = gasdyn('M_from_A', As, g, 'sup');
    S = gasdyn('normal_shock', M1, g);
    Me = gasdyn('M_from_A', ratio_e*S.p02_p01, g, 'sub');
    pe = S.p02_p01/gasdyn('p0_p', Me, g);
end

function x = bisect(f, lo, hi)
    flo = f(lo);
    for it = 1:200
        x = 0.5*(lo + hi); fx = f(x);
        if (hi - lo) < 1e-13, return; end
        if sign(fx) == sign(flo), lo = x; flo = fx; else, hi = x; end
    end
end
