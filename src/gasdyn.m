function out = gasdyn(fn, varargin)
%GASDYN  Perfect-gas relations used throughout the project (gamma passed explicitly).
%
%   p0_p   = gasdyn('p0_p',   M, g)        stagnation over static pressure
%   T0_T   = gasdyn('T0_T',   M, g)
%   A_As   = gasdyn('A_As',   M, g)        area over sonic area
%   M      = gasdyn('M_from_A', A_As, g, 'sup' | 'sub')
%   M      = gasdyn('M_from_p', p0_p, g)
%   nu     = gasdyn('nu', M, g)            Prandtl-Meyer angle, radians
%   M      = gasdyn('M_from_nu', nu, g)
%   mu     = gasdyn('mu', M)               Mach angle, radians
%   S      = gasdyn('normal_shock', M1, g) S.M2, S.p2_p1, S.T2_T1, S.rho2_rho1, S.p02_p01
%
%   Isentropic flow, normal shocks and Prandtl-Meyer expansion (Anderson,
%   Modern Compressible Flow, chapters 3-5). Inverse functions use bracketed
%   bisection, so they always converge and need no starting guess.

    switch fn
        case 'p0_p'
            [M, g] = varargin{:}; out = (1 + (g-1)/2*M.^2).^(g/(g-1));
        case 'T0_T'
            [M, g] = varargin{:}; out = 1 + (g-1)/2*M.^2;
        case 'A_As'
            [M, g] = varargin{:};
            out = 1./M.*((2/(g+1))*(1 + (g-1)/2*M.^2)).^((g+1)/(2*(g-1)));
        case 'M_from_A'
            [a, g, branch] = varargin{:};
            out = zeros(size(a));
            for k = 1:numel(a)
                if strcmp(branch, 'sup'), lo = 1; hi = 50; else, lo = 1e-8; hi = 1; end
                f = @(M) gasdyn('A_As', M, g) - a(k);
                out(k) = bisect(f, lo, hi);
            end
        case 'M_from_p'
            [r, g] = varargin{:};
            out = sqrt(2/(g-1)*(r.^((g-1)/g) - 1));
        case 'nu'
            [M, g] = varargin{:};
            c = sqrt((g+1)/(g-1));
            out = c*atan(sqrt((M.^2 - 1)/c^2)) - atan(sqrt(M.^2 - 1));
        case 'M_from_nu'
            [nu, g] = varargin{:};
            out = zeros(size(nu));
            for k = 1:numel(nu)
                if nu(k) <= 0, out(k) = 1; continue; end
                out(k) = bisect(@(M) gasdyn('nu', M, g) - nu(k), 1, 100);
            end
        case 'mu'
            M = varargin{1}; out = asin(1./M);
        case 'normal_shock'
            [M1, g] = varargin{:};
            S.M2 = sqrt((1 + (g-1)/2*M1.^2)./(g*M1.^2 - (g-1)/2));
            S.p2_p1 = 1 + 2*g/(g+1)*(M1.^2 - 1);
            S.rho2_rho1 = (g+1)*M1.^2./((g-1)*M1.^2 + 2);
            S.T2_T1 = S.p2_p1./S.rho2_rho1;
            S.p02_p01 = S.p2_p1.*gasdyn('p0_p', S.M2, g)./gasdyn('p0_p', M1, g);
            out = S;
        otherwise
            error('gasdyn:fn', 'Unknown relation %s', fn);
    end
end

function x = bisect(f, lo, hi)
    flo = f(lo);
    for it = 1:200
        x = 0.5*(lo + hi); fx = f(x);
        if fx == 0 || (hi - lo) < 1e-14*max(1, abs(x)), return; end
        if sign(fx) == sign(flo), lo = x; flo = fx; else, hi = x; end
    end
end
