function N = moc_nozzle(Me, g, n, hs)
%MOC_NOZZLE  Minimum-length planar nozzle by the method of characteristics.
%
%   N = moc_nozzle(Me, g, n, hs)
%
%   Me   design exit Mach number      g   ratio of specific heats
%   n    number of characteristics from the throat corner
%   hs   throat half-height (the axis is y = 0)
%
%   The flow at the throat is uniform and sonic. A sharp corner at the throat
%   turns the wall by theta_max = nu(Me)/2, spreading an expansion fan of n
%   C- waves. Each wave reflects from the axis as a C+ wave, and the wall is
%   shaped to cancel every reflected wave where it arrives, leaving uniform
%   flow at Me (Anderson, Modern Compressible Flow, section 11.7).
%
%   Riemann invariants: K- = theta + nu along C-, K+ = theta - nu along C+.
%   With fan angles theta_k = theta_max k/n, C- wave k carries K- = 2 theta_k
%   and C+ wave j (starting on the axis, theta = 0) carries K+ = -2 theta_j, so
%   point (j,k) has  theta = theta_k - theta_j,  nu = theta_k + theta_j.
%   Positions come from straight-line segments whose slopes are the averages
%   of theta -/+ mu at their two ends.
%
%   N.xw, N.yw   wall contour from the throat corner to the exit
%   N.Mw         Mach number along the wall (the inviscid wall pressure follows)
%   N.ratio      exit half-height / throat half-height (theory: A/A*(Me))
%   N.theta_max  initial wall angle, rad
%   N.X, N.Y, N.Mnet  the characteristic net (NaN where no point exists)

    nu_e = gasdyn('nu', Me, g);
    tmax = nu_e/2;
    th = tmax*(1:n)/n;                              % fan angles
    Mfan = gasdyn('M_from_nu', th, g);              % each fan wave at the corner: theta = nu = th(k)
    mufan = gasdyn('mu', Mfan);

    [X, Y, TH, NU, MU] = deal(NaN(n, n));
    [xw, yw, tw, Mw] = deal(zeros(n + 1, 1));
    xw(1) = 0; yw(1) = hs; tw(1) = tmax; Mw(1) = gasdyn('M_from_nu', tmax, g);
    for j = 1:n
        for k = j:n
            TH(j,k) = th(k) - th(j); NU(j,k) = th(k) + th(j);
            MU(j,k) = gasdyn('mu', gasdyn('M_from_nu', NU(j,k), g));
            % the C- through this point comes from the corner (j = 1) or from point (j-1,k)
            if j == 1
                x1 = 0; y1 = hs; s1 = tan(0.5*((th(k) - mufan(k)) + (TH(j,k) - MU(j,k))));
            else
                x1 = X(j-1,k); y1 = Y(j-1,k); s1 = tan(0.5*((TH(j-1,k) - MU(j-1,k)) + (TH(j,k) - MU(j,k))));
            end
            if k == j                                % on the axis
                X(j,k) = x1 - y1/s1; Y(j,k) = 0;
            else                                     % meet the C+ from point (j,k-1)
                x2 = X(j,k-1); y2 = Y(j,k-1);
                s2 = tan(0.5*((TH(j,k-1) + MU(j,k-1)) + (TH(j,k) + MU(j,k))));
                X(j,k) = (y2 - y1 + s1*x1 - s2*x2)/(s1 - s2);
                Y(j,k) = y1 + s1*(X(j,k) - x1);
            end
        end
        % wall point: the C+ leaves (j,n) unchanged and meets the wall
        tw(j+1) = TH(j,n); Mw(j+1) = gasdyn('M_from_nu', NU(j,n), g);
        s2 = tan(TH(j,n) + MU(j,n));
        sw = tan(0.5*(tw(j) + tw(j+1)));
        x2 = X(j,n); y2 = Y(j,n);
        xw(j+1) = (y2 - yw(j) + sw*xw(j) - s2*x2)/(sw - s2);
        yw(j+1) = yw(j) + sw*(xw(j+1) - xw(j));
    end
    N.xw = xw; N.yw = yw; N.thw = tw; N.Mw = Mw;
    N.ratio = yw(end)/hs; N.theta_max = tmax;
    N.X = X; N.Y = Y; N.Mnet = arrayfun(@(v) ifnan_mach(v, g), NU);
    N.length = xw(end);
end

function M = ifnan_mach(nu, g)
    if isnan(nu), M = NaN; else, M = gasdyn('M_from_nu', nu, g); end
end
