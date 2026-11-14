function p_p0 = wall_pressure_theory(z, x)
%WALL_PRESSURE_THEORY  Inviscid wall pressure along the nozzle, p/p0, at positions x.
%
%   Converging part (x < 0): quasi-1D subsonic flow at the local area ratio.
%   Diverging part: the Mach number the method of characteristics puts on the
%   wall (the wall just behind the corner sees the fully turned flow,
%   nu = theta_max), linearly interpolated between wall points.
%   Beyond the exit the design exit pressure is held.

    g = z.gamma; hs = z.throat_half_height;
    p_p0 = zeros(size(x));
    for k = 1:numel(x)
        if x(k) < 0
            y = interp1(z.x, z.y, x(k));
            M = gasdyn('M_from_A', y/hs, g, 'sub');
        else
            M = interp1(z.N.xw, z.N.Mw, min(x(k), z.x_exit));
        end
        p_p0(k) = 1/gasdyn('p0_p', M, g);
    end
end
