function z = load_nozzle(file)
%LOAD_NOZZLE  Read data/nozzle.json and build the nozzle contour.
%
%   z = load_nozzle()
%
%   z.x, z.y      inner wall of the upper half, from the chamber inlet to the
%                 exit (x = 0 at the throat, y = 0 on the axis), metres
%   z.N           the method-of-characteristics design (moc_nozzle)
%   z.ratio       exit/throat height ratio of the built contour
%   z.x_exit, z.y_exit, z.x_inlet
%   Converging part: a cosine from the chamber half-height down to the
%   throat, with zero slope at both ends.

    if nargin < 1
        file = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'nozzle.json');
    end
    z = jsondecode(fileread(file));
    hs = z.throat_half_height; Hc = z.chamber_half_height; Lc = z.converging_length;
    z.N = moc_nozzle(z.exit_mach, z.gamma, z.n_characteristics, hs);
    xc = linspace(-Lc, 0, 41)';
    yc = hs + (Hc - hs)*(1 + cos(pi*(xc + Lc)/Lc))/2;
    z.x = [xc; z.N.xw(2:end)];
    z.y = [yc; z.N.yw(2:end)];
    z.ratio = z.N.ratio;
    z.x_exit = z.N.xw(end); z.y_exit = z.N.yw(end); z.x_inlet = -Lc;
    z.A_ratio_exact = gasdyn('A_As', z.exit_mach, z.gamma);
end
