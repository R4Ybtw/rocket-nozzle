function [L, Mj, bj] = shock_cell_spacing(pc, pa, bs, g)
%SHOCK_CELL_SPACING  Shock-cell length of an under-expanded planar jet (linear theory).
%
%   [L, Mj, bj] = shock_cell_spacing(pc, pa, bs, g)
%
%   pc, pa   chamber stagnation and ambient pressure
%   bs       full throat height (twice the half-height)
%   Mj       fully expanded jet Mach number: isentropic from pc to pa
%   bj       fully expanded jet height: bs * A/A*(Mj)
%   L        shock-cell length  2 bj sqrt(Mj^2 - 1)
%
%   Derivation: small pressure disturbances in a uniform supersonic jet obey
%   (M^2 - 1) phi_xx = phi_yy. With no cross-flow on the axis and ambient
%   pressure (phi_x = 0) on the jet boundary y = bj/2, the lowest mode is
%   cos(pi y/bj) cos(k x) with k = pi/(bj sqrt(M^2 - 1)), whose wavelength
%   2 pi/k is the cell length above. (The same argument for a round jet gives
%   Prandtl's 1.306 D sqrt(M^2 - 1).) Linear theory assumes a weak mismatch;
%   at pe/pa = 2 expect it to be approximate, which the CFD measures.

    Mj = gasdyn('M_from_p', pc/pa, g);
    bj = bs*gasdyn('A_As', Mj, g);
    L = 2*bj*sqrt(Mj^2 - 1);
end
