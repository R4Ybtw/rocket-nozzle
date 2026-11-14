function P = nozzle_performance(ratio_e, pc, pa, T0, g, Rg)
%NOZZLE_PERFORMANCE  Ideal thrust coefficient and mass flow of the planar nozzle.
%
%   P = nozzle_performance(ratio_e, pc, pa, T0, g, Rg)
%
%   Isentropic, shock-free flow to a supersonic exit (valid from just above the
%   shock-at-exit back pressure down to vacuum). Per unit depth, with throat
%   height h*:
%     mdot = p0 A* sqrt(g/(R T0)) (2/(g+1))^((g+1)/(2(g-1)))
%     F    = mdot Ve + (pe - pa) Ae
%     C_F  = F / (p0 A*)
%   P.CF and P.CF_formula compute C_F two ways (momentum from the exit state, and
%   the closed form); the tests require them to agree.
%   P.mdot_per_area = mdot/A* (kg/s per m^2 of throat), P.Me, P.pe, P.Ve.

    Me = gasdyn('M_from_A', ratio_e, g, 'sup');
    pe = pc/gasdyn('p0_p', Me, g);
    Te = T0/gasdyn('T0_T', Me, g);
    Ve = Me*sqrt(g*Rg*Te);
    mA = pc*sqrt(g/(Rg*T0))*(2/(g+1))^((g+1)/(2*(g-1)));
    P.Me = Me; P.pe = pe; P.Ve = Ve; P.mdot_per_area = mA;
    P.CF = (mA*Ve + (pe - pa)*ratio_e)/pc;
    P.CF_formula = sqrt(2*g^2/(g-1)*(2/(g+1))^((g+1)/(g-1))*(1 - (pe/pc)^((g-1)/g))) + (pe - pa)/pc*ratio_e;
end
