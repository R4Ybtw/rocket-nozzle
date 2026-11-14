function [ok, n_pass, n_total] = test_nozzle()
%TEST_NOZZLE  Verification checks against exact values and independent calculations.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..');
addpath(fullfile(root, 'src'));
r = []; g = 1.4;

%% 1. Gas dynamics against textbook values (M = 2, gamma = 1.4)
r(end+1) = chk('p0/p at M 2', gasdyn('p0_p', 2, g), 7.8244, 1e-4);
r(end+1) = chk('A/A* at M 2', gasdyn('A_As', 2, g), 1.6875, 1e-4);
r(end+1) = chk('Prandtl-Meyer angle at M 2 (deg)', gasdyn('nu', 2, g)*180/pi, 26.380, 1e-4);
S = gasdyn('normal_shock', 2, g);
r(end+1) = chk('Normal shock M 2: downstream Mach', S.M2, 0.57735, 1e-4);
r(end+1) = chk('Normal shock M 2: total-pressure ratio', S.p02_p01, 0.72087, 1e-4);
r(end+1) = chk('Inverse A/A*, supersonic branch', gasdyn('M_from_A', gasdyn('A_As', 3, g), g, 'sup'), 3, 1e-10);
r(end+1) = chk('Inverse A/A*, subsonic branch', gasdyn('M_from_A', gasdyn('A_As', 0.3, g), g, 'sub'), 0.3, 1e-10);
r(end+1) = chk('Inverse Prandtl-Meyer', gasdyn('M_from_nu', gasdyn('nu', 2.5, g), g), 2.5, 1e-10);

%% 2. Method of characteristics: the exit height ratio must converge to A/A*
ex = gasdyn('A_As', 3, g); ns = [10 20 40 80]; e = zeros(size(ns));
for k = 1:4, N = moc_nozzle(3, g, ns(k), 1); e(k) = abs(N.ratio/ex - 1); end
ord = log2(e(1:3)./e(2:4));
fprintf('   MOC exit-ratio error %.2e %.2e %.2e %.2e, observed order %.2f %.2f %.2f\n', e, ord);
r(end+1) = chk('MOC error falls every time n doubles', all(diff(e) < 0), 1, 0);
r(end+1) = chk('MOC exit ratio within 0.3% at 40 characteristics', e(3) < 3e-3, 1, 0);
r(end+1) = chk('MOC observed order between 1 and 2', ord(end) > 1 && ord(end) < 2, 1, 0);
N = moc_nozzle(3, g, 40, 1);
r(end+1) = chk('MOC initial wall angle = nu(Me)/2', N.theta_max, gasdyn('nu', 3, g)/2, 1e-12);
r(end+1) = chk('MOC exit wall Mach = design Mach', N.Mw(end), 3, 1e-9);
r(end+1) = chk('MOC wall rises monotonically', all(diff(N.yw) > 0), 1, 0);

%% 3. Quasi-1D regimes
Q = quasi1d(ex, 0.6, g);
r(end+1) = chk('Shock inside the nozzle at pb/p0 = 0.6', strcmp(Q.regime, 'shock inside'), 1, 0);
r(end+1) = chk('Shock solution reproduces the back pressure', Q.pe_p0, 0.6, 1e-8);
r(end+1) = chk('Exit subsonic behind the shock', Q.Me < 1, 1, 0);
Q = quasi1d(ex, 1/gasdyn('p0_p', 3, g), g);
r(end+1) = chk('Design regime at pe = pb', strcmp(Q.regime, 'design'), 1, 0);
Q2 = quasi1d(ex, Q.pb_limits(2)*1.0001, g);
r(end+1) = chk('Shock at the exit at the limiting back pressure', Q2.A_shock, ex, 1e-3);

%% 4. Thrust coefficient: momentum balance against the closed form
P = nozzle_performance(ex, 3.722e6, 101325, 300, g, 287.05);
r(end+1) = chk('C_F: momentum = closed form', P.CF, P.CF_formula, 1e-10);
r(end+1) = chk('C_F at the design point (closed form 1.4525)', P.CF, 1.4525, 1e-3);
Pv = nozzle_performance(ex, 3.722e6, 0, 300, g, 287.05);
r(end+1) = chk('Vacuum C_F exceeds sea-level C_F', Pv.CF > P.CF, 1, 0);

%% 5. Exit-plane integration on a uniform design-point profile
hs = 0.005; ye = ex*hs; pc = 3.722e6; T0 = 300; Rg = 287.05;
pe = pc/gasdyn('p0_p', 3, g); Te = T0/gasdyn('T0_T', 3, g); rho = pe/(Rg*Te); u = 3*sqrt(g*Rg*Te);
y = linspace(0, ye, 50)';
E = exit_plane(struct('s', y, 'p', pe + 0*y, 'rho', rho + 0*y, 'u', u + 0*y), 101325, pc, hs, g);
r(end+1) = chk('Exit plane: mass flow = choked mass flow', E.mdot, P.mdot_per_area*2*hs, 1e-10);
r(end+1) = chk('Exit plane: C_F = ideal C_F', E.CF, P.CF, 1e-10);
r(end+1) = chk('Exit plane: mean Mach', E.M_mean, 3, 1e-10);

%% 6. Shock cells
[L, Mj] = shock_cell_spacing(7.44e6, 101325, 0.01, g);
r(end+1) = chk('Fully expanded Mach matches pc/pa', gasdyn('p0_p', Mj, g), 7.44e6/101325, 1e-10);
x = linspace(0, 1.2, 2000)'; p = 101325*(1 + 0.5*cos(2*pi*(x - 0.1)/0.4));
Lm = cell_spacing_from_axis(x, p, 101325, 0.05, L);
r(end+1) = chk('Cell spacing from a periodic axis pressure', Lm, 0.4, 2e-3);

%% 7. Wall pressure along the built nozzle
z = load_nozzle();
r(end+1) = chk('Wall pressure at the exit = exit pressure', wall_pressure_theory(z, z.x_exit), ...
               1/gasdyn('p0_p', z.N.Mw(end), g), 1e-10);
r(end+1) = chk('Chamber inlet wall Mach is subsonic', wall_pressure_theory(z, z.x_inlet) > 0.528, 1, 0);

%% 8. Readers
C1 = {'Pressure along exit', [], [], []; 'Model Y [mm]', 'Pressure [Pa]', 'Density (Fluid) [kg/m^3]', 'Velocity (X) [m/s]'; 0, 1e5, 1.2, 600; 20, 1.1e5, 1.3, 610};
T = xyplot_table({C1});
r(end+1) = chk('XY table, one sheet: mm -> m', T.s(2), 0.02, 1e-12);
r(end+1) = chk('XY table, one sheet: velocity column', T.u(2), 610, 1e-12);
Cp = {'Model Y [m]', 'Pressure [Pa]'; 0, 1e5; 0.02, 1.1e5};
Cu = {'Model Y [m]', 'Velocity (X) [m/s]'; 0, 600; 0.01, 605; 0.02, 610};
T = xyplot_table({Cp, Cu});
r(end+1) = chk('XY table, one sheet per parameter: merged', T.u(2), 610, 1e-12);
Ct = {'Model X [m]', 'Total Pressure [Pa]', 'Pressure [Pa]'; 0, 9, 1; 1, 9, 2};
T = xyplot_table({Ct});
r(end+1) = chk('XY table: static pressure, not total', T.p(2), 2, 1e-12);
f = [tempname '.csv'];
fid = fopen(f, 'w'); fprintf(fid, 'case,mesh,mdot_inlet,span_m,converged\ndesign,coarse,0.43,0.01,\nover,coarse,,,\nunder,coarse,-0.8,,0\n'); fclose(fid);
R = read_runs(f); delete(f);
r(end+1) = chk('Runs reader: blank rows skipped, span inherited', R.n + R.span(2), 2.01, 1e-12);
r(end+1) = chk('Runs reader: sign of mass flow ignored', R.mdot(2), 0.8, 1e-12);
G = mesh_study(2.6, 1.4, 1.1);
r(end+1) = chk('Mesh study: Richardson extrapolation', G.f_ext, 1, 1e-12);

ok = all(r); n_pass = nnz(r); n_total = numel(r);
if ok, fprintf('\nALL %d CHECKS PASS\n', n_total); else, fprintf('\n%d OF %d CHECKS FAILED\n', n_total - n_pass, n_total); end
end

function ok = chk(name, got, want, reltol)
    if want == 0, e = abs(got); else, e = abs(got - want)/abs(want); end
    ok = e <= reltol;
    s = 'FAIL'; if ok, s = 'PASS'; end
    fprintf('%-50s got %10.5g  want %10.5g  -> %s\n', name, got, want, s);
end
