%% plan_nozzle.m — the design, every Flow Simulation setting, and the predictions
%
% Run from the project folder:  run('src/plan_nozzle.m')
% Also writes the contour for the CAD: data/results/contour.csv (metres) and
% cad/nozzle_wall_mm.txt (X Y Z in millimetres, for Curve Through XYZ Points).

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..');
addpath(here);
z = load_nozzle();
g = z.gamma; hs = z.throat_half_height; pa = z.p_ambient; mm = 1000;

fprintf('==== %s: exit Mach %.2f, gamma %.2f, T0 %.0f K\n', z.name, z.exit_mach, g, z.T0);
fprintf('Throat half-height %.1f mm, exit half-height %.2f mm, throat to exit %.1f mm\n', ...
        mm*hs, mm*z.y_exit, mm*z.x_exit);
fprintf('Height ratio %.4f with %d characteristics (exact A/A* %.4f, %+.2f%%); initial wall angle %.2f deg\n\n', ...
        z.ratio, z.n_characteristics, z.A_ratio_exact, 100*(z.ratio/z.A_ratio_exact - 1), z.N.theta_max*180/pi);

% Contour files
fid = fopen(fullfile(root, 'data', 'results', 'contour.csv'), 'w');
fprintf(fid, 'x_m,y_m\n'); fprintf(fid, '%.8f,%.8f\n', [z.x z.y]'); fclose(fid);
fid = fopen(fullfile(root, 'cad', 'nozzle_wall_mm.txt'), 'w');
fprintf(fid, '%.5f\t%.5f\t0\n', mm*[z.x z.y]'); fclose(fid);

% Domain and mesh
x_lid = z.x_inlet - z.lid_thickness;
dom = [x_lid - 0.01, z.x_exit + 1.2, 0, 0.3];
fprintf('Computational domain (2D, XY plane), metres: X %.4f to %.4f, Y %.1f to %.1f (Y min symmetry)\n', dom);
nb = [131 30; 262 60; 524 120];
names = {'coarse', 'medium', 'fine'};
fprintf('Basic mesh: coarse %d x %d, medium %d x %d, fine %d x %d (cells 10, 5, 2.5 mm)\n', nb');
fprintf('Local mesh 1, nozzle: region X %.4f to %.4f, Y 0 to %.4f, level 4 (0.63 / 0.31 / 0.16 mm; %.0f / %.0f / %.0f cells across the throat half-height)\n', ...
        x_lid, z.x_exit + 0.005, z.y_exit + z.wall_thickness, hs/0.000625, hs/0.0003125, hs/0.00015625);
fprintf('Local mesh 2, plume:  region X %.4f to %.4f, Y 0 to 0.06, level 3 (1.25 / 0.63 / 0.31 mm)\n\n', z.x_exit, 0.6);

% Cases
fprintf('Cases (ambient %.0f Pa, chamber total temperature %.0f K):\n', pa, z.T0);
fprintf('  %-7s %10s %8s %8s %-14s %8s %8s\n', 'case', 'pc [MPa]', 'pc/pa', 'pe/pa', 'regime', 'C_F', 'Me');
cs = fieldnames(z.cases);
for k = 1:numel(cs)
    pc = z.cases.(cs{k});
    Q = quasi1d(z.ratio, pa/pc, g);
    Pf = nozzle_performance(z.ratio, pc, pa, z.T0, g, z.R_gas);
    fprintf('  %-7s %10.3f %8.2f %8.3f %-14s %8.4f %8.3f\n', cs{k}, pc/1e6, pc/pa, Pf.pe/pa, Q.regime, Pf.CF, Pf.Me);
end
Q = quasi1d(z.ratio, pa/z.cases.design, g);
fprintf('Normal shock at the exit when pc = %.3f MPa; choking needs pc above %.4f MPa\n', ...
        pa/Q.pb_limits(2)/1e6, pa/Q.pb_limits(1)/1e6);
mdot = nozzle_performance(z.ratio, z.cases.design, pa, z.T0, g, z.R_gas).mdot_per_area*2*hs;
fprintf('Design mass flow: %.3f kg/s per metre of depth (whole nozzle); %.4f kg/s for the half model and a %.0f mm slab\n\n', ...
        mdot, mdot/2*z.depth, mm*z.depth);

% Separation (over-expanded) and shock cells (under-expanded)
pc = z.cases.over;
xs = linspace(0, z.x_exit, 400);
pw = pc*wall_pressure_theory(z, xs);
i = find(pw <= 0.4*pa, 1);
if isempty(i)
    fprintf('Over-expanded: wall pressure stays above 0.4 pa; Summerfield predicts no separation\n');
else
    fprintf('Over-expanded: inviscid wall pressure falls to 0.4 pa at x = %.1f mm (%.0f%% of the divergent length): Summerfield separation\n', ...
            mm*xs(i), 100*xs(i)/z.x_exit);
end
[L, Mj, bj] = shock_cell_spacing(z.cases.under, pa, 2*hs, g);
fprintf('Under-expanded: fully expanded Mach %.3f, jet height %.1f mm, shock-cell length %.0f mm (linear theory)\n', Mj, mm*bj, mm*L);
