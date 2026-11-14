%% analyze_nozzle.m — every verification and validation step, with verdicts
%
% Run from the project folder:  run('src/analyze_nozzle.m')
% Uses whatever results exist; steps without data are reported as pending.
% Writes data/results/summary.md and figures/nozzle_profiles.png.
%
%   V1  Code verification       gas-dynamic relations, MOC, quasi-1D, integrals (test suite)
%   V2  MOC convergence         exit height ratio against A/A* as characteristics double
%   V3  Theory                  regime, exit pressure, thrust coefficient for the three cases
%   V4  Iterative convergence   which runs reached a steady solution
%   V5  Mass conservation       inlet mass flow = exit-plane mass flow = choked theory
%   V6  Mesh convergence        design case, three meshes: thrust coefficient
%   V7  Design point            exit Mach, thrust coefficient, wall pressure against MOC
%   V8  Over-expanded           separation location against the Summerfield criterion
%   V9  Under-expanded          shock-cell length against linear theory

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..');
addpath(here); addpath(fullfile(root, 'tests'));
res = fullfile(root, 'data', 'results');
z = load_nozzle(); g = z.gamma; hs = z.throat_half_height; pa = z.p_ambient; Rg = z.R_gas;
out = {}; st = {}; mm = 1000;
num = @(v, f) ifelse_str(isfinite(v), sprintf(f, v), '–');

%% V1
[ok, np, nt] = test_nozzle();
if ~ok, error('analyze_nozzle:tests', '%d of %d checks failed.', nt - np, nt); end
out = rpt(out, '\n## V1. Code verification');
out = rpt(out, '\nAll %d checks in `tests/test_nozzle.m` pass: isentropic, normal-shock and Prandtl–Meyer relations against textbook values; quasi-1D shock solutions reproduce their back pressure; the thrust coefficient from the exit momentum equals the closed form; the exit-plane integrals recover the choked mass flow exactly.', nt);
st(end+1, :) = {'V1 Code verification', sprintf('pass: %d checks', nt)};

%% V2
out = rpt(out, '\n## V2. Method of characteristics: convergence');
ex = z.A_ratio_exact; ns = [10 20 40 80 160]; e = zeros(size(ns));
for k = 1:numel(ns), N = moc_nozzle(z.exit_mach, g, ns(k), 1); e(k) = N.ratio/ex - 1; end
ord = log2(abs(e(1:end-1)./e(2:end)));
out = rpt(out, '');
out = rpt(out, '| Characteristics | Exit height ratio error | Observed order |');
out = rpt(out, '|---|---|---|');
out = rpt(out, '| %d | %+.3f%% | |', ns(1), 100*e(1));
for k = 2:numel(ns), out = rpt(out, '| %d | %+.3f%% | %.2f |', ns(k), 100*e(k), ord(k-1)); end
out = rpt(out, '\nThe exit-to-throat height ratio converges to the exact A/A* = %.4f. The built contour uses %d characteristics (%+.2f%%). The observed order, about 1.3, reflects the straight-segment construction around the centred expansion at the sharp throat corner.', ...
          ex, z.n_characteristics, 100*(z.ratio/ex - 1));
st(end+1, :) = {'V2 MOC convergence', sprintf('error %+.2f%% at %d characteristics, order %.2f', 100*(z.ratio/ex - 1), z.n_characteristics, ord(end))};

%% V3
out = rpt(out, '\n## V3. Theory for the three cases');
cs = {'design', 'over', 'under'}; TH = struct();
out = rpt(out, '');
out = rpt(out, '| Case | Chamber pressure (MPa) | pc/pa | pe/pa | Regime (inviscid) | Ideal C_F | Prediction |');
out = rpt(out, '|---|---|---|---|---|---|---|');
for k = 1:3
    pc = z.cases.(cs{k});
    Q = quasi1d(z.ratio, pa/pc, g); Pf = nozzle_performance(z.ratio, pc, pa, z.T0, g, Rg);
    TH.(cs{k}) = struct('pc', pc, 'Q', Q, 'P', Pf);
    note = '–';
    if strcmp(cs{k}, 'over')
        xs = linspace(0, z.x_exit, 2000); pw = pc*wall_pressure_theory(z, xs);
        i = find(pw <= 0.4*pa, 1); TH.over.x_sep = NaN;
        if ~isempty(i), TH.over.x_sep = xs(i); note = sprintf('Summerfield separation at x = %.1f mm', mm*xs(i)); end
    elseif strcmp(cs{k}, 'under')
        [L, Mj, bj] = shock_cell_spacing(pc, pa, 2*hs, g); TH.under.L = L;
        note = sprintf('shock cells %.0f mm long (M_j %.2f)', mm*L, Mj);
    end
    out = rpt(out, '| %s | %.3f | %.2f | %.3f | %s | %.4f | %s |', cs{k}, pc/1e6, pc/pa, Pf.pe/pa, Q.regime, Pf.CF, note);
end
st(end+1, :) = {'V3 Theory', sprintf('ideal C_F %.4f / %.4f / %.4f (design / over / under)', TH.design.P.CF, TH.over.P.CF, TH.under.P.CF)};

%% V4, V5
R = read_runs(fullfile(res, 'runs.csv'));
out = rpt(out, '\n## V4. Iterative convergence');
if R.n == 0
    out = rpt(out, '\nPending: no runs recorded in data/results/runs.csv.');
    st(end+1, :) = {'V4 Iterative convergence', 'pending'};
else
    out = rpt(out, '\n%d of %d runs converged.', nnz(R.converged), R.n);
    st(end+1, :) = {'V4 Iterative convergence', sprintf('%d of %d converged', nnz(R.converged), R.n)};
end
ncv = convert_exports(root);
if ncv > 0, out = rpt(out, '\n%d XY-plot export(s) converted.', ncv); end
out = rpt(out, '\n## V5. Mass conservation');
rows = {};
for k = 1:R.n
    P = read_profile(fullfile(res, sprintf('%s_%s_exit.csv', R.case{k}, R.mesh{k})));
    if P.n == 0 || ~isfield(TH, R.case{k}), continue; end
    pc = TH.(R.case{k}).pc;
    E = exit_plane(P, pa, pc, hs, g);
    m_in = 2*R.mdot(k)/R.span(k);                        % whole nozzle, per metre of depth
    m_th = TH.(R.case{k}).P.mdot_per_area*2*hs;
    rows(end+1, :) = {R.case{k}, R.mesh{k}, m_in, E.mdot, m_th, E}; %#ok<SAGROW>
end
if isempty(rows)
    out = rpt(out, '\nPending: needs runs with an exit-plane profile.');
    st(end+1, :) = {'V5 Mass conservation', 'pending'};
else
    out = rpt(out, '');
    out = rpt(out, '| Case | Mesh | Inlet (kg/s per m) | Exit plane | Choked theory | Exit vs inlet | Inlet vs theory |');
    out = rpt(out, '|---|---|---|---|---|---|---|');
    worst = 0;
    for k = 1:size(rows, 1)
        a = rows{k,4}/rows{k,3} - 1; b = rows{k,3}/rows{k,5} - 1; worst = max([worst abs(a) abs(b)]);
        out = rpt(out, '| %s | %s | %.3f | %.3f | %.3f | %+.2f%% | %+.2f%% |', rows{k,1:5}, 100*a, 100*b);
    end
    out = rpt(out, '\nVerdict: **%s** (largest difference %.2f%%, tolerance 1%%).', ifelse_str(worst <= 0.01, 'pass', 'fail'), 100*worst);
    st(end+1, :) = {'V5 Mass conservation', sprintf('%s: largest difference %.2f%%', ifelse_str(worst <= 0.01, 'pass', 'fail'), 100*worst)};
end

%% V6
out = rpt(out, '\n## V6. Mesh convergence (design case)');
names = {'coarse', 'medium', 'fine'}; cf = NaN(1, 3);
for m = 1:3
    i = find(strcmp(R.case, 'design') & strcmp(R.mesh, names{m}) & R.converged, 1);
    P = read_profile(fullfile(res, sprintf('design_%s_exit.csv', names{m})));
    if ~isempty(i) && P.n > 0, E = exit_plane(P, pa, TH.design.pc, hs, g); cf(m) = E.CF; end
end
rec = ''; Ucf = NaN;
if all(isfinite(cf))
    G = mesh_study(cf(1), cf(2), cf(3)); rel = G.U./cf;
    k = find(rel <= 0.01, 1); if isempty(k), k = 3; end
    rec = names{k}; Ucf = G.U(k);
    out = rpt(out, '');
    out = rpt(out, '| | Coarse | Medium | Fine | Behaviour | Observed order | Extrapolated |');
    out = rpt(out, '|---|---|---|---|---|---|---|');
    out = rpt(out, '| C_F | %.4f | %.4f | %.4f | %s | %s | %s |', cf, G.type, num(G.p, '%.2f'), num(G.f_ext, '%.4f'));
    out = rpt(out, '| Uncertainty | %.2f%% | %.2f%% | %.2f%% | | | |', 100*rel);
    out = rpt(out, '\nVerdict: the **%s** mesh (C_F uncertainty %.2f%%) is used for the three cases.', rec, 100*rel(k));
    st(end+1, :) = {'V6 Mesh convergence', sprintf('%s; %s mesh, C_F uncertainty %.2f%%', G.type, rec, 100*rel(k))};
else
    out = rpt(out, '\nPending: needs converged design runs with exit profiles on the coarse, medium and fine meshes.');
    st(end+1, :) = {'V6 Mesh convergence', 'pending'};
end
pick = @(c) pick_mesh(res, c, rec);

%% V7
out = rpt(out, '\n## V7. Design point');
[mD, P, W] = pick('design');
if P.n > 0
    E = exit_plane(P, pa, TH.design.pc, hs, g);
    Mdes = z.N.Mw(end);
    out = rpt(out, '\nMesh: %s. Mass-weighted exit Mach %.3f (design %.3f, %+.1f%%); C_F %.4f against ideal %.4f (%+.2f%%)%s.', ...
              mD, E.M_mean, Mdes, 100*(E.M_mean/Mdes - 1), E.CF, TH.design.P.CF, 100*(E.CF/TH.design.P.CF - 1), ...
              ifelse_str(isfinite(Ucf), sprintf(', mesh uncertainty %.4f', Ucf), ''));
    s7 = sprintf('exit Mach %.3f, C_F %+.2f%% vs ideal', E.M_mean, 100*(E.CF/TH.design.P.CF - 1));
    if W.n > 0
        d = W.p/TH.design.pc - wall_pressure_theory(z, W.s);
        in = W.s > 0.05*z.x_exit & W.s < 0.98*z.x_exit;
        out = rpt(out, 'Wall pressure against the method of characteristics: RMS difference %.4f p0 over the divergent wall (exit pressure %.4f p0).', ...
                  sqrt(mean(d(in).^2)), 1/gasdyn('p0_p', Mdes, g));
    end
    st(end+1, :) = {'V7 Design point', s7};
else
    out = rpt(out, '\nPending: needs the design-case exit profile.');
    st(end+1, :) = {'V7 Design point', 'pending'};
end

%% V8
out = rpt(out, '\n## V8. Over-expanded: separation');
[mO, P, W] = pick('over');
if W.n > 0
    pin = TH.over.pc*wall_pressure_theory(z, W.s);
    div = W.s > 0.05*z.x_exit & W.s <= z.x_exit;
    i = find(div & W.p > 1.25*pin, 1);
    if isempty(i)
        out = rpt(out, '\nMesh: %s. The wall pressure follows the attached-flow value to the exit: no separation inside the nozzle.', mO);
        s8 = 'no separation inside the nozzle';
    else
        xs = W.s(i); ratio = pin(i)/pa;
        out = rpt(out, '\nMesh: %s. The wall pressure leaves the attached-flow value at x = %.1f mm (%.0f%% of the divergent length), where the attached-flow wall pressure is %.2f of ambient.', ...
                  mO, mm*xs, 100*xs/z.x_exit, ratio);
        out = rpt(out, 'Summerfield criterion: separation where the wall pressure falls to about 0.4 of ambient (measured nozzles: 0.2–0.5), predicted at x = %.1f mm.', mm*TH.over.x_sep);
        s8 = sprintf('separation at %.1f mm (p/pa %.2f); Summerfield %.1f mm', mm*xs, ratio, mm*TH.over.x_sep);
    end
    if P.n > 0
        E = exit_plane(P, pa, TH.over.pc, hs, g);
        out = rpt(out, 'C_F %.4f against %.4f for attached inviscid flow: separation %s the pressure loss of over-expansion.', E.CF, TH.over.P.CF, ...
                  ifelse_str(E.CF > TH.over.P.CF, 'recovers part of', 'does not recover'));
    end
    st(end+1, :) = {'V8 Over-expanded', s8};
else
    out = rpt(out, '\nPending: needs the over-expanded wall profile.');
    st(end+1, :) = {'V8 Over-expanded', 'pending'};
end

%% V9
out = rpt(out, '\n## V9. Under-expanded: shock cells');
[mU, P, ~, A] = pick('under');
if A.n > 0
    [Lc, xpk] = cell_spacing_from_axis(A.s, A.p, pa, z.x_exit, TH.under.L);
    out = rpt(out, '\nMesh: %s. Axis-pressure peaks at x = %s mm: mean shock-cell length %s mm against %.0f mm from linear theory (%s).', mU, ...
              strjoin(arrayfun(@(v) sprintf('%.0f', mm*v), xpk, 'UniformOutput', false), ', '), num(mm*Lc, '%.0f'), mm*TH.under.L, ...
              ifelse_str(isfinite(Lc), sprintf('%+.0f%%', 100*(Lc/TH.under.L - 1)), 'fewer than two peaks'));
    st(end+1, :) = {'V9 Under-expanded', sprintf('cell length %s mm vs %.0f mm theory', num(mm*Lc, '%.0f'), mm*TH.under.L)};
    if P.n > 0
        E = exit_plane(P, pa, TH.under.pc, hs, g);
        out = rpt(out, 'C_F %.4f against ideal %.4f.', E.CF, TH.under.P.CF);
    end
else
    out = rpt(out, '\nPending: needs the under-expanded axis profile.');
    st(end+1, :) = {'V9 Under-expanded', 'pending'};
end

%% Figure: wall pressure, all cases with data
fig = figure; hold on; grid on;
xs = linspace(z.x_inlet, z.x_exit, 400);
plot(mm*xs, wall_pressure_theory(z, xs), 'k-', 'DisplayName', 'Inviscid (quasi-1D + MOC)');
any_w = false;
for k = 1:3
    [~, ~, W] = pick(cs{k});
    if W.n > 0, plot(mm*W.s, W.p/TH.(cs{k}).pc, '-', 'DisplayName', ['Flow Sim, ' cs{k}]); any_w = true; end
end
plot(mm*[0 z.x_exit], 0.4*pa/TH.over.pc*[1 1], 'k:', 'DisplayName', '0.4 p_a, over-expanded case');
set(gca, 'YScale', 'log'); xlabel('x [mm] (throat at 0)'); ylabel('Wall pressure / chamber pressure'); legend('Location', 'southwest');
if any_w, print(fig, fullfile(root, 'figures', 'nozzle_profiles.png'), '-dpng', '-r110'); end

%% Summary
head = {'# Verification and validation summary', '', ...
        sprintf('Generated by `src/analyze_nozzle.m` on %s. %s; cold air, T0 = %.0f K, ambient %.0f Pa.', datestr(now, 'yyyy-mm-dd'), ...
                z.name, z.T0, pa), '', '| Step | Result |', '|---|---|'};
for k = 1:size(st, 1), head{end+1} = sprintf('| %s | %s |', st{k,1}, st{k,2}); end %#ok<SAGROW>
out = [head, out];
fid = fopen(fullfile(res, 'summary.md'), 'w', 'n', 'UTF-8'); fprintf(fid, '%s\n', out{:}); fclose(fid);
fprintf('\n'); for k = 1:size(st, 1), fprintf('%-26s %s\n', st{k,1}, st{k,2}); end
fprintf('\nSummary written to data/results/summary.md\n');
