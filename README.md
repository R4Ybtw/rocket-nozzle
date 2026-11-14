# Planar rocket nozzle: method of characteristics and compressible CFD

Does compressible CFD reproduce the design, over-expanded and under-expanded behaviour of a nozzle designed by the method of characteristics? The same nozzle behaves differently as the ambient pressure changes: matched at its design altitude, shocking and separating when the outside pressure is too high, and forming shock diamonds when it is too low. That is why rocket plumes look different at sea level and at altitude.

| Item | Value |
|---|---|
| Nozzle | Planar (two-dimensional) minimum-length nozzle for Mach 3, designed with 40 characteristics; throat height 10 mm, exit height 42.2 mm, 84 mm from throat to exit |
| Gas | Cold air, γ = 1.4, chamber total temperature 300 K, so no combustion is modelled |
| Cases | Chamber pressure 3.72 MPa (design for 101.3 kPa ambient), 1.00 MPa (over-expanded) and 7.44 MPa (under-expanded) |

A planar nozzle keeps the method of characteristics exact in two dimensions, and the exit-to-throat height ratio must equal the isentropic area ratio, which makes the design itself testable. Planar nozzles are flown on the F-22, whose two-dimensional nozzles vector thrust.

*Status: the design and theory tools are complete and verified (V1–V3). Flow Simulation results will be added as the runs are completed.*

## Verification and validation

`src/analyze_nozzle.m` runs every step on whatever data exists and writes the verdicts to [`data/results/summary.md`](data/results/summary.md).

| Step | Question | Method and acceptance |
|---|---|---|
| V1 Code verification | Are the gas-dynamic tools correct? | Isentropic, normal-shock and Prandtl–Meyer relations against textbook values; quasi-1D shock solutions reproduce their back pressure; thrust coefficient by momentum equals the closed form; 36 checks |
| V2 MOC convergence | Does the characteristic net converge to the exact design? | Exit height ratio against A/A* = 4.2346 as the characteristics double (−0.26% at 40, −0.04% at 160) |
| V3 Theory | What should each case do? | Quasi-1D regime, exit pressure and ideal thrust coefficient; Summerfield separation point; linear-theory shock-cell length |
| V4 Iterative convergence | Is each run converged? | Goal convergence before the iteration limit |
| V5 Mass conservation | Is the flow choked at the right mass flow? | Inlet mass flow, exit-plane mass flow and the choked-flow formula within 1% |
| V6 Mesh convergence | How large is the discretisation error? | Three meshes, every cell halved, design case; grid convergence index on the thrust coefficient |
| V7 Design point | Does CFD reproduce the design? | Exit Mach number, thrust coefficient against 1.4525, wall pressure against the characteristic solution |
| V8 Over-expanded | Where does the flow separate? | Departure of the wall pressure from attached flow against the Summerfield criterion (separation near 0.4 of ambient; 0.2–0.5 in tests) |
| V9 Under-expanded | How long are the shock cells? | Axis-pressure peaks against linear theory, L = 2 b_j √(M_j² − 1), derived in `src/shock_cell_spacing.m` |

## Predictions

| Case | pc/pa | pe/pa | Inviscid regime | Ideal C_F | Prediction |
|---|---|---|---|---|---|
| Design | 36.7 | 1.00 | Matched | 1.4525 | Uniform Mach 3 exit |
| Over-expanded | 9.87 | 0.27 | Oblique shocks at the exit | 1.1396 | Summerfield separation at x = 55 mm (66% of the divergent length) |
| Under-expanded | 73.4 | 2.01 | Expansion fans at the exit | 1.5100 | Shock cells 441 mm long, fully expanded Mach 3.47 |

## Simulation setup

| Item | Setting |
|---|---|
| Solver | SolidWorks Flow Simulation, 2D (XY plane), external, steady, compressible, high-Mach option |
| Model | Upper half of the nozzle with symmetry on the axis; total pressure and temperature on the inner face of the chamber lid; still air at 101.3 kPa outside |
| Domain | From the chamber lid to 1.2 m behind the exit, 0.3 m from the axis |
| Mesh | Basic cells of 10, 5 and 2.5 mm; nozzle region refined 4 levels and plume region 3 levels, so every cell halves between meshes (8, 16 and 32 cells across the throat half-height) |
| Outputs | Mass flow at the inlet; pressure, density and axial velocity across the exit plane; pressure along the wall and the jet axis |

## Repository layout

| Path | Contents |
|---|---|
| `src/analyze_nozzle.m` | Runs V1–V9 and writes `data/results/summary.md` and the figure |
| `src/plan_nozzle.m` | The design, every Flow Simulation setting and the predictions; writes the contour |
| `src/moc_nozzle.m`, `src/load_nozzle.m` | Method-of-characteristics design and the full contour |
| `src/gasdyn.m`, `src/quasi1d.m`, `src/nozzle_performance.m`, `src/shock_cell_spacing.m`, `src/wall_pressure_theory.m` | Theory |
| `src/exit_plane.m`, `src/cell_spacing_from_axis.m`, `src/mesh_study.m` | CFD post-processing |
| `src/convert_exports.m`, `src/xyplot_table.m`, `src/read_profile.m`, `src/read_runs.m` | Readers for the Flow Simulation exports |
| `tests/test_nozzle.m` | Verification suite |
| `data/nozzle.json` | Design inputs and the three cases |
| `data/results/` | Contour, runs table, converted profiles and `summary.md` |
| `cad/` | STEP of the half nozzle with its chamber lid, the wall curve, and the script that builds them |
| `docs/run_log.md` | Run-by-run record of settings |

## Reproducing

```matlab
run('src/plan_nozzle.m')      % design, contour, settings, predictions
run('src/analyze_nozzle.m')   % tests, then V1–V9
```

The STEP file is rebuilt with `python make_nozzle_step.py` in `cad/` (CadQuery). The MATLAB code also runs in GNU Octave, except the conversion of the Excel exports, which needs MATLAB.

## Limitations

- **Cold air:** γ is constant at 1.4; a hot rocket exhaust has lower γ and changing composition.
- **Planar geometry:** a round nozzle has different separation and shock-cell behaviour; the planar formulas used here do not transfer directly.
- **Steady RANS:** nozzle separation is unsteady in reality and can switch between patterns; a steady solution gives one representative state.
- **Sharp throat:** the minimum-length design turns the flow at a sharp corner, which real nozzles round off.

## References

- J. D. Anderson, *Modern Compressible Flow*, McGraw-Hill: method of characteristics and the minimum-length nozzle (chapter 11), isentropic and shock relations.
- G. P. Sutton and O. Biblarz, *Rocket Propulsion Elements*, Wiley: thrust coefficient.
- R. H. Schmucker, *Status of flow separation prediction in liquid propellant rocket nozzles*, NASA TM X-64890 (1974): the Summerfield criterion (separation where wall pressure falls about 60% below ambient) and its measured range.
- I. B. Celik et al., "Procedure for estimation and reporting of uncertainty due to discretization in CFD applications," *Journal of Fluids Engineering* 130(7), 078001 (2008).

## Licence

Code: MIT (see `LICENSE`). Text and figures: CC BY 4.0.
