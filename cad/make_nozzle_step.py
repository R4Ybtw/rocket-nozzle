"""Build the half-nozzle body as a STEP file for SolidWorks Flow Simulation.

    pip install cadquery
    python make_nozzle_step.py          (run from the cad/ folder)

Reads ../data/results/contour.csv (written by src/plan_nozzle.m) and
../data/nozzle.json. The body is the upper wall of the planar nozzle plus the
lid that closes the chamber: the lid's inner face (x = inlet) carries the inlet
boundary condition, and the axis y = 0 is the symmetry plane. Millimetres;
x = 0 at the throat; extruded along +Z by the slab depth.
"""
import csv
import json
import cadquery as cq

z = json.load(open("../data/nozzle.json"))
pts = [(1000 * float(r["x_m"]), 1000 * float(r["y_m"])) for r in csv.DictReader(open("../data/results/contour.csv"))]
tw, tl, depth = 1000 * z["wall_thickness"], 1000 * z["lid_thickness"], 1000 * z["depth"]
x_in, x_ex, y_ex = pts[0][0], pts[-1][0], pts[-1][1]
x_lid, y_top = x_in - tl, y_ex + tw

outline = [(x_lid, 0.0), (x_in, 0.0)] + pts + [(x_ex, y_top), (x_lid, y_top)]
body = cq.Workplane("XY").polyline(outline).close().extrude(depth)
cq.exporters.export(body, "nozzle_half_2d.step")
s = body.val()
bb = s.BoundingBox()
print(f"nozzle_half_2d: x {bb.xmin:.2f}..{bb.xmax:.2f} mm, y {bb.ymin:.2f}..{bb.ymax:.2f} mm, depth {bb.zlen:.1f} mm, "
      f"exit half-height {y_ex:.3f} mm, valid={s.isValid()}")
