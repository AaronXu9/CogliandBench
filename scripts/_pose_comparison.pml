# PyMOL pose comparison render template.
#
# Invoked as:
#   pymol -cq _pose_comparison.pml -- <receptor.pdb> <crystal.sdf> <predicted_path> <mode> <output.png>
#
# mode = "sdf"        -> predicted is an SDF in the same frame as the crystal
#                        (e.g., GNINA); the 4th positional arg is ignored.
# mode = "cif:<chain>" -> predicted is a full AF3 CIF (protein + ligand) in AF3's
#                        own frame. <chain> is the auth_asym_id of the ligand
#                        chain to display (the AF3 system may have multiple
#                        HETATM chains; we filter to the one matching the
#                        crystal ligand by atom count).
#
# Output: receptor cartoon (grey80) + crystal ligand (green) + predicted ligand (magenta).

import sys
from pymol import cmd

args = sys.argv[1:]
sep = args.index("--") if "--" in args else -1
assert sep >= 0, "Expected -- separator before script args"
receptor, crystal, predicted, mode, output = args[sep + 1:sep + 6]

cmd.bg_color("white")

cmd.load(receptor, "receptor")
cmd.load(crystal, "crystal_lig")

if mode == "sdf":
    # Predicted is already in the crystal frame
    cmd.load(predicted, "pred_lig")
elif mode.startswith("cif:"):
    chain_id = mode.split(":", 1)[1]
    # Load the full AF3 prediction (protein + HETATM ligands)
    cmd.load(predicted, "pred_full")
    # Align AF3's protein to the crystal receptor; this carries everything along
    cmd.align("pred_full and polymer", "receptor and polymer")
    # Extract only the HETATMs on the specified chain (the matching ligand)
    sel = f"pred_full and hetatm and chain {chain_id}"
    cmd.create("pred_lig", sel)
    cmd.delete("pred_full")
else:
    raise ValueError(f"Unknown mode: {mode!r}")

cmd.hide("everything")
cmd.show("cartoon", "receptor")
cmd.color("grey80", "receptor")

cmd.show("sticks", "crystal_lig")
cmd.color("green", "crystal_lig")
cmd.util.cnc("crystal_lig")

cmd.show("sticks", "pred_lig")
cmd.color("magenta", "pred_lig")
cmd.util.cnc("pred_lig")

# Frame on the union of crystal and predicted ligand so both are visible
cmd.orient("crystal_lig or pred_lig")
cmd.zoom("crystal_lig or pred_lig", 6)

cmd.set("ray_opaque_background", 0)
cmd.set("ray_shadows", 0)
cmd.set("specular", 0.2)
cmd.set("cartoon_transparency", 0.3)

cmd.ray(1200, 900)
cmd.png(output, dpi=150)
