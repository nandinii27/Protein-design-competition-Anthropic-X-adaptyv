"""
Publication figures for one designed pMHC complex.

    pip install pymol-open-source
    python render_design.py md_in/pmhc_g12d_06.pdb figs/

Produces three images on a MATCHED camera, which is the point: the complex,
the target alone, and the complex with contacted peptide residues coloured.
Same viewpoint across all three means the reader compares structures rather
than compensating for a rotation.

Chains: A binder, B HLA alpha1-alpha2, C beta2m, D peptide.
Prints the contacted peptide residues, which is the declared epitope.

The craft here is entirely in the scene: muted greys for context, one
saturated colour for the subject, one accent for the single residue the
whole argument is about, shadows off, depth cue off, antialias up.
"""
import sys
from pathlib import Path

import pymol
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd

src = Path(sys.argv[1] if len(sys.argv) > 1 else "pmhc_g12d_06.pdb")
figs = Path(sys.argv[2] if len(sys.argv) > 2 else "figs")
figs.mkdir(parents=True, exist_ok=True)
stem = src.stem

cmd.load(str(src), "cx")
cmd.hide("everything"); cmd.bg_color("white")
cmd.show("cartoon", "polymer")
cmd.color("grey75", "chain B"); cmd.color("grey88", "chain C")
cmd.color("0xE8A33D", "chain D"); cmd.color("0x2E86AB", "chain A")
cmd.set("cartoon_side_chain_helper", 1)
cmd.show("sticks", "chain D and sidechain")
cmd.color("0xC1272D", "chain D and resi 6")
cmd.show("spheres", "chain D and resi 6 and sidechain")
cmd.set("sphere_scale", 0.32, "chain D and resi 6")
for k, v in [("ray_shadows",0),("antialias",2),("depth_cue",0),("ambient",0.20),
             ("specular",0.2),("ray_opaque_background",1),("cartoon_fancy_helices",1)]:
    cmd.set(k, v)

cmd.orient("chain B+D"); cmd.turn("x", 90); cmd.zoom("polymer", 2)
view = cmd.get_view()

cmd.ray(1800,1350); cmd.png(str(figs / f"{stem}_1_complex.png"), dpi=300)

cmd.hide("everything", "chain A")          # same camera, binder removed
cmd.set_view(view)
cmd.ray(1800,1350); cmd.png(str(figs / f"{stem}_2_target.png"), dpi=300)

# contact residues on the peptide, from the complex
cmd.show("cartoon", "chain A")
cmd.set_view(view)
sele = "chain D within 4.5 of chain A"
cmd.select("epi", sele)
print("peptide residues contacted:",
      sorted({int(a.resi) for a in cmd.get_model("epi").atom}))
cmd.color("0x7B3294", "byres epi")
cmd.color("0xC1272D", "chain D and resi 6")
cmd.ray(1800,1350); cmd.png(str(figs / f"{stem}_3_epitope.png"), dpi=300)
print("done")
