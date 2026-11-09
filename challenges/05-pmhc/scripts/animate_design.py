"""
Cinematic animation of a designed pMHC complex.

    pip install prody pymol-open-source      # plus ffmpeg on PATH
    python animate_design.py md_in/pmhc_g12d_06.pdb out.mp4

Three stages, which is the same shape as any molecular animation:

  MOTION   Normal mode analysis. An elastic network model of the complex
           gives its genuine low-frequency collective modes; the lowest one
           is the slowest, largest-amplitude breathing the structure has.
           This is real physics, not a decorative wobble, and unlike MD it
           costs seconds. When a trajectory exists, swap it in here.

  CAMERA   Two views captured with the SAME rotation matrix - wide, then
           tight on the epitope - so elements 9-17 of the view matrix can be
           interpolated exactly. Easing is smootherstep, never linear: a
           linear zoom reads as mechanical. A slow yaw runs underneath.

  LABELS   Burned in with ffmpeg rather than PyMOL, for real typography and
           backing plates. Colons must be escaped or drawtext silently
           truncates the string at the first one.

Chains: A binder, B HLA alpha1-alpha2, C beta2m, D peptide.
"""
from __future__ import annotations

import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

FRAMES = 120
FPS = 20
W, H = 1200, 900
MODE_RMSD = 2.0          # angstroms of deformation at full amplitude
N_STATES = 41            # odd, so the midpoint is the undeformed structure
CYCLE = 40               # frames per breath
YAW = 22.0               # degrees of rotation across the whole shot

FONT_BOLD = "/usr/share/fonts/opentype/inter/InterDisplay-Bold.otf"
FONT_MED = "/usr/share/fonts/opentype/inter/InterDisplay-Medium.otf"


def normal_modes(src: Path, dst: Path) -> None:
    from prody import ANM, confProDy, extendModel, parsePDB, traverseMode, writePDB

    confProDy(verbosity="none")
    protein = parsePDB(str(src)).select("protein")
    anm = ANM("complex")
    anm.buildHessian(protein.select("calpha"), cutoff=15.0)
    anm.calcModes(n_modes=6)

    extended, selection = extendModel(anm, protein.select("calpha"), protein, norm=True)
    ensemble = traverseMode(extended[0], selection,
                            n_steps=(N_STATES - 1) // 2, rmsd=MODE_RMSD)

    out = selection.copy()
    out.delCoordset(0)
    for coords in ensemble.getCoordsets():
        out.addCoordset(coords)
    writePDB(str(dst), out)


def render(states: Path, frames: Path) -> None:
    import pymol

    pymol.finish_launching(["pymol", "-qc"])
    from pymol import cmd

    cmd.load(str(states), "m")
    cmd.dss("m", state=(N_STATES + 1) // 2)   # SS from the UNDEFORMED state, or
    cmd.hide("everything")                     # every helix renders as a tube
    cmd.bg_color("white")
    cmd.show("cartoon", "polymer")
    cmd.color("grey75", "chain B")
    cmd.color("grey88", "chain C")
    cmd.color("0xE8A33D", "chain D")
    cmd.color("0x2E86AB", "chain A")
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.show("sticks", "chain D and sidechain")
    cmd.color("0xC1272D", "chain D and resi 6")
    cmd.show("spheres", "chain D and resi 6 and sidechain")
    cmd.set("sphere_scale", 0.32, "chain D and resi 6")
    for key, value in [("ray_shadows", 0), ("antialias", 1), ("depth_cue", 0),
                       ("ambient", 0.20), ("specular", 0.2),
                       ("ray_opaque_background", 1), ("cartoon_fancy_helices", 1),
                       ("all_states", 0)]:
        cmd.set(key, value)

    cmd.orient("chain B+D")
    cmd.turn("x", 90)
    cmd.zoom("polymer", 2)
    wide = list(cmd.get_view())
    cmd.zoom("chain D and resi 4-9", 13)
    close = list(cmd.get_view())
    assert max(abs(wide[i] - close[i]) for i in range(9)) < 1e-6, "rotation drifted"

    def smooth(t: float) -> float:
        t = min(1.0, max(0.0, t))
        return t * t * t * (t * (t * 6 - 15) + 10)

    for f in range(FRAMES):
        blend = smooth((f - FRAMES * 0.20) / (FRAMES * 0.50))
        cmd.set_view(wide[:9] + [wide[i] + (close[i] - wide[i]) * blend
                                 for i in range(9, 18)])
        cmd.turn("y", YAW * f / FRAMES)
        state = (N_STATES + 1) // 2 + round((N_STATES - 1) / 2
                                            * math.sin(2 * math.pi * f / CYCLE))
        cmd.frame(max(1, min(N_STATES, state)))
        cmd.ray(W, H)
        cmd.png(str(frames / f"f{f:04d}.png"))


def escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace(":", "\\:").replace("'", "\\'")


def encode(frames: Path, out: Path, cards, footer: str) -> None:
    def drawtext(font, text, colour, size, y, t0, t1):
        return (f"drawtext=fontfile={font}:text='{escape(text)}':fontcolor=0x{colour}"
                f":fontsize={size}:x=56:y={y}:box=1:boxcolor=white@0.82"
                f":boxborderw=14:enable='between(t,{t0},{t1})'")

    parts = ["[0:v]tpad=stop_mode=clone:stop_duration=2"]
    for head, sub, t0, t1, colour in cards:
        parts.append(drawtext(FONT_BOLD, head, colour, 44, 48, t0, t1))
        parts.append(drawtext(FONT_MED, sub, "5A5A5A", 28, 112, t0, t1))
    parts.append(drawtext(FONT_MED, footer, "8A8A8A", 22, 822, 0, 99))

    script = frames / "filter.txt"
    script.write_text(",\n".join(parts) + "\n[v]")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS),
                    "-i", str(frames / "f%04d.png"),
                    "-filter_complex_script", str(script), "-map", "[v]",
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                    "-movflags", "+faststart", str(out)], check=True)


def main() -> None:
    src = Path(sys.argv[1])
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "animation.mp4")
    work = Path(tempfile.mkdtemp(prefix="anim_"))
    try:
        states = work / "modes.pdb"
        normal_modes(src, states)
        render(states, work)
        encode(work, out, [
            ("Designing a binder to a buried mutation",
             "HLA-A*11:01   ·   KRAS G12D neoantigen   ·   VVVGADGVGK",
             0, 2.3, "1A1A1A"),
            ("Epitope: peptide residues G4 to G9",
             "the solvent-exposed central bulge", 2.5, 4.8, "1A1A1A"),
            ("Asp6: the G12D mutation",
             "8.8 A² exposed. One carboxylate oxygen fully buried.",
             5.0, 8.0, "C1272D"),
        ], "motion: lowest-frequency normal mode, elastic network model")
        print(f"wrote {out}")
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
