# ================================================================
# Epitope Pipeline — Spin Video Generator
# ================================================================
# Usage: Open your epitope PML session first, then run this script:
#   @scripts/spin_video.pml
#
# Output: saves frames + final .mp4 next to the PML file you loaded.
#   e.g. smpdl3b_epitope.pml -> smpdl3b_epitope_spin/
#
# Settings — edit these before running:
#   DRAFT:  36 frames, 480x270  (~1 min render, 1.2 sec video)
#   FINAL: 180 frames, 1920x1080 (~75 min render, 6 sec video)
# ================================================================

# --- Ray trace quality ---
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set specular, 0.4
set ambient, 0.35
set ray_opaque_background, on
bg_color black

# --- Render frames ---
python
import os
import subprocess
from pymol import cmd

# === EDIT THESE ===
n_frames = 36          # Quick test — 1.2 sec spin
width = 346            # 28% wider than 270 draft portrait
height = 480           # Portrait
# ==================

# Auto-detect protein name and source directory
protein_obj = None
for name in cmd.get_names("objects"):
    obj_type = cmd.get_type(name)
    if obj_type == "object:molecule":
        protein_obj = name
        break

from datetime import datetime
base_name = protein_obj or "spin"
timestamp = datetime.now().strftime("%y%m%d_%H%M")

# Get the directory where the PML session was loaded from
# PyMOL's get_session returns the session path, but we can use
# the log_file setting or just check where the hidden PDB lives
base_dir = None
if protein_obj:
    # The PML loads ".{name}.pdb" from its own directory
    # Walk up from cwd looking for it
    check = os.getcwd()
    for _ in range(5):
        candidate = os.path.join(check, ".{}.pdb".format(protein_obj.lower()))
        if os.path.exists(candidate):
            base_dir = check
            break
        # Also check Structures/zone_sessions for bispecific
        for sub in ["", "Structures", "Structures/zone_sessions", "zone_sessions"]:
            candidate = os.path.join(check, sub, ".{}.pdb".format(protein_obj.lower())) if sub else candidate
            if os.path.exists(candidate):
                base_dir = os.path.join(check, sub) if sub else check
                break
        if base_dir:
            break
        check = os.path.dirname(check)

if not base_dir:
    base_dir = os.getcwd()
    print("Could not find PDB source directory, saving to: {}".format(base_dir))

output_dir = os.path.join(base_dir, "{}_spin_draft_{}".format(base_name, timestamp))
os.makedirs(output_dir, exist_ok=True)

if protein_obj is None:
    print("ERROR: No molecule object found!")
else:
    degrees_per_frame = 360.0 / n_frames
    print("Spinning: {} ({} frames at {}x{})".format(protein_obj, n_frames, width, height))
    print("Output: {}".format(output_dir))

    for i in range(n_frames):
        cmd.rotate("y", degrees_per_frame, protein_obj)
        cmd.ray(width, height)
        cmd.png(os.path.join(output_dir, "frame_{:04d}.png".format(i)))
        print("Frame {}/{} done".format(i + 1, n_frames))

    # Auto-stitch with ffmpeg
    mp4_name = "{}_spin_draft_{}.mp4".format(base_name, timestamp)
    mp4_full = os.path.join(output_dir, mp4_name)
    # Use unique /tmp symlink to avoid spaces-in-path issues with ffmpeg
    tmp_link = "/tmp/pymol_spin_{}".format(timestamp)
    try:
        if os.path.islink(tmp_link):
            os.unlink(tmp_link)
        os.symlink(output_dir, tmp_link)
        subprocess.run([
            "/opt/homebrew/bin/ffmpeg", "-y", "-framerate", "30",
            "-i", os.path.join(tmp_link, "frame_%04d.png"),
            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
            os.path.join(tmp_link, mp4_name),
        ], check=True, capture_output=True)
        print("\nVideo saved: {}".format(mp4_full))
    except FileNotFoundError:
        print("\nffmpeg not found at /opt/homebrew/bin/ffmpeg")
        print("Stitch manually:")
        print("  ffmpeg -framerate 30 -i /tmp/pymol_spin_{}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 output.mp4".format(timestamp))
    except Exception as e:
        print("\nffmpeg failed: {}".format(e))
python end
