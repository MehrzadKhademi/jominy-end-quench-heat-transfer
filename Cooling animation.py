# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 21:44:21 2025

@author: Windows-24H2
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# -----------------------------
# Load NEW simulation results
# -----------------------------
# Put these .npy files in the same folder as this script,
# or replace with full paths.
x = np.load("x.npy")                 # [m] shape (nx,)
t = np.load("t.npy")                 # [s] shape (nt,)
T_hist = np.load("T_hist.npy")       # [°C] shape (nt, nx)

nx = x.size
nt = t.size
print(f"Loaded: nx={nx}, nt={nt}")

# -----------------------------
# Animation settings
# -----------------------------
pause_seconds = 5.0   # viewer pause before cooling begins (visual hold only)
fps = 15              # output frames per second
step = 5              # take every 'step' simulation time steps (dt_anim = step * dt)

# Build frame indices: repeat t=0 frame for pause, then show simulation frames
frames_sim = list(range(0, nt, step))
pause_frames = int(round(pause_seconds * fps))
frames = [0] * pause_frames + frames_sim

# Color limits (fixed across frames)
T_min = float(np.min(T_hist))
T_max = float(np.max(T_hist))

# -----------------------------
# Figure setup (1D heatmap)
# -----------------------------
fig, ax = plt.subplots(figsize=(9, 2.2))

im = ax.imshow(
    T_hist[0, :][np.newaxis, :],
    aspect="auto",
    origin="lower",
    extent=[x[0] * 1000.0, x[-1] * 1000.0, 0, 1],  # x in mm
    vmin=T_min,
    vmax=T_max
)

ax.set_xlabel("Position along bar [mm]")
ax.set_yticks([])
ax.set_title("Cooling of 50 mm bar (Option-B convection at x=0, insulated at x=L)")

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Temperature [°C]")

time_text = ax.text(
    0.02, 0.78, "t = 0.00 s (water jet OFF)",
    transform=ax.transAxes,
    color="white"
)

# -----------------------------
# Frame update
# -----------------------------
def update(frame_idx):
    n = frames[frame_idx]
    im.set_data(T_hist[n, :][np.newaxis, :])

    if frame_idx < pause_frames:
        time_text.set_text("t = 0.00 s (water jet OFF)")
    else:
        time_text.set_text(f"t = {t[n]:.2f} s (water jet ON)")

    return im, time_text

interval_ms = 1000.0 / fps

ani = FuncAnimation(
    fig,
    update,
    frames=len(frames),
    interval=interval_ms,
    blit=True
)

plt.tight_layout()

# -----------------------------
# Save output
# -----------------------------
# MP4 requires ffmpeg installed and available on PATH.
ani.save("jominy_cooling_new.mp4", writer="ffmpeg", fps=fps)

# GIF option (slower + bigger file):
# ani.save("jominy_cooling_new.gif", writer="pillow", fps=fps)

plt.show()
