"""
CMG Momentum Envelope Visualization
------------------------------------
Python conversion of the MATLAB script that plots the angular-momentum
envelope of four CMG (Control Moment Gyroscope) cluster configurations:
Rooftop, Box, Pyramid, and Scissor Pair (6 CMGs).

For each configuration two views are produced:
  1. External Envelope   -> convex hull of all achievable momentum states
  2. Internal Singularities -> every individual gimbal-sign-combination
                                surface, shown semi-transparently so the
                                internal singular surfaces are visible.

Resolution / clarity improvements over the original MATLAB version:
  - Denser spherical sampling grid (200 x 100 instead of 80 x 40)
  - Custom light-source shading (matplotlib LightSource) to approximate
    MATLAB's camlight + gouraud lighting
  - Anti-aliased, high-DPI PNG output (300 DPI)
  - Larger figure size and tighter layout
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D projection)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

# --------------------------------------------------------------------------
# Configuration definitions
# --------------------------------------------------------------------------
d2r = np.pi / 180.0
beta_pyr = 54.74 * d2r

configs = [
    (
        "Rooftop",
        np.array([
            [0, np.sin(np.pi / 4), np.cos(np.pi / 4)],
            [0, -np.sin(np.pi / 4), np.cos(np.pi / 4)],
            [0, np.sin(np.pi / 4), np.cos(np.pi / 4)],
            [0, -np.sin(np.pi / 4), np.cos(np.pi / 4)],
        ]),
        "independent",
    ),
    (
        "Box",
        np.array([
            [1, 0, 0],
            [-1, 0, 0],
            [0, 1, 0],
            [0, -1, 0],
        ]),
        "independent",
    ),
    (
        "Pyramid",
        np.array([
            [np.sin(beta_pyr), 0, np.cos(beta_pyr)],
            [0, np.sin(beta_pyr), np.cos(beta_pyr)],
            [-np.sin(beta_pyr), 0, np.cos(beta_pyr)],
            [0, -np.sin(beta_pyr), np.cos(beta_pyr)],
        ]),
        "independent",
    ),
    (
        "Scissor Pair (6 CMGs)",
        None,
        "scissor",
    ),
]

# --------------------------------------------------------------------------
# Higher-resolution spherical sampling grid
# --------------------------------------------------------------------------
N_AZ, N_EL = 200, 100
az = np.linspace(0, 2 * np.pi, N_AZ)
el = np.linspace(-np.pi / 2, np.pi / 2, N_EL)
AZ, EL = np.meshgrid(az, el)  # shape (N_EL, N_AZ), matches MATLAB meshgrid(az,el)

ux = np.cos(EL) * np.cos(AZ)
uy = np.cos(EL) * np.sin(AZ)
uz = np.sin(EL)

ls = LightSource(azdeg=315, altdeg=45)
CMAP = cm.turbo


def projected_surface(gx, gy, gz, scale=1.0):
    """Project the unit sphere onto the plane orthogonal to gimbal axis g,
    re-normalize, and scale (mirrors the MATLAB projection step)."""
    udotg = ux * gx + uy * gy + uz * gz
    px = ux - udotg * gx
    py = uy - udotg * gy
    pz = uz - udotg * gz
    norm_p = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    norm_p[norm_p == 0] = 1e-6
    return scale * px / norm_p, scale * py / norm_p, scale * pz / norm_p


def sign_combinations(n):
    """All +/-1 sign combinations for n gimbals, as MATLAB's dec2bin(0:2^n-1)."""
    bits = ((np.arange(2 ** n)[:, None] >> np.arange(n - 1, -1, -1)) & 1)
    return bits * 2 - 1


def plot_internal_surfaces(ax, hx_list, hy_list, hz_list):
    """Plot each gimbal-sign-combination surface semi-transparently."""
    for hx, hy, hz in zip(hx_list, hy_list, hz_list):
        norm = plt.Normalize(vmin=-1.5, vmax=1.5)
        colors = CMAP(norm(hz))
        rgb = ls.shade_rgb(colors[..., :3], hz, blend_mode="soft")
        ax.plot_surface(
            hx, hy, hz,
            facecolors=rgb,
            rstride=2, cstride=2,
            linewidth=0, antialiased=True,
            shade=False, alpha=0.18,
        )


def plot_external_hull(ax, hx_all, hy_all, hz_all):
    """Compute and plot the convex hull as a clean shaded manifold."""
    pts = np.column_stack([hx_all, hy_all, hz_all])
    hull = ConvexHull(pts)
    verts = pts[hull.simplices]

    # Color each triangular facet by the mean z (height) of its vertices
    facet_z = verts[:, :, 2].mean(axis=1)
    norm = plt.Normalize(vmin=facet_z.min(), vmax=facet_z.max())
    facecolors = CMAP(norm(facet_z))

    poly = Poly3DCollection(verts, facecolor=facecolors, edgecolor="none",
                             alpha=0.97, linewidths=0)
    poly.set_zsort("average")
    ax.add_collection3d(poly)

    lim = np.max(np.abs(pts)) * 1.05
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)


def style_axis(ax, title):
    ax.set_xlabel("$h_x$", fontsize=11)
    ax.set_ylabel("$h_y$", fontsize=11)
    ax.set_zlabel("$h_z$", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.view_init(elev=30, azim=45)
    ax.xaxis.pane.set_alpha(0.05)
    ax.yaxis.pane.set_alpha(0.05)
    ax.zaxis.pane.set_alpha(0.05)
    ax.grid(True, alpha=0.15)
    ax.set_box_aspect([1, 1, 1])


# --------------------------------------------------------------------------
# Main loop: build one figure per configuration
# --------------------------------------------------------------------------
output_files = []

for config_name, g, ctype in configs:
    fig = plt.figure(figsize=(16, 8), dpi=300)
    ax_ext = fig.add_subplot(1, 2, 1, projection="3d")
    ax_int = fig.add_subplot(1, 2, 2, projection="3d")

    hx_all, hy_all, hz_all = [], [], []
    hx_list, hy_list, hz_list = [], [], []

    if ctype == "independent":
        h_base = [projected_surface(*g[i]) for i in range(4)]
        eps = sign_combinations(4)  # 16 x 4

        for row in eps:
            hx = sum(row[i] * h_base[i][0] for i in range(4))
            hy = sum(row[i] * h_base[i][1] for i in range(4))
            hz = sum(row[i] * h_base[i][2] for i in range(4))

            hx_all.append(hx.ravel())
            hy_all.append(hy.ravel())
            hz_all.append(hz.ravel())
            hx_list.append(hx); hy_list.append(hy); hz_list.append(hz)

    elif ctype == "scissor":
        beta_s = 45 * d2r
        gs = np.array([
            [np.sin(beta_s), 0, np.cos(beta_s)],
            [0, np.sin(beta_s), np.cos(beta_s)],
            [-np.sin(beta_s), 0, np.cos(beta_s)],
        ])
        h_base = [projected_surface(*gs[i], scale=2.0) for i in range(3)]
        eps = sign_combinations(3)  # 8 x 3

        for row in eps:
            hx = sum(row[i] * h_base[i][0] for i in range(3))
            hy = sum(row[i] * h_base[i][1] for i in range(3))
            hz = sum(row[i] * h_base[i][2] for i in range(3))

            hx_all.append(hx.ravel())
            hy_all.append(hy.ravel())
            hz_all.append(hz.ravel())
            hx_list.append(hx); hy_list.append(hy); hz_list.append(hz)

    hx_all = np.concatenate(hx_all)
    hy_all = np.concatenate(hy_all)
    hz_all = np.concatenate(hz_all)

    plot_external_hull(ax_ext, hx_all, hy_all, hz_all)
    style_axis(ax_ext, f"{config_name} - External Envelope")

    plot_internal_surfaces(ax_int, hx_list, hy_list, hz_list)
    lim = np.max(np.abs(np.column_stack([hx_all, hy_all, hz_all]))) * 1.05
    ax_int.set_xlim(-lim, lim); ax_int.set_ylim(-lim, lim); ax_int.set_zlim(-lim, lim)
    style_axis(ax_int, f"{config_name} - Internal Singularities")

    fig.tight_layout()

    fname = f"/mnt/user-data/outputs/cmg_{config_name.split()[0].lower()}.png"
    fig.savefig(fname, dpi=300, bbox_inches="tight", facecolor="white")
    output_files.append(fname)
    plt.close(fig)
    print(f"Saved: {fname}")

print("\nAll figures generated:", output_files)