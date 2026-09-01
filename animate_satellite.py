"""
animate_satellite.py
====================
3D CubeSat attitude animation with HCMG cluster.

USAGE
-----
  python animate_satellite.py                        # uses built-in demo
  python animate_satellite.py my_simulation.csv      # your MATLAB export
  python animate_satellite.py my_simulation.csv 3    # skip every 3 rows (faster)

CSV FORMAT
----------
Required columns:
  t          – time (s)
  Ex,Ey,Ez,n – quaternion components (scalar-last, body->ECI)
  rx,ry,rz   – ECI position (m)  — used for nadir arrow direction

Optional columns (HCMG):
  delta1,delta2,delta3,delta4  – gimbal angles (rad), one per CMG
  omega_fw                     – flywheel spin speed (rad/s)

Export from MATLAB:
  T = table(t', Ex', Ey', Ez', n', rx', ry', rz', ...
            delta1', delta2', delta3', delta4', omega_fw', ...
            'VariableNames', {'t','Ex','Ey','Ez','n', ...
                              'rx','ry','rz', ...
                              'delta1','delta2','delta3','delta4','omega_fw'});
  writetable(T, 'my_simulation.csv');

BACKEND / WINDOW FIX
--------------------
If you get an error about Qt5 or no window appears:
  pip install PyQt5          # recommended on Windows
  OR
  pip install PyQt6
Tkinter (TkAgg) comes bundled with the standard Python Windows installer.
If it is missing, re-run the Python installer and tick "tcl/tk and IDLE".

SAVE TO FILE (instead of interactive window)
--------------------------------------------
Set SAVE_ANIMATION = True below and choose a filename.
Requires ffmpeg for .mp4:  https://ffmpeg.org/download.html
Or use .gif (no ffmpeg needed, but large file).
"""

import sys
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
#  BACKEND SELECTION  — edit FORCE_BACKEND to pin one, or leave '' for auto
# ─────────────────────────────────────────────────────────────────────────────
FORCE_BACKEND  = ''        # e.g. 'TkAgg', 'Qt5Agg', 'Qt6Agg' — or '' for auto
SAVE_ANIMATION = False     # True -> save to file instead of showing window
SAVE_FILENAME  = 'satellite_animation.gif'   # .mp4 needs ffmpeg; .gif always works

import matplotlib
if FORCE_BACKEND:
    matplotlib.use(FORCE_BACKEND)
else:
    # Try in order; TkAgg first because Tkinter ships with Windows Python
    _candidates = ['TkAgg', 'Qt5Agg', 'Qt6Agg', 'WXAgg', 'GTK3Agg', 'MacOSX']
    _dep_check  = {
        'TkAgg':   'tkinter',
        'Qt5Agg':  'PyQt5',
        'Qt6Agg':  'PyQt6',
        'WXAgg':   'wx',
        'GTK3Agg': 'gi',
        'MacOSX':  None,
    }
    _chosen = None
    for _b in _candidates:
        try:
            dep = _dep_check.get(_b)
            if dep:
                import importlib
                importlib.import_module(dep)
            matplotlib.use(_b)
            _chosen = _b
            break
        except Exception:
            continue
    if _chosen is None:
        matplotlib.use('Agg')
        print("WARNING: No interactive backend found.")
        print("  pip install PyQt5   <- recommended fix on Windows")
        print("  Setting SAVE_ANIMATION=True automatically.")
        SAVE_ANIMATION = True

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ─────────────────────────────────────────────────────────────────────────────
#  SATELLITE GEOMETRY  (edit freely)
# ─────────────────────────────────────────────────────────────────────────────
SAT = dict(
    # 12U CubeSat: 2U x 2U x 3U  ->  0.2 m x 0.2 m x 0.3 m
    Lx = 0.20,    # body dimension along body-X  (m)
    Ly = 0.20,    # body dimension along body-Y  (m)
    Lz = 0.30,    # body dimension along body-Z  (m)

    # Centre of Mass offset from geometric centre (m)
    CoM = np.array([0.005, -0.003, 0.010]),

    # Centre of Pressure offset from geometric centre (m)
    CoP = np.array([0.001,  0.002, 0.015]),
)

# HCMG cluster — 4-CMG pyramid
HCMG = dict(
    skew_angle        = np.deg2rad(54.74),  # pyramid half-angle
    flywheel_radius   = 0.030,              # disk radius (m)
    flywheel_thickness= 0.008,              # disk thickness (m)
    gimbal_arm_length = 0.025,              # arm half-length (m)
    cluster_scale     = 0.55,              # scale to fit inside satellite
)

# ─────────────────────────────────────────────────────────────────────────────
#  ANIMATION SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
ARROW_LEN   = 0.28
BG_COLOR    = '#0d0d0d'
INTERVAL_MS = 30          # ms between frames

# ─────────────────────────────────────────────────────────────────────────────
#  MATH HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def quat_to_rotm(q):
    """q = [Ex, Ey, Ez, n] scalar-last -> 3x3 rotation matrix (body->ECI)."""
    Ex, Ey, Ez, n = q
    return np.array([
        [1-2*(Ey**2+Ez**2),  2*(Ex*Ey - n*Ez),   2*(Ex*Ez + n*Ey)],
        [2*(Ex*Ey + n*Ez),   1-2*(Ex**2+Ez**2),   2*(Ey*Ez - n*Ex)],
        [2*(Ex*Ez - n*Ey),   2*(Ey*Ez + n*Ex),   1-2*(Ex**2+Ey**2)],
    ])

def rodrigues(axis, angle):
    """Rotation matrix via Rodrigues formula."""
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle)*K + (1 - np.cos(angle))*(K @ K)

def make_disk_verts(radius, thickness, n_pts=24):
    """Flat cylinder (disk) faces, axis along local Z. Returns list of (k,3) arrays."""
    theta = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    top = np.stack([radius*cos_t, radius*sin_t, np.full(n_pts, thickness/2)], axis=1)
    bot = np.stack([radius*cos_t, radius*sin_t, np.full(n_pts, -thickness/2)], axis=1)
    tc  = np.array([0., 0.,  thickness/2])
    bc  = np.array([0., 0., -thickness/2])
    polys = []
    for i in range(n_pts):
        j = (i+1) % n_pts
        polys.append(np.array([top[i], top[j], bot[j], bot[i]]))  # side
        polys.append(np.array([tc, top[i], top[j]]))               # top cap
        polys.append(np.array([bc, bot[j], bot[i]]))               # bot cap
    return polys

# ─────────────────────────────────────────────────────────────────────────────
#  SATELLITE BOX GEOMETRY
# ─────────────────────────────────────────────────────────────────────────────
def build_sat_geometry(sat):
    hx, hy, hz = sat['Lx']/2, sat['Ly']/2, sat['Lz']/2
    gc = np.array([[-hx,-hy,-hz],[ hx,-hy,-hz],[ hx, hy,-hz],[-hx, hy,-hz],
                   [-hx,-hy, hz],[ hx,-hy, hz],[ hx, hy, hz],[-hx, hy, hz]])
    corners = gc - sat['CoM']        # shift origin to CoM

    face_idx = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
    face_verts  = [corners[idx] for idx in face_idx]
    face_colors = ['#3a3a88','#cc5500','#777777','#777777','#999999','#999999']
    return face_verts, face_colors, gc, corners

# ─────────────────────────────────────────────────────────────────────────────
#  HCMG CLUSTER
# ─────────────────────────────────────────────────────────────────────────────
def build_cmg_bases(hcmg):
    """Returns list of 4 CMG dicts (azimuth, nominal spin axis, gimbal axis)."""
    beta     = hcmg['skew_angle']
    azimuths = [0, np.pi/2, np.pi, 3*np.pi/2]
    cmgs = []
    for az in azimuths:
        spin_nom = np.array([np.sin(beta)*np.cos(az),
                             np.sin(beta)*np.sin(az),
                             np.cos(beta)])
        gimbal_ax = np.array([-np.sin(az), np.cos(az), 0.0])
        cmgs.append({'azimuth': az, 'spin_axis_nominal': spin_nom,
                     'gimbal_axis': gimbal_ax})
    return cmgs

# ─────────────────────────────────────────────────────────────────────────────
#  MAIN ANIMATION
# ─────────────────────────────────────────────────────────────────────────────
def run_animation(df, skip=1):
    df = df.iloc[::skip].reset_index(drop=True)

    has_hcmg = all(c in df.columns for c in ['delta1','delta2','delta3','delta4'])
    has_fw   = 'omega_fw' in df.columns

    sat       = SAT
    hcmg_cfg  = HCMG
    face_verts, face_colors, corners_gc, corners_com = build_sat_geometry(sat)
    cmg_bases = build_cmg_bases(hcmg_cfg)

    sc    = hcmg_cfg['cluster_scale']
    r_fw  = hcmg_cfg['flywheel_radius']    * sc
    t_fw  = hcmg_cfg['flywheel_thickness'] * sc
    arm   = hcmg_cfg['gimbal_arm_length']  * sc
    beta  = hcmg_cfg['skew_angle']

    disk_polys_local = make_disk_verts(r_fw, t_fw)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(11, 8), facecolor=BG_COLOR)
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors='#666666', labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#2a2a2a')
    for spine in ax.spines.values():
        spine.set_color('#333333')
    ax.xaxis.label.set_color('#aaaaaa'); ax.yaxis.label.set_color('#aaaaaa')
    ax.zaxis.label.set_color('#aaaaaa')
    ax.set_xlabel('X_I', fontsize=8); ax.set_ylabel('Y_I', fontsize=8)
    ax.set_zlabel('Z_I', fontsize=8)
    ax.set_title('12U CubeSat + HCMG Attitude Animation',
                 color='white', fontsize=13, pad=10)

    L   = ARROW_LEN
    lim = 0.48
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)

    # Inertial frame (faint dashed)
    for vec, lbl in [([L,0,0],'X_I'),([0,L,0],'Y_I'),([0,0,L],'Z_I')]:
        ax.quiver(0,0,0,*vec, color='#444444', linewidth=0.8,
                  arrow_length_ratio=0.12, alpha=0.5, linestyle='--')
        ax.text(vec[0]+0.025, vec[1]+0.01, vec[2]+0.01, lbl,
                color='#666666', fontsize=7)

    # ── Satellite body ────────────────────────────────────────────────────────
    body_patch = Poly3DCollection(face_verts, facecolors=face_colors,
                                  edgecolors='white', linewidths=0.7, alpha=0.82)
    ax.add_collection3d(body_patch)

    # ── Body axes (stored as list so we can replace quivers) ─────────────────
    # matplotlib Quiver3D segments can be updated in-place
    qx = ax.quiver(0,0,0,L,0,0, color='#ff4444', linewidth=2.2,
                   arrow_length_ratio=0.12)
    qy = ax.quiver(0,0,0,0,L,0, color='#44ff44', linewidth=2.2,
                   arrow_length_ratio=0.12)
    qz = ax.quiver(0,0,0,0,0,L, color='#4488ff', linewidth=2.2,
                   arrow_length_ratio=0.12)
    tx = ax.text(L+0.03, 0, 0, 'X_B', color='#ff4444', fontsize=9, fontweight='bold')
    ty = ax.text(0, L+0.03, 0, 'Y_B', color='#44ff44', fontsize=9, fontweight='bold')
    tz = ax.text(0, 0, L+0.03, 'Z_B', color='#4488ff', fontsize=9, fontweight='bold')

    # ── Nadir arrow ───────────────────────────────────────────────────────────
    qn = ax.quiver(0,0,0, 0,0,-L, color='cyan', linewidth=2.5,
                   arrow_length_ratio=0.14)
    tn = ax.text(0, 0, -L-0.06, 'Earth', color='cyan', fontsize=9,
                 fontweight='bold')

    # ── CoM marker (fixed at origin — body rotates around it) ────────────────
    ax.plot([0],[0],[0], marker='+', color='yellow', markersize=14,
            markeredgewidth=2.5, zorder=10, linestyle='none')
    ax.text(0.025, 0.0, 0.025, 'CoM', color='yellow', fontsize=8,
            fontweight='bold')

    # ── CoP marker (moves with body) ─────────────────────────────────────────
    cop_body = sat['CoP'] - sat['CoM']   # CoP relative to CoM in body frame
    cop_pt,  = ax.plot([cop_body[0]], [cop_body[1]], [cop_body[2]],
                       marker='x', color='magenta', markersize=10,
                       markeredgewidth=2.0, zorder=10, linestyle='none')
    cop_lbl   = ax.text(cop_body[0]+0.025, cop_body[1], cop_body[2]+0.025,
                        'CoP', color='magenta', fontsize=8)

    # Geometric centre marker (moves with body)
    gc_body = -sat['CoM']
    gc_pt,  = ax.plot([gc_body[0]], [gc_body[1]], [gc_body[2]],
                      marker='+', color='#aaaaaa', markersize=8,
                      markeredgewidth=1.5, zorder=9, linestyle='none', alpha=0.6)
    gc_lbl   = ax.text(gc_body[0]+0.02, gc_body[1], gc_body[2]+0.02,
                       'GC', color='#888888', fontsize=7, alpha=0.8)

    # ── HCMG cluster patches ─────────────────────────────────────────────────
    cmg_fw_patches = []
    cmg_arm_lines  = []

    if has_hcmg:
        for i in range(4):
            fw_patch = Poly3DCollection(disk_polys_local,
                                        facecolors='#b8860b',
                                        edgecolors='#ffe066',
                                        linewidths=0.3, alpha=0.85)
            ax.add_collection3d(fw_patch)
            cmg_fw_patches.append(fw_patch)

            gline, = ax.plot([],[],[], color='#ddaa00', linewidth=2.0,
                             alpha=0.9, solid_capstyle='round')
            cmg_arm_lines.append(gline)

    # ── HUD text ─────────────────────────────────────────────────────────────
    t_txt  = ax.text2D(0.02, 0.97, '', transform=ax.transAxes,
                       color='white', fontsize=10)
    al_txt = ax.text2D(0.02, 0.92, '', transform=ax.transAxes,
                       color='white', fontsize=10, fontweight='bold')
    fw_txt = ax.text2D(0.02, 0.87, '', transform=ax.transAxes,
                       color='#aaaaaa', fontsize=9)
    ax.text2D(0.02, 0.04,
              f"CoM: [{sat['CoM'][0]*1e3:.1f}, {sat['CoM'][1]*1e3:.1f}, "
              f"{sat['CoM'][2]*1e3:.1f}] mm\n"
              f"CoP: [{sat['CoP'][0]*1e3:.1f}, {sat['CoP'][1]*1e3:.1f}, "
              f"{sat['CoP'][2]*1e3:.1f}] mm",
              transform=ax.transAxes, color='#bbbbbb', fontsize=8,
              verticalalignment='bottom')

    # ── UPDATE ────────────────────────────────────────────────────────────────
    def update(frame):
        row = df.iloc[frame]

        # --- Attitude ---------------------------------------------------------
        q  = np.array([row['Ex'], row['Ey'], row['Ez'], row['n']], dtype=float)
        q /= np.linalg.norm(q) + 1e-15
        R  = quat_to_rotm(q)      # body -> ECI

        # --- Body faces -------------------------------------------------------
        body_patch.set_verts([(R @ fv.T).T for fv in face_verts])

        # --- Body axes --------------------------------------------------------
        xb = R[:,0]; yb = R[:,1]; zb = R[:,2]
        qx.set_segments([[[0,0,0], (L*xb).tolist()]])
        qy.set_segments([[[0,0,0], (L*yb).tolist()]])
        qz.set_segments([[[0,0,0], (L*zb).tolist()]])
        tx.set_position_3d((L+0.03)*xb)
        ty.set_position_3d((L+0.03)*yb)
        tz.set_position_3d((L+0.03)*zb)

        # --- Nadir arrow ------------------------------------------------------
        r_eci = np.array([row['rx'], row['ry'], row['rz']], dtype=float)
        nr    = np.linalg.norm(r_eci)
        nadir_I = -r_eci/nr if nr > 1e3 else np.array([0.,0.,-1.])
        qn.set_segments([[[0,0,0], (L*nadir_I).tolist()]])
        tn.set_position_3d((L+0.07)*nadir_I)

        # --- CoP / GC markers (rotate with body) ------------------------------
        cop_eci = R @ cop_body
        cop_pt.set_data_3d([cop_eci[0]], [cop_eci[1]], [cop_eci[2]])
        cop_lbl.set_position_3d(cop_eci + np.array([0.025, 0, 0.025]))

        gc_eci = R @ gc_body
        gc_pt.set_data_3d([gc_eci[0]], [gc_eci[1]], [gc_eci[2]])
        gc_lbl.set_position_3d(gc_eci + np.array([0.02, 0, 0.02]))

        # --- Alignment indicator (Z_body -> nadir) ---------------------------
        cos_a     = float(np.clip(np.dot(zb, nadir_I), -1, 1))
        angle_deg = np.degrees(np.arccos(cos_a))
        if angle_deg < 10:
            al_txt.set_color('#00ff66');  al_txt.set_text(f'Z->Earth: {angle_deg:.1f} deg  OK')
        elif angle_deg < 30:
            al_txt.set_color('#ffff00');  al_txt.set_text(f'Z->Earth: {angle_deg:.1f} deg  ~')
        else:
            al_txt.set_color('#ff4444');  al_txt.set_text(f'Z->Earth: {angle_deg:.1f} deg  X')

        t_txt.set_text(f't = {row["t"]:.1f} s')

        # --- HCMG cluster -----------------------------------------------------
        if has_hcmg:
            deltas = [float(row.get('delta1', 0)), float(row.get('delta2', 0)),
                      float(row.get('delta3', 0)), float(row.get('delta4', 0))]
            fw_speed = float(row['omega_fw']) if has_fw else 500.0
            if has_fw:
                fw_txt.set_text(f'omega_fw = {fw_speed:.0f} rad/s')

            # Visual flywheel self-spin accumulator
            dt_total = float(row['t'])
            fw_phi   = (fw_speed * dt_total * 0.005) % (2*np.pi)   # visual only

            for i, (cmg, fw_patch, gline) in enumerate(
                    zip(cmg_bases, cmg_fw_patches, cmg_arm_lines)):

                delta_i  = deltas[i]
                g_ax     = cmg['gimbal_axis']        # unit vector, body frame
                spin_nom = cmg['spin_axis_nominal']  # nominal flywheel axis

                # 1) Gimbal rotation -> actual flywheel axis
                R_gim  = rodrigues(g_ax, delta_i)
                fw_ax  = R_gim @ spin_nom

                # 2) Flywheel self-spin around fw_ax
                R_spin = rodrigues(fw_ax, fw_phi)

                # 3) Combined local->body
                R_cmg  = R_spin @ R_gim

                # CMG centre position in body frame (tip of pyramid)
                az  = cmg['azimuth']
                c_body = np.array([np.sin(beta)*np.cos(az),
                                   np.sin(beta)*np.sin(az),
                                   np.cos(beta)]) * 0.11 * sc

                # Rotate disk verts: local->body->ECI
                new_polys = []
                for poly in disk_polys_local:
                    pts_eci = (R @ (R_cmg @ poly.T + c_body[:,None])).T
                    new_polys.append(pts_eci)
                fw_patch.set_verts(new_polys)

                # Gimbal arm endpoints (body -> ECI)
                p0 = R @ (c_body - arm * g_ax)
                p1 = R @ (c_body + arm * g_ax)
                gline.set_data_3d([p0[0],p1[0]], [p0[1],p1[1]], [p0[2],p1[2]])

        return ([body_patch, qx, qy, qz, tx, ty, tz,
                 qn, tn, cop_pt, cop_lbl, gc_pt, gc_lbl,
                 t_txt, al_txt, fw_txt]
                + cmg_fw_patches + cmg_arm_lines)

    anim = animation.FuncAnimation(
        fig, update,
        frames=len(df),
        interval=INTERVAL_MS,
        blit=False,
    )

    if SAVE_ANIMATION:
        print(f"Saving animation to '{SAVE_FILENAME}' ...")
        if SAVE_FILENAME.endswith('.gif'):
            writer = animation.PillowWriter(fps=30)
        else:
            writer = animation.FFMpegWriter(fps=30, bitrate=1800)
        anim.save(SAVE_FILENAME, writer=writer)
        print("Saved.")
    else:
        plt.tight_layout()
        plt.show()

    return anim


# ─────────────────────────────────────────────────────────────────────────────
#  ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    csv_path = sys.argv[1] if len(sys.argv) > 1 else None
    skip     = int(sys.argv[2]) if len(sys.argv) > 2 else 1

    if csv_path is None:
        print("No CSV supplied — running built-in demo.")
        N = 400
        t = np.linspace(0, 100, N)
        angle = np.where(t < 60, np.pi/3*(t/60)**2, np.pi/3)
        Ex = 0.002*np.sin(0.3*t)
        Ey = np.sin(angle/2)
        Ez = 0.001*np.cos(0.5*t)
        n  = np.cos(angle/2)
        norms = np.sqrt(Ex**2 + Ey**2 + Ez**2 + n**2)
        Ex /= norms; Ey /= norms; Ez /= norms; n /= norms

        R_orb = 6878e3
        om    = np.sqrt(3.986e14 / R_orb**3)
        df = pd.DataFrame({
            't': t,
            'Ex': Ex, 'Ey': Ey, 'Ez': Ez, 'n': n,
            'rx': R_orb*np.cos(om*t),
            'ry': R_orb*np.sin(om*t),
            'rz': R_orb*0.05*np.sin(0.5*om*t),
            'delta1':  0.5*np.sin(0.05*t),
            'delta2': -0.5*np.sin(0.05*t + np.pi/2),
            'delta3':  0.3*np.cos(0.05*t + np.pi/4),
            'delta4': -0.3*np.cos(0.05*t - np.pi/4),
            'omega_fw': 500 + 30*np.sin(0.08*t),
        })
    else:
        print(f"Loading: {csv_path}")
        df = pd.read_csv(csv_path)
        required = ['t', 'Ex', 'Ey', 'Ez', 'n', 'rx', 'ry', 'rz']
        missing  = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"CSV is missing required columns: {missing}")
        print(f"  {len(df)} rows loaded  (skip={skip})")

    run_animation(df, skip=skip)