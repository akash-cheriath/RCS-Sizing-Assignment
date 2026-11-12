"""plot_vehicle.py - figures of the RCS sizing study: the preferred layout on the vehicle, all architectures in 3D,
and a requirement-compliance grid.

Reads config/config.cfg (with its include files) and the CSV files the C++ tool wrote to results/, and writes
three images to figures/:

  preferred_layout.png    1 x 3: (a) side view of the vehicle with its dimensions and CoGs, (b) side view of the
                          rank-1 layout with its ring station and exhaust directions, (c) the same layout in 3D
  all_architectures.png   every architecture in 3D (thruster positions and exhaust directions), in rank order
  requirements_table.png  requirements R-01 ... R-16 (rows) against every architecture in rank order (columns):
                          green = met, yellow = partly met, red = not met. Thrust-dependent requirements are checked
                          with the baseline thruster (thrust, min_pulse_time of baseline_thruster.cfg)

Usage:   python plot_vehicle.py                  (project folder = the folder above this script)
         python plot_vehicle.py <project folder>
Needs:   matplotlib, numpy

Geometry (vehicle.cfg). x = 0 is the heat-shield max-diameter plane, which is also the top of the service module.
  capsule_height_in_combination   capsule height above the service module: the nose (capsule top face) is at this x
  capsule_height                  whole capsule; the difference is the heat shield, which sits inside the SM
"""
import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402

ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")

# ---- colours: categorical slots of the reference palette (one hue per thruster function), status colours for
#      the requirement table (always with an icon and a word), neutral chrome
SLOT = ["#2a78d6", "#eb6834", "#1baf7a"]
STATUS = {"met": ("#0ca30c", "✓ Met"), "partly": ("#fab219", "◐ Partly met"), "not": ("#d03b3b", "✗ Not met")}
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
CAP_FILL, SM_FILL, EDGE = "#eceae4", "#f4f3ef", "#8f8d86"
plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"], "font.size": 8.5,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9.5,
    "xtick.color": MUTED, "ytick.color": INK2, "axes.grid": False,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8,
})

# thruster type -> colour, legend label. Colour follows the function: +X force, -X force, roll
TYPE = {
    "AFT": (SLOT[0], "AFT: exhaust aft, force +X"),
    "FWD": (SLOT[1], "FWD: exhaust forward, force -X"),
    "NOSE": (SLOT[1], "NOSE: exhaust forward, force -X"),
    "ROLLP": (SLOT[2], "ROLL+ / ROLL-: exhaust sideways (roll)"),
    "ROLLN": (SLOT[2], None),
}
SHORT = {"AFT": "AFT", "FWD": "FWD", "NOSE": "NOSE", "ROLLP": "ROLL+", "ROLLN": "ROLL-"}
ARROW = 0.6                                          # [m] drawn length of an exhaust arrow (side view)


# ------------------------------------------------------------------------------------------ input
def read_config(path, cfg=None):
    """key = value  # comment;  'include = file' reads that file too (same folder)."""
    cfg = {} if cfg is None else cfg
    folder = os.path.dirname(path)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.split("#", 1)[0]
            if "=" not in line:
                continue
            k, v = (s.strip() for s in line.split("=", 1))
            if k == "include":
                read_config(os.path.join(folder, v), cfg)
            elif k:
                cfg[k] = v
    return cfg


def rows(name):
    with open(os.path.join(RES, name), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


cfg = read_config(os.path.join(ROOT, "config", "config.cfg"))


def num(key, default=None):
    if key in cfg:
        return float(cfg[key])
    if default is None:
        raise KeyError(f"config key '{key}' not found")
    return default


# ---- geometry (body frame: origin on the heat-shield max-diameter plane = top of the SM, +X to the nose)
R0 = 0.5 * num("capsule_base_diameter")
HALF = math.radians(num("sidewall_half_angle"))
H_CAP = num("capsule_height")
H_TOP = min(num("capsule_height_in_combination", num("capsule_height_nyx", H_CAP)), H_CAP)   # nose = capsule top
HS = H_CAP - H_TOP                                   # heat-shield depth below x = 0 (inside the SM)
H_SM = num("service_module_height")
R_SM = 0.5 * num("service_module_diameter")
T_MIN = num("min_pulse_time")                        # [ms] baseline thruster
F_BASE = num("thrust")
STUCK_MIN = num("stuck_on_ratio_min", 1.5)


def wall_r(x):
    return R0 - x * math.tan(HALF)


def d(x):
    """distance from the nose (capsule top face) [m]"""
    return H_TOP - x


def heat_shield_x(z):
    """x of the heat-shield surface: spherical cap of depth HS through the capsule base edge"""
    z = np.asarray(z, dtype=float)
    if HS <= 1e-9:
        return np.zeros_like(z)
    rs = (R0 ** 2 + HS ** 2) / (2 * HS)
    return -HS + (rs - np.sqrt(np.maximum(rs ** 2 - z ** 2, 0)))


R_TOP = wall_r(H_TOP)

# ---- CoGs: combined and capsule from the config, service module from the mass balance of the C++ check
mp = rows("mass_properties_check.csv")[0]
COGS = [
    ("Combined", num("combined_mass"), np.array([num("combined_cog_x"), num("combined_cog_y", 0.0), num("combined_cog_z", 0.0)]),
     dict(marker="o", ms=6.5, mfc=INK, mec=SURFACE, mew=0.8)),
    ("Capsule", num("capsule_mass"), np.array([num("capsule_cog_x"), num("capsule_cog_y", 0.0), num("capsule_cog_z", 0.0)]),
     dict(marker="D", ms=5.5, mfc=SURFACE, mec=INK, mew=1.2)),
    ("Service module", float(mp["sm_mass_kg"]), np.array([float(mp["sm_cg_x_m"]), float(mp["sm_cg_y_m"]), float(mp["sm_cg_z_m"])]),
     dict(marker="s", ms=5.5, mfc=SURFACE, mec=INK2, mew=1.2)),
]

# ---- results
ranking = sorted(rows("architecture_ranking.csv"), key=lambda r: int(r["rank"]))
n_arch = len(ranking)
summary = {r["layout"]: r for r in rows("layout_summary.csv")}
cmp_rows = rows("comparison.csv")
thr_rows = rows("layout_thrusters.csv")
drivers = {}
for r in rows("thrust_drivers.csv"):
    drivers.setdefault(r["layout"], []).append(r)


def col(t, a, b):
    return float(t[a] if a in t else t[b])


def thrusters(layout):
    out = []
    for t in thr_rows:
        if t["layout"] != layout:
            continue
        f = np.array([col(t, "dir_x", "dx"), col(t, "dir_y", "dy"), col(t, "dir_z", "dz")])   # unit FORCE direction
        out.append(dict(name=t["name"], type=t["type"], pod=int(t["pod"]),
                        pos=np.array([float(t["x"]), float(t["y"]), float(t["z"])]), exhaust=-f / np.linalg.norm(f)))
    return out


def rings(thr):
    """thrusters grouped by axial station, nose first: [(x, [thrusters])]"""
    g = {}
    for t in thr:
        g.setdefault(round(float(t["pos"][0]), 4), []).append(t)
    return sorted(g.items(), key=lambda kv: -kv[0])


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name), dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", os.path.join("figures", name))


# ------------------------------------------------------------------------------------------ drawing helpers
def draw_body(ax):
    """side view: horizontal = body Z, vertical = distance from the nose (downwards)"""
    ax.add_patch(Rectangle((-R_SM, d(0)), 2 * R_SM, H_SM, fc=SM_FILL, ec=EDGE, lw=1.0, zorder=1))
    ax.add_patch(Polygon([(-R0, d(0)), (-R_TOP, d(H_TOP)), (R_TOP, d(H_TOP)), (R0, d(0))],
                         closed=True, fc=CAP_FILL, ec=EDGE, lw=1.0, zorder=1))
    if HS > 1e-6:                                    # heat shield, hidden inside the service module
        zz = np.linspace(-R0, R0, 80)
        ax.plot(zz, d(heat_shield_x(zz)), color=EDGE, lw=0.9, ls=(0, (3, 2)), zorder=2)
    ax.plot([0, 0], [d(H_TOP) - 0.25, d(-H_SM) + 0.25], color=AXIS, lw=0.7, ls=(0, (8, 3, 2, 3)), zorder=1)


def dim(ax, p0, p1, text, offset=0.0, rotation=0, fs=7.5):
    """dimension line between p0 and p1 (data coordinates) with the value written beside it"""
    ax.annotate("", xy=p0, xytext=p1, arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=0.8, mutation_scale=7,
                                                      shrinkA=0, shrinkB=0))
    mx, my = 0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])
    if rotation:
        ax.text(mx - offset, my, text, rotation=90, ha="right", va="center", fontsize=fs, color=INK)
    else:
        ax.text(mx, my - offset, text, ha="center", va="bottom", fontsize=fs, color=INK)


def ext(ax, z0, z1, dd):
    """thin extension line from the drawing to a dimension line"""
    ax.plot([z0, z1], [dd, dd], color=AXIS, lw=0.6, zorder=0)


def side_axes(ax, zlim):
    ax.set_xlim(*zlim)
    ax.set_ylim(d(-H_SM) + 0.5, d(H_TOP) - 0.7)      # distance from the nose grows downwards
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def body_3d(ax):
    th = np.linspace(0, 2 * math.pi, 73)
    T, X = np.meshgrid(th, np.array([0.0, H_TOP]))
    R = R0 - X * math.tan(HALF)
    ax.plot_surface(R * np.cos(T), R * np.sin(T), X, color="#b9b6ad", alpha=0.16, linewidth=0, shade=False)
    T3, X3 = np.meshgrid(th, np.array([-H_SM, 0.0]))
    ax.plot_surface(R_SM * np.cos(T3), R_SM * np.sin(T3), X3, color="#d6d4cc", alpha=0.10, linewidth=0, shade=False)
    for x, r in ((0.0, R0), (H_TOP, R_TOP), (-H_SM, R_SM)):
        ax.plot(r * np.cos(th), r * np.sin(th), np.full_like(th, x), color=EDGE, lw=0.6)
    for ang in np.radians([0, 90, 180, 270]):        # a few generator lines for the shape
        ax.plot([R0 * math.cos(ang), R_TOP * math.cos(ang)], [R0 * math.sin(ang), R_TOP * math.sin(ang)], [0, H_TOP],
                color=AXIS, lw=0.4)
        ax.plot([R_SM * math.cos(ang)] * 2, [R_SM * math.sin(ang)] * 2, [-H_SM, 0], color=AXIS, lw=0.4)


def layout_3d(ax, thr, arrow=0.85, lw=1.3):
    """thruster positions and exhaust arrows of one layout on the vehicle; plot axes = (body Y, body Z, body X)"""
    body_3d(ax)
    P = np.array([t["pos"] for t in thr])
    E = np.array([t["exhaust"] for t in thr])
    for typ in ("AFT", "FWD", "NOSE", "ROLLP", "ROLLN"):
        m = np.array([t["type"] == typ for t in thr])
        if m.any():
            ax.quiver(P[m, 1], P[m, 2], P[m, 0], E[m, 1], E[m, 2], E[m, 0], length=arrow, normalize=True,
                      color=TYPE[typ][0], linewidth=lw, arrow_length_ratio=0.3)
    ax.scatter(P[:, 1], P[:, 2], P[:, 0], s=6, color=INK2, depthshade=False)
    lim = max(R0, R_SM) + 0.3
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-H_SM - 1.2, H_TOP + 0.9)            # extra room below: the SM bottom ellipse tilts into view
    ax.set_box_aspect((2 * lim, 2 * lim, H_TOP + H_SM + 2.1), zoom=1.15)
    ax.view_init(elev=16, azim=35)
    ax.set_axis_off()


def exhaust_arrows(ax, thr):
    """side view of the exhaust directions; far-side jets (y < 0) are drawn faded"""
    seen = set()
    for near in (True, False):
        for t in thr:
            if (t["pos"][1] > -1e-6) != near:
                continue
            z0, d0 = t["pos"][2], d(t["pos"][0])
            e = t["exhaust"]
            z1, d1 = z0 + ARROW * e[2], d0 - ARROW * e[0]
            key = tuple(round(v, 3) for v in (z0, d0, z1, d1))
            if key in seen:
                continue
            seen.add(key)
            a = 1.0 if near else 0.35
            ax.annotate("", xy=(z1, d1), xytext=(z0, d0), zorder=5 if near else 3,
                        arrowprops=dict(arrowstyle="-|>", color=TYPE[t["type"]][0], lw=1.5, alpha=a,
                                        mutation_scale=9, shrinkA=0, shrinkB=0))
            ax.plot(z0, d0, ls="", marker="o", ms=3.4, color=INK2, alpha=a, zorder=6 if near else 3)


# ------------------------------------------------------------------------------------------ figure 1: preferred layout
def fig_preferred():
    best = ranking[0]["layout"]
    thr = thrusters(best)
    cogs = [c for c in COGS if c[0] != "Service module"]           # combined and capsule CoGs only
    zlim = (-R0 - 1.4, R0 + 1.3)
    s = 0.62                                         # scale [inch / m]
    w2d = s * (zlim[1] - zlim[0])
    h2d = s * (H_TOP + H_SM + 1.2)
    w3d = 0.8 * h2d
    L, G, T, B = 0.1, 0.25, 0.15, 0.75               # margins [inch]; B leaves room for legend + subcaption
    W, H = L + 2 * w2d + w3d + 2 * G + 0.1, T + h2d + B
    fig = plt.figure(figsize=(W, H))
    x0 = [L, L + w2d + G, L + 2 * w2d + 2 * G]      # left edge of each panel [inch]
    wid = [w2d, w2d, w3d]

    def subcaption(k, text):
        fig.text((x0[k] + 0.5 * wid[k]) / W, 0.12 / H, text, fontsize=9, color=INK, ha="center", va="bottom")

    # (a) vehicle with dimensions and CoGs
    ax = fig.add_axes([x0[0] / W, B / H, w2d / W, h2d / H])
    side_axes(ax, zlim)
    draw_body(ax)
    rmax = max(R0, R_SM)
    z1, z2 = -rmax - 0.45, -rmax - 1.0
    dim(ax, (z1, d(H_TOP)), (z1, d(0)), f"{H_TOP:.2f} m", 0.07, 90)              # capsule in combination
    dim(ax, (z1, d(0)), (z1, d(-H_SM)), f"{H_SM:.2f} m", 0.07, 90)               # service module
    dim(ax, (z2, d(H_TOP)), (z2, d(-HS)), f"{H_CAP:.2f} m", 0.07, 90)            # whole capsule incl. heat shield
    for dd, zf, zt in ((d(H_TOP), -R_TOP, z2), (d(0), -R0, z1), (d(-HS), -0.1, z2), (d(-H_SM), -R_SM, z1)):
        ext(ax, zf, zt - 0.07, dd)
    dim(ax, (-R0, d(0) - 0.55), (R0, d(0) - 0.55), f"Ø {2 * R0:.2f} m", 0.06)    # capsule max diameter
    for zz in (-R0, R0):
        ax.plot([zz, zz], [d(0) - 0.05, d(0) - 0.62], color=AXIS, lw=0.6, zorder=0)
    # CoGs: markers, with their distance from the nose on the right
    zr = R0 + 0.35
    for name, m, p, st in cogs:
        ax.plot(p[2], d(p[0]), ls="", zorder=6, **dict(st, ms=st["ms"] * 0.85))
    for k, (name, m, p, st) in enumerate(cogs):
        zc = zr + 0.45 * k
        dim(ax, (zc, d(H_TOP)), (zc, d(p[0])), f"{d(p[0]):.2f} m", 0.07, 90, fs=7)
        ext(ax, p[2] + 0.08, zc + 0.05, d(p[0]))
    ext(ax, R_TOP, zr + 0.5, d(H_TOP))
    for name, m, p, st in cogs:
        ax.plot([], [], ls="", label=f"CoG {name.lower()}", **dict(st, ms=st["ms"] * 0.85))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=6.5, handletextpad=0.3,
              columnspacing=0.8)
    subcaption(0, "(a) Capsule Model with Service Module")

    # (b) preferred layout, side view
    ax = fig.add_axes([x0[1] / W, B / H, w2d / W, h2d / H])
    side_axes(ax, zlim)
    draw_body(ax)
    exhaust_arrows(ax, thr)
    for k, (x, ts) in enumerate(rings(thr)):
        r = max(math.hypot(t["pos"][1], t["pos"][2]) for t in ts)
        ax.plot([-r, r], [d(x)] * 2, color=MUTED, lw=0.6, ls=(0, (2, 2)), zorder=2)
        zc = R0 + 0.45 + 0.45 * k
        dim(ax, (zc, d(H_TOP)), (zc, d(x)), f"{d(x):.2f} m", 0.07, 90, fs=7)
        ext(ax, r + 0.1, zc + 0.05, d(x))
    ext(ax, R_TOP, R0 + 0.5, d(H_TOP))
    for name, m, p, st in cogs:
        ax.plot(p[2], d(p[0]), ls="", zorder=6, **dict(st, ms=st["ms"] * 0.7))
    h = [plt.Line2D([], [], color=TYPE[t][0], lw=2, marker=">", ms=4) for t in ("AFT", "FWD", "ROLLP")]
    ax.legend(h, ["aft exhaust", "forward exhaust", "roll exhaust"], loc="upper center", bbox_to_anchor=(0.5, 0.0),
              ncol=3, fontsize=6.5, handletextpad=0.3, columnspacing=0.8)
    subcaption(1, f"(b) {best} - 2D view")

    # (c) preferred layout in 3D
    ax = fig.add_axes([x0[2] / W, (B - 0.55) / H, w3d / W, (h2d + 0.5) / H], projection="3d")
    layout_3d(ax, thr, arrow=0.9, lw=1.5)
    subcaption(2, f"(c) {best} - 3D view")
    save(fig, "preferred_layout.png")


# ------------------------------------------------------------------------------------------ figure 2: all architectures
def fig_all_architectures():
    lays = [r["layout"] for r in ranking]
    ncol = 6
    nrow = math.ceil(len(lays) / ncol)
    cw, ch = 1.05, 1.35                              # cell size [inch]; 6 columns fit a 16 cm A4 text width
    fig = plt.figure(figsize=(ncol * cw, nrow * ch))
    for k, lay in enumerate(lays):
        ax = fig.add_subplot(nrow, ncol, k + 1, projection="3d")
        layout_3d(ax, thrusters(lay), arrow=0.9, lw=1.0)
        ax.set_title(lay, fontsize=8, y=0.995)
    fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
    save(fig, "all_architectures.png")


# ------------------------------------------------------------------------------------------ requirement compliance grid
REQ_IDS = [f"R-{k:02d}" for k in range(1, 17)]        # numbering of the report's numerical requirements table


def requirement_status():
    """status 'met' / 'partly' / 'not' of R-01 ... R-16 for every architecture, in rank order.
    Thrust-dependent requirements are checked with the baseline thruster (F_BASE, T_MIN)."""
    prox = {r["id"]: r["status"] for r in rows("proximity_checks.csv")}

    def prox_status(ids):                            # PASS everywhere = met, any FLAG = partly, any FAIL = not
        st = [prox.get(i, "PASS") for i in ids]
        return "not" if "FAIL" in st else "partly" if "FLAG" in st else "met"

    kos = prox_status(["MIS-04a", "MIS-04c", "RCS-15b", "RCS-15c"])     # KOS transfer + missed-burn cases
    departure = prox_status(["RCS-15a", "RCS-15d"])                     # nominal departure cone + burn budget
    mib = num("mib_tol")                                                # 25 % allocation of the IDSS limits
    aft_ok = num("cant_aft") >= math.degrees(HALF) + 12.0 - 1e-9
    fwd_ok = num("cant_fwd") >= num("cant_fwd_min_plume", 20.0) - 1e-9

    out = {}
    for r in ranking:
        lay, ft = r["layout"], r["fault_tolerance"]
        drv = drivers.get(lay, [])

        def need(prefix):                            # thrust per driver row (relaxed value for relaxed layouts)
            return [float(x["F_relaxed_N"] if ft == "relaxed 2FT" else x["F_N"]) for x in drv
                    if x["driver"].startswith(prefix)]

        def lower(prefix):                           # every case of the driver must accept the baseline thrust
            ok = [F_BASE >= v - 1e-9 for v in need(prefix)]
            return "met" if ok and all(ok) else "partly" if any(ok) else "not"

        rot = min(need("single-jet rotation pulse") or [1e300])         # allowed thrust, 25 % of IDSS rate
        lat = min(need("single lateral pulse") or [1e300])              # allowed thrust, 25 % of IDSS lateral rate
        types = {t["type"] for t in thrusters(lay)}
        plume = (aft_ok or "AFT" not in types) and (fwd_ok or not types & {"FWD", "NOSE"})
        out[lay] = {
            "R-01": lower("detumble"),
            "R-02": "met",                           # Δv budgeted; no time or thrust limit
            "R-03": "met",
            "R-04": "met" if summary.get(lay, {}).get("all_axes_1FT") == "yes" else "partly",
            "R-05": kos,
            "R-06": "partly" if F_BASE <= lat else "not",                   # path on the axis; tracking not shown
            "R-07": "partly" if F_BASE <= min(rot, lat) / mib else "not",   # one pulse within the full IDSS limits;
                                                                            # contact not simulated
            "R-08": "met" if F_BASE <= min(rot, lat) else "not",
            "R-09": departure,
            "R-10": lower("CAM"),
            "R-11": lower("deorbit burn"),
            "R-12": lower("bank reversal"),
            "R-13": {"strict 2FT": "met", "relaxed 2FT": "partly"}.get(ft, "not"),
            "R-14": "met" if float(r["stuck_on_ratio"]) >= STUCK_MIN else "not",
            "R-15": "met" if plume else "not",
            "R-16": "met",                           # margins applied in every calculation
        }
    return out


def fig_requirements_grid():
    status = requirement_status()
    lays = [r["layout"] for r in ranking]
    nr, nc = len(REQ_IDS), len(lays)
    cell, gap = 0.30, 0.06                           # cell size [inch]; gap between cells [fraction of a cell]
    L, R, T, B = 0.55, 0.1, 0.6, 0.55               # margins [inch]
    W, H = L + nc * cell + R, T + nr * cell + B
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([L / W, B / H, nc * cell / W, nr * cell / H])
    ax.set_xlim(0, nc)
    ax.set_ylim(nr, 0)
    ax.axis("off")
    for j, lay in enumerate(lays):
        for i, req in enumerate(REQ_IDS):
            ax.add_patch(Rectangle((j + gap / 2, i + gap / 2), 1 - gap, 1 - gap, fc=STATUS[status[lay][req]][0], lw=0))
        ax.text(j + 0.5, -0.15, lay, ha="center", va="bottom", rotation=90, fontsize=7.5, color=INK)
    for i, req in enumerate(REQ_IDS):
        ax.text(-0.15, i + 0.5, req, ha="right", va="center", fontsize=7.5, color=INK)
    h = [Rectangle((0, 0), 1, 1, fc=STATUS[k][0], lw=0) for k in ("met", "partly", "not")]
    fig.legend(h, [STATUS[k][1] for k in ("met", "partly", "not")], loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, 0.0), fontsize=8, handlelength=1.2, handleheight=1.2)
    save(fig, "requirements_table.png")


if __name__ == "__main__":
    print("project:", ROOT)
    print("preferred layout (rank 1):", ranking[0]["layout"])
    fig_preferred()
    fig_all_architectures()
    fig_requirements_grid()
