# ==============================================================
# TASK 2 — GDOP TRADE STUDY
# Original code unchanged through Block 5.
# Plot section replaced with improved publication-quality figures.
# ==============================================================


# ==============================================================
# BLOCK 1 — MISSION REQUIREMENTS
# ==============================================================

import numpy as np

REQUIRED_ACCURACY = 200.0      # m
CEP_REQUIREMENT = 0.59 * REQUIRED_ACCURACY
TIMING_ERROR_NS   = 100.0      # ns
C                 = 299792458.0  # m/s

sigma_R    = C * TIMING_ERROR_NS * 1e-9
sigma_TDOA = np.sqrt(2) * sigma_R
GDOP_MAX   = REQUIRED_ACCURACY / sigma_TDOA

print("=" * 60)
print("MISSION REQUIREMENTS")
print("=" * 60)
print(f"Required Accuracy      : {REQUIRED_ACCURACY:.1f} m")
print(f"Equivalent CEP Requirement : {CEP_REQUIREMENT:.1f} m")
print(f"Timing Error           : {TIMING_ERROR_NS:.1f} ns")
print("\n")
print(f"Range Error            : {sigma_R:.2f} m")
print(f"TDOA Error             : {sigma_TDOA:.2f} m")
print(f"Maximum Allowed GDOP   : {GDOP_MAX:.2f}")
print("\n")
print("Design Requirement:")
print(f"GDOP must be < {GDOP_MAX:.2f}")


# ==============================================================
# BLOCK 2 — EMITTER + CLUSTER CENTER
# ==============================================================

import orekit_jpype as orekit
orekit.initVM()

from org.orekit.data import DataContext, DirectoryCrawler
from java.io import File

data_path = "C:/Users/Sanjana/Desktop/MissionDesign/orekit-data-main"
manager = DataContext.getDefault().getDataProvidersManager()
manager.addProvider(DirectoryCrawler(File(data_path)))

from org.orekit.frames import FramesFactory
from org.orekit.frames import TopocentricFrame
from org.orekit.utils import Constants
from org.orekit.utils import IERSConventions
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.bodies import GeodeticPoint

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

earth = OneAxisEllipsoid(
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    Constants.WGS84_EARTH_FLATTENING,
    ITRF
)

EMITTER_LAT = 35.0
EMITTER_LON = 75.0
EMITTER_ALT = 0.0

emitter_gp = GeodeticPoint(
    np.radians(EMITTER_LAT),
    np.radians(EMITTER_LON),
    EMITTER_ALT
)

emitter_ecef = earth.transform(emitter_gp)
re = np.array([emitter_ecef.getX(), emitter_ecef.getY(), emitter_ecef.getZ()])

topo = TopocentricFrame(earth, emitter_gp, "Emitter")

CLUSTER_ALTITUDE = 550e3

cluster_gp = GeodeticPoint(
    np.radians(EMITTER_LAT),
    np.radians(EMITTER_LON),
    CLUSTER_ALTITUDE
)

cluster_ecef = earth.transform(cluster_gp)
rc = np.array([cluster_ecef.getX(), cluster_ecef.getY(), cluster_ecef.getZ()])

cluster_range = np.linalg.norm(rc - re)

print("\n")
print("=" * 60)
print("EMITTER")
print("=" * 60)
print(f"Latitude   : {EMITTER_LAT:.3f} deg")
print(f"Longitude  : {EMITTER_LON:.3f} deg")
print(f"Altitude   : {EMITTER_ALT/1000:.3f} km")
print("\nECEF Coordinates [m]")
print(f"X = {re[0]:.3f}")
print(f"Y = {re[1]:.3f}")
print(f"Z = {re[2]:.3f}")

print("\n")
print("=" * 60)
print("CLUSTER CENTER")
print("=" * 60)
print(f"Latitude   : {EMITTER_LAT:.3f} deg")
print(f"Longitude  : {EMITTER_LON:.3f} deg")
print(f"Altitude   : {CLUSTER_ALTITUDE/1000:.3f} km")
print("\nECEF Coordinates [m]")
print(f"X = {rc[0]:.3f}")
print(f"Y = {rc[1]:.3f}")
print(f"Z = {rc[2]:.3f}")

print("\n")
print("=" * 60)
print("GEOMETRY CHECK")
print("=" * 60)
print(f"Emitter → Cluster Center Distance = {cluster_range/1000:.3f} km")


# ==============================================================
# BLOCK 3 — SATELLITE GEOMETRIES AND BASELINE TRADE STUDY
# ==============================================================

import pandas as pd

baseline_values = [
    1e3, 5e3, 10e3, 20e3, 50e3,
    100e3, 150e3, 200e3, 250e3, 300e3
]

results = []

lat0 = np.radians(EMITTER_LAT)
lon0 = np.radians(EMITTER_LON)

east_hat = np.array([-np.sin(lon0),  np.cos(lon0),  0.0])
north_hat = np.array([
    -np.sin(lat0)*np.cos(lon0),
    -np.sin(lat0)*np.sin(lon0),
     np.cos(lat0)
])
up_hat = np.array([
     np.cos(lat0)*np.cos(lon0),
     np.cos(lat0)*np.sin(lon0),
     np.sin(lat0)
])

def generate_geometry(name, s):
    if name == "Linear":
        A = (-s/2,  0.0,  0.0)
        B = ( 0.0,  0.0,  0.0)
        C = ( s/2,  0.0,  0.0)
    elif name == "Equilateral":
        h = np.sqrt(3)/2 * s
        A = ( 0.0,   2*h/3, 0.0)
        B = (-s/2,  -h/3,   0.0)
        C = ( s/2,  -h/3,   0.0)
    elif name == "Isosceles":
        h = 1.2 * s
        A = ( 0.0,   2*h/3, 0.0)
        B = (-s/2,  -h/3,   0.0)
        C = ( s/2,  -h/3,   0.0)
    elif name == "Right Triangle":
        A = (-s/3,  -s/3,  0.0)
        B = ( 2*s/3,-s/3,  0.0)
        C = (-s/3,   2*s/3, 0.0)
    else:
        raise ValueError("Unknown Geometry")
    return {"A": A, "B": B, "C": C}

geometry_list = ["Linear", "Equilateral", "Isosceles", "Right Triangle"]

for geometry in geometry_list:
    for baseline in baseline_values:
        sats    = generate_geometry(geometry, baseline)
        sat_xyz = {}
        for sat_name, (east, north, up) in sats.items():
            r_sat = rc + east*east_hat + north*north_hat + up*up_hat
            sat_xyz[sat_name] = r_sat
        AB = np.linalg.norm(sat_xyz["A"] - sat_xyz["B"]) / 1000
        BC = np.linalg.norm(sat_xyz["B"] - sat_xyz["C"]) / 1000
        AC = np.linalg.norm(sat_xyz["A"] - sat_xyz["C"]) / 1000
        for sat_name in ["A", "B", "C"]:
            xyz = sat_xyz[sat_name]
            gp  = earth.transform(Vector3D(xyz[0], xyz[1], xyz[2]), ITRF, None)
            results.append({
                "Geometry"     : geometry,
                "Baseline_km"  : baseline/1000,
                "Satellite"    : sat_name,
                "Latitude_deg" : np.degrees(gp.getLatitude()),
                "Longitude_deg": np.degrees(gp.getLongitude()),
                "Altitude_km"  : gp.getAltitude()/1000,
                "X_m"          : xyz[0],
                "Y_m"          : xyz[1],
                "Z_m"          : xyz[2],
                "AB_km"        : AB,
                "BC_km"        : BC,
                "AC_km"        : AC
            })

df_states = pd.DataFrame(results)
pd.set_option("display.max_rows",    None)
pd.set_option("display.max_columns", None)

csv_filename = "Geometry_Baseline_Trade_Study1.csv"
df_states.to_csv(csv_filename, index=False)

print("\n")
print("=" * 80)
print("CSV FILE SAVED")
print("=" * 80)
print(f"File Name : {csv_filename}")
print(f"Total Rows: {len(df_states)}")


# ==============================================================
# BLOCK 4 — 2D TDOA GDOP ANALYSIS
# ==============================================================

gdop_results = []

grouped = df_states.groupby(
    ["Geometry", "Baseline_km"]
)

for (geometry, baseline), group in grouped:

    satA = group[group["Satellite"] == "A"].iloc[0]
    satB = group[group["Satellite"] == "B"].iloc[0]
    satC = group[group["Satellite"] == "C"].iloc[0]

    rA = np.array([
        satA["X_m"],
        satA["Y_m"],
        satA["Z_m"]
    ])

    rB = np.array([
        satB["X_m"],
        satB["Y_m"],
        satB["Z_m"]
    ])

    rC = np.array([
        satC["X_m"],
        satC["Y_m"],
        satC["Z_m"]
    ])

    # ---------------------------------------------------------
    # LOS VECTORS
    # ---------------------------------------------------------

    uA = (rA - re) / np.linalg.norm(rA - re)
    uB = (rB - re) / np.linalg.norm(rB - re)
    uC = (rC - re) / np.linalg.norm(rC - re)

    # ---------------------------------------------------------
    # EN PROJECTIONS
    # ---------------------------------------------------------

    uA_E = np.dot(uA, east_hat)
    uA_N = np.dot(uA, north_hat)

    uB_E = np.dot(uB, east_hat)
    uB_N = np.dot(uB, north_hat)

    uC_E = np.dot(uC, east_hat)
    uC_N = np.dot(uC, north_hat)

    # ---------------------------------------------------------
    # 2D TDOA GEOMETRY MATRIX
    # ---------------------------------------------------------

    H = np.array([
        [
            uB_E - uA_E,
            uB_N - uA_N
        ],
        [
            uC_E - uA_E,
            uC_N - uA_N
        ]
    ])

    # ---------------------------------------------------------
    # CHECK SINGULARITY
    # ---------------------------------------------------------

    detH = np.linalg.det(H)

    if abs(detH) < 1e-12:

        gdop = np.inf
        position_error = np.inf
        CEP = np.inf

    else:

        P = sigma_TDOA**2 * np.linalg.inv(H.T @ H)

        gdop = np.sqrt(
            np.trace(P)
        ) / sigma_TDOA

        position_error = np.sqrt(
            np.trace(P)
        )
        # ---------------------------------------------------------
        # CEP CALCULATION
        # ---------------------------------------------------------

        CEP = 0.59 * position_error
    meets_requirement = (
        position_error <= REQUIRED_ACCURACY
    )

    gdop_results.append({

    "Geometry":
        geometry,

    "Baseline_km":
        baseline,

    "GDOP":
        gdop,

    "Position_Error_m":
        position_error,

    "CEP_m":
        CEP,

    "Meets_Requirement":
        meets_requirement

})
df_gdop = pd.DataFrame(
    gdop_results
).sort_values(
    ["Geometry", "Baseline_km"]
)

df_gdop.to_csv(
    "GDOP_Trade_Study.csv",
    index=False
)

# ==============================================================
# BLOCK 5 — FINAL TRADE STUDY SUMMARY
# ==============================================================

feasible = df_gdop[
    df_gdop["Meets_Requirement"] == True
]

summary_results = []

for geometry in df_gdop["Geometry"].unique():

    geom = feasible[
        feasible["Geometry"] == geometry
    ]

    if len(geom) == 0:

        summary_results.append({
            "Geometry"        : geometry,
            "Baseline_km"     : np.nan,
            "GDOP"            : np.nan,
            "Position_Error_m": np.nan,
            "CEP_m"           : np.nan
        })

    else:

        # BEST RESULT FOR THIS GEOMETRY
        best = geom.loc[
            geom["Position_Error_m"].idxmin()
        ]

        summary_results.append({
            "Geometry"        : best["Geometry"],
            "Baseline_km"     : best["Baseline_km"],
            "GDOP"            : best["GDOP"],
            "Position_Error_m": best["Position_Error_m"],
            "CEP_m"           : best["CEP_m"]
        })

df_summary = pd.DataFrame(summary_results)

print("\n")
print("=" * 80)
print("SUMMARY TABLE")
print("=" * 80)
print(df_summary)

df_summary.to_csv(
    "Task2_Final_Summary.csv",
    index=False
)

print("\nTask2_Final_Summary.csv saved")
# ==============================================================
# AUTOMATIC SELECTION OF OPTIMAL CONFIGURATION
# ==============================================================

selected_config = feasible.loc[
    feasible["Position_Error_m"].idxmin()
]

OPTIMAL_GEOMETRY = selected_config["Geometry"]
OPTIMAL_BASELINE = selected_config["Baseline_km"]

selected_cluster = df_states[
    (df_states["Geometry"] == OPTIMAL_GEOMETRY) &
    (np.isclose(df_states["Baseline_km"], OPTIMAL_BASELINE))
]


# Best overall
best_overall  = feasible.loc[feasible["Position_Error_m"].idxmin()]
best_geometry = best_overall["Geometry"]
best_baseline = best_overall["Baseline_km"]

best_cluster = df_states[
    (df_states["Geometry"] == best_geometry) &
    (np.isclose(df_states["Baseline_km"], best_baseline))
]

print("\n")
print("=" * 80)
print("BEST OVERALL CONFIGURATION")
print("=" * 80)
print(f"Geometry        : {best_geometry}")
print(f"Baseline        : {best_baseline:.1f} km")
print(f"GDOP            : {best_overall['GDOP']:.3f}")
print(f"Position Error  : {best_overall['Position_Error_m']:.2f} m")
print(f"CEP             : {best_overall['CEP_m']:.2f} m")
print()
print("Satellite States")
print("-" * 80)
print(best_cluster[
    ["Satellite","Latitude_deg","Longitude_deg","Altitude_km","X_m","Y_m","Z_m"]
].to_string(index=False))


# ==============================================================
# IMPROVED PLOTS
# ==============================================================

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---------------------------------------------------------------
# GLOBAL STYLE
# ---------------------------------------------------------------

plt.rcParams.update({
    'font.family'      : 'sans-serif',
    'font.size'        : 11,
    'axes.titlesize'   : 13,
    'axes.labelsize'   : 12,
    'legend.fontsize'  : 10,
    'xtick.labelsize'  : 10,
    'ytick.labelsize'  : 10,
    'axes.grid'        : True,
    'grid.alpha'       : 0.35,
    'grid.linestyle'   : '--',
    'figure.facecolor' : 'white',
    'axes.facecolor'   : '#F8F8F8',
    'axes.spines.top'  : False,
    'axes.spines.right': False,
})

# ---------------------------------------------------------------
# GEOMETRY STYLE MAP  (consistent across both plots)
# ---------------------------------------------------------------

geo_style = {
    'Linear': {
        'color': '#CC3333', 'marker': 'X',
        'linestyle': '--',  'linewidth': 1.8,
        'markersize': 7,    'zorder': 2
    },
    'Equilateral': {
        'color': '#2266BB', 'marker': 's',
        'linestyle': '-',   'linewidth': 2.0,
        'markersize': 7,    'zorder': 3
    },
    'Isosceles': {
        'color': '#228855', 'marker': '^',
        'linestyle': '-',   'linewidth': 2.0,
        'markersize': 7,    'zorder': 3
    },
    'Right Triangle': {
        'color': '#CC7700', 'marker': 'o',
        'linestyle': '-',   'linewidth': 2.5,
        'markersize': 8,    'zorder': 4
    },
}

ORDER = ['Linear', 'Equilateral', 'Isosceles', 'Right Triangle']


# ==============================================================
# PLOT 1 — GDOP vs BASELINE (improved)
# ==============================================================

fig1, ax1 = plt.subplots(figsize=(11, 7))

# Feasible zone shading
ax1.axhspan(0, GDOP_MAX, alpha=0.09, color='#22BB55', zorder=0)

# Requirement line
ax1.axhline(
    y=GDOP_MAX,
    color='black', linestyle='--', linewidth=2.0, zorder=5,
    label=f'Requirement  (GDOP = {GDOP_MAX:.2f})'
)

# Plot all geometries
for geometry in ORDER:
    geom_df = df_gdop[df_gdop['Geometry'] == geometry].sort_values('Baseline_km')
    s = geo_style[geometry]
    ax1.plot(
        geom_df['Baseline_km'], geom_df['GDOP'],
        color=s['color'],  marker=s['marker'],
        linestyle=s['linestyle'], linewidth=s['linewidth'],
        markersize=s['markersize'], zorder=s['zorder'],
        label=geometry, clip_on=True
    )

# ── Optimal point star ──────────────────────────────────────────
optimal_gdop = df_gdop[
    (df_gdop['Geometry'] == OPTIMAL_GEOMETRY) &
    (np.isclose(df_gdop['Baseline_km'], OPTIMAL_BASELINE))
]

if len(optimal_gdop):
    gv = optimal_gdop['GDOP'].values[0]
    ax1.plot(OPTIMAL_BASELINE, gv, '*',
             color='#CC7700', markersize=18, zorder=6,
             markeredgecolor='#333333', markeredgewidth=0.6)
    ax1.annotate(
        f'Optimal\n{OPTIMAL_GEOMETRY} @ {OPTIMAL_BASELINE:.0f} km\nGDOP = {gv:.2f}',
        xy=(200, gv), xytext=(210, gv + 6),
        fontsize=9, color='#884400', fontweight='bold',
        arrowprops=dict(
            arrowstyle='->', color='#CC7700', lw=1.5,
            connectionstyle='arc3,rad=0.15'
        ),
        bbox=dict(boxstyle='round,pad=0.35', facecolor='#FFF8E7',
                  edgecolor='#CC7700', alpha=0.92)
    )

# ── Feasible zone text ──────────────────────────────────────────
ax1.text(2, 0.5, 'Feasible zone  (GDOP < 4.72)',
         fontsize=9, color='#116633', fontstyle='italic', alpha=0.85)

# ── Clipping note ───────────────────────────────────────────────
ax1.text(2, 47.8,
         '✱  Linear geometry: GDOP > 10⁶ at small baselines (extends off chart)',
         fontsize=8.5, color='#CC3333', fontstyle='italic')

# ── Axes ────────────────────────────────────────────────────────
ax1.set_xlim(0, 315)
ax1.set_ylim(0, 50)
ax1.set_xlabel('Baseline (km)', fontsize=12)
ax1.set_ylabel('GDOP', fontsize=12)
ax1.set_title(
    'GDOP vs Baseline — Satellite Cluster Geometry Trade Study\n'
    'Emitter: 35°N, 75°E  |  Cluster altitude: 550 km  |  σ_timing = 100 ns',
    fontsize=13, pad=14
)

# ── Legend (add feasible patch) ─────────────────────────────────
feasible_patch = mpatches.Patch(
    facecolor='#22BB55', alpha=0.22, label='Feasible zone'
)
h, l = ax1.get_legend_handles_labels()
ax1.legend(handles=h + [feasible_patch], labels=l + ['Feasible zone'],
           loc='upper right', framealpha=0.92, edgecolor='#CCCCCC',
           borderpad=0.8, labelspacing=0.5)

plt.tight_layout()
plt.savefig('GDOP_vs_Baseline.png', dpi=300, bbox_inches='tight')
plt.show()


# ==============================================================
# PLOT 2 — POSITION ERROR vs BASELINE
# ==============================================================

fig2, ax2 = plt.subplots(figsize=(11, 7))

# --------------------------------------------------------------
# Feasible zone
# --------------------------------------------------------------

ax2.axhspan(
    0,
    REQUIRED_ACCURACY,
    alpha=0.09,
    color='#22BB55',
    zorder=0
)

# --------------------------------------------------------------
# Requirement line
# --------------------------------------------------------------

ax2.axhline(
    y=REQUIRED_ACCURACY,
    color='black',
    linestyle='--',
    linewidth=2.0,
    zorder=5,
    label=f'Requirement ({REQUIRED_ACCURACY:.0f} m)'
)

# --------------------------------------------------------------
# Plot all geometries
# --------------------------------------------------------------

for geometry in ORDER:

    geom_df = (
        df_gdop[df_gdop['Geometry'] == geometry]
        .sort_values('Baseline_km')
    )

    s = geo_style[geometry]

    ax2.plot(
        geom_df['Baseline_km'],
        geom_df['Position_Error_m'],
        color=s['color'],
        marker=s['marker'],
        linestyle=s['linestyle'],
        linewidth=s['linewidth'],
        markersize=s['markersize'],
        zorder=s['zorder'],
        label=geometry
    )

# --------------------------------------------------------------
# Optimal configuration
# --------------------------------------------------------------

optimal_error = df_gdop[
    (df_gdop['Geometry'] == OPTIMAL_GEOMETRY) &
    (np.isclose(df_gdop['Baseline_km'], OPTIMAL_BASELINE))
]

if len(optimal_error):

    ev = optimal_error['Position_Error_m'].values[0]

    ax2.plot(
        OPTIMAL_BASELINE,
        ev,
        marker='*',
        color='#CC7700',
        markersize=20,
        zorder=10,
        markeredgecolor='black',
        markeredgewidth=0.8
    )

    ax2.annotate(
        f'Optimal Configuration\n'
        f'{OPTIMAL_GEOMETRY}\n'
        f'Baseline = {OPTIMAL_BASELINE:.0f} km\n'
        f'Error = {ev:.1f} m',
        xy=(OPTIMAL_BASELINE, ev),
        xytext=(OPTIMAL_BASELINE - 90, ev + 55),
        fontsize=9,
        color='#884400',
        fontweight='bold',
        arrowprops=dict(
            arrowstyle='->',
            color='#CC7700',
            lw=1.5
        ),
        bbox=dict(
            boxstyle='round,pad=0.4',
            facecolor='#FFF8E7',
            edgecolor='#CC7700',
            alpha=0.95
        )
    )

    # ----------------------------------------------------------
    # Margin annotation
    # ----------------------------------------------------------

    margin_pct = (
        (REQUIRED_ACCURACY - ev)
        / REQUIRED_ACCURACY
        * 100
    )

    x_margin = min(OPTIMAL_BASELINE + 15, 300)

    ax2.annotate(
        '',
        xy=(x_margin, ev),
        xytext=(x_margin, REQUIRED_ACCURACY),
        arrowprops=dict(
            arrowstyle='<->',
            color='#777777',
            lw=1.2
        )
    )

    ax2.text(
        x_margin + 5,
        (ev + REQUIRED_ACCURACY) / 2,
        f'{margin_pct:.0f}%\nmargin',
        fontsize=8.5,
        color='#555555',
        va='center'
    )

# --------------------------------------------------------------
# Notes
# --------------------------------------------------------------

ax2.text(
    2,
    8,
    f'Feasible zone (Error < {REQUIRED_ACCURACY:.0f} m)',
    fontsize=9,
    color='#116633',
    fontstyle='italic'
)

ax2.text(
    2,
    477,
    'Linear geometry exhibits extremely large errors at small baselines',
    fontsize=8.5,
    color='#CC3333',
    fontstyle='italic'
)

# --------------------------------------------------------------
# Axes
# --------------------------------------------------------------

ax2.set_xlim(0, 315)
ax2.set_ylim(0, 500)

ax2.set_xlabel(
    'Baseline (km)',
    fontsize=12
)

ax2.set_ylabel(
    'Position Error (m)',
    fontsize=12
)

ax2.set_title(
    'Position Error vs Baseline — Satellite Cluster Geometry Trade Study\n'
    'Emitter: 35°N, 75°E  |  Cluster Altitude: 550 km  |  σtiming = 100 ns',
    fontsize=13,
    pad=14
)

ax2.grid(
    True,
    linestyle='--',
    alpha=0.4
)

# --------------------------------------------------------------
# Legend
# --------------------------------------------------------------

feasible_patch = mpatches.Patch(
    facecolor='#22BB55',
    alpha=0.22,
    label='Feasible Zone'
)

handles, labels = ax2.get_legend_handles_labels()

ax2.legend(
    handles + [feasible_patch],
    labels + ['Feasible Zone'],
    loc='upper right',
    framealpha=0.92
)

# --------------------------------------------------------------
# Save
# --------------------------------------------------------------

plt.tight_layout()

plt.savefig(
    'Position_Error_vs_Baseline.png',
    dpi=300,
    bbox_inches='tight'
)

plt.show()
# ==============================================================
# PLOT 3 — CEP vs BASELINE
# ==============================================================

fig3, ax3 = plt.subplots(figsize=(11, 7))

# ---------------------------------------------------------------
# Mission Requirement
# ---------------------------------------------------------------

ax3.axhspan(
    0,
    CEP_REQUIREMENT,
    alpha=0.09,
    color='#22BB55',
    zorder=0
)

ax3.axhline(
    y=CEP_REQUIREMENT,
    color='black',
    linestyle='--',
    linewidth=2.0,
    label='Requirement (200 m)'
)

# ---------------------------------------------------------------
# Plot all geometries
# ---------------------------------------------------------------

for geometry in ORDER:

    geom_df = df_gdop[
        df_gdop['Geometry'] == geometry
    ].sort_values('Baseline_km')

    s = geo_style[geometry]

    ax3.plot(
        geom_df['Baseline_km'],
        geom_df['CEP_m'],
        color=s['color'],
        marker=s['marker'],
        linestyle=s['linestyle'],
        linewidth=s['linewidth'],
        markersize=s['markersize'],
        zorder=s['zorder'],
        label=geometry
    )



optimal_cep = df_gdop[
    (df_gdop['Geometry'] == OPTIMAL_GEOMETRY) &
    (np.isclose(df_gdop['Baseline_km'], OPTIMAL_BASELINE))
]

if len(optimal_cep):

    cep_val = optimal_cep['CEP_m'].values[0]

    ax3.plot(
        OPTIMAL_BASELINE,
        cep_val,
        '*',
        color='#CC7700',
        markersize=18,
        markeredgecolor='#333333',
        markeredgewidth=0.6,
        zorder=6
    )

    ax3.annotate(
        f'{OPTIMAL_GEOMETRY} @ {OPTIMAL_BASELINE:.0f} km',
        xy=(200, cep_val),
        xytext=(210, cep_val + 20),
        fontsize=9,
        color='#884400',
        fontweight='bold',
        arrowprops=dict(
            arrowstyle='->',
            color='#CC7700',
            lw=1.5
        ),
        bbox=dict(
            boxstyle='round,pad=0.35',
            facecolor='#FFF8E7',
            edgecolor='#CC7700',
            alpha=0.92
        )
    )

# ---------------------------------------------------------------
# Axes
# ---------------------------------------------------------------

ax3.set_xlim(0, 315)
ax3.set_ylim(0, 250)

ax3.set_xlabel('Baseline (km)')
ax3.set_ylabel('CEP (m)')

ax3.set_title(
    'CEP vs Baseline — Satellite Cluster Geometry Trade Study\n'
    'Emitter: 35°N, 75°E | Cluster altitude: 550 km | σ_timing = 100 ns'
)

# ---------------------------------------------------------------
# Legend
# ---------------------------------------------------------------

feasible_patch3 = mpatches.Patch(
    facecolor='#22BB55',
    alpha=0.22,
    label='Feasible zone'
)

h3, l3 = ax3.get_legend_handles_labels()

ax3.legend(
    handles=h3 + [feasible_patch3],
    labels=l3 + ['Feasible zone'],
    loc='upper right',
    framealpha=0.92,
    edgecolor='#CCCCCC'
)

plt.tight_layout()

plt.savefig(
    'CEP_vs_Baseline.png',
    dpi=300,
    bbox_inches='tight'
)

plt.tight_layout()
plt.savefig('CEP_vs_Baseline.png', dpi=300, bbox_inches='tight')
plt.show()
# ==============================================================
# END
# ==============================================================



