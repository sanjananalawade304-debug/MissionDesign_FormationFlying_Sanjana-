# ==============================================================
# TASK 2 — GDOP TRADE STUDY WITH OFF-NADIR CLUSTER OFFSET
# ==============================================================
# CORRECTIONS APPLIED:
#   1. Cluster center offset from emitter (off-nadir geometry)
#   2. ENU basis recomputed at cluster center (not emitter)
#   3. Elevation angle check added
#   4. Offset trade study block added
#   5. Plot y-axis range widened for readability
#   6. Best configuration output at end
# ==============================================================


# ==============================================================
# BLOCK 1 — MISSION REQUIREMENTS
# ==============================================================

import numpy as np

REQUIRED_ACCURACY = 200.0        # m
TIMING_ERROR_NS   = 10.0        # ns
C                 = 299792458.0  # m/s

sigma_R    = C * TIMING_ERROR_NS * 1e-9
sigma_TDOA = np.sqrt(2) * sigma_R
GDOP_MAX   = REQUIRED_ACCURACY / sigma_TDOA

print("=" * 60)
print("MISSION REQUIREMENTS")
print("=" * 60)
print(f"Required Accuracy      : {REQUIRED_ACCURACY:.1f} m")
print(f"Timing Error           : {TIMING_ERROR_NS:.1f} ns")
print(f"Range Error            : {sigma_R:.2f} m")
print(f"TDOA Error             : {sigma_TDOA:.2f} m")
print(f"Maximum Allowed GDOP   : {GDOP_MAX:.2f}")
print(f"\nDesign Requirement: GDOP must be < {GDOP_MAX:.2f}")


# ==============================================================
# BLOCK 2 — EMITTER + CLUSTER CENTER (OFF-NADIR)
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

ITRF = FramesFactory.getITRF(
    IERSConventions.IERS_2010,
    True
)

earth = OneAxisEllipsoid(
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    Constants.WGS84_EARTH_FLATTENING,
    ITRF
)

# ----------------------------------------------------------
# EMITTER LOCATION
# ----------------------------------------------------------

EMITTER_LAT = 35.0
EMITTER_LON = 75.0
EMITTER_ALT = 0.0

emitter_gp = GeodeticPoint(
    np.radians(EMITTER_LAT),
    np.radians(EMITTER_LON),
    EMITTER_ALT
)

emitter_ecef = earth.transform(emitter_gp)
re = np.array([
    emitter_ecef.getX(),
    emitter_ecef.getY(),
    emitter_ecef.getZ()
])

topo = TopocentricFrame(earth, emitter_gp, "Emitter")

# ----------------------------------------------------------
# UP VECTOR AT EMITTER  (used for elevation angle check)
# ----------------------------------------------------------

lat_e = np.radians(EMITTER_LAT)
lon_e = np.radians(EMITTER_LON)

up_hat_emitter = np.array([
    np.cos(lat_e) * np.cos(lon_e),
    np.cos(lat_e) * np.sin(lon_e),
    np.sin(lat_e)
])

# ----------------------------------------------------------
# CLUSTER OFFSET PARAMETERS  ← KEY CORRECTION
# ----------------------------------------------------------
#
# Set the cluster center OFFSET from the emitter in lat/lon.
# A non-zero offset means satellites are no longer directly
# overhead, which increases angular diversity and lowers GDOP.
#
# Rule of thumb:
#   3° lat ≈ 333 km offset  →  elevation angle ≈ 59°
#   5° lat ≈ 555 km offset  →  elevation angle ≈ 45° (near-optimal)
# ----------------------------------------------------------

CLUSTER_ALTITUDE        = 550e3   # m
CLUSTER_LAT_OFFSET_DEG  = 0.0    # degrees north offset (~333 km)
CLUSTER_LON_OFFSET_DEG  = 0.0    # degrees east  offset

CLUSTER_LAT = EMITTER_LAT + CLUSTER_LAT_OFFSET_DEG
CLUSTER_LON = EMITTER_LON + CLUSTER_LON_OFFSET_DEG

cluster_gp = GeodeticPoint(
    np.radians(CLUSTER_LAT),
    np.radians(CLUSTER_LON),
    CLUSTER_ALTITUDE
)

cluster_ecef = earth.transform(cluster_gp)
rc = np.array([
    cluster_ecef.getX(),
    cluster_ecef.getY(),
    cluster_ecef.getZ()
])

# ----------------------------------------------------------
# GEOMETRY CHECK — slant range + elevation angle
# ----------------------------------------------------------

los_vec             = rc - re
los_range           = np.linalg.norm(los_vec)
los_unit            = los_vec / los_range
elevation_angle_deg = np.degrees(
    np.arcsin(np.dot(los_unit, up_hat_emitter))
)

print("\n")
print("=" * 60)
print("EMITTER")
print("=" * 60)
print(f"Latitude   : {EMITTER_LAT:.3f} deg")
print(f"Longitude  : {EMITTER_LON:.3f} deg")
print(f"Altitude   : {EMITTER_ALT / 1000:.3f} km")
print(f"\nECEF Coordinates [m]")
print(f"X = {re[0]:.3f}")
print(f"Y = {re[1]:.3f}")
print(f"Z = {re[2]:.3f}")

print("\n")
print("=" * 60)
print("CLUSTER CENTER  (OFF-NADIR)")
print("=" * 60)
print(f"Latitude   : {CLUSTER_LAT:.3f} deg  (offset: +{CLUSTER_LAT_OFFSET_DEG:.1f}°)")
print(f"Longitude  : {CLUSTER_LON:.3f} deg  (offset: +{CLUSTER_LON_OFFSET_DEG:.1f}°)")
print(f"Altitude   : {CLUSTER_ALTITUDE / 1000:.3f} km")
print(f"\nECEF Coordinates [m]")
print(f"X = {rc[0]:.3f}")
print(f"Y = {rc[1]:.3f}")
print(f"Z = {rc[2]:.3f}")

print("\n")
print("=" * 60)
print("GEOMETRY CHECK")
print("=" * 60)
print(f"Emitter → Cluster Slant Range : {los_range / 1000:.3f} km")
print(f"Elevation Angle (from emitter) : {elevation_angle_deg:.2f}°")
print(f"  (90° = nadir / directly overhead — worst GDOP)")
print(f"  (45° = oblique — near-optimal GDOP)")


# ==============================================================
# BLOCK 3 — SATELLITE GEOMETRIES AND BASELINE TRADE STUDY
# ==============================================================

import pandas as pd

baseline_values = [
    1e3,
    5e3,
    10e3,
    20e3,
    50e3,
    100e3,
    150e3,
    200e3,
    250e3,
    300e3,
    400e3,
    500e3
]

results = []

# ----------------------------------------------------------
# ENU BASIS AT CLUSTER CENTER  ← CORRECTED
# Previously used EMITTER lat/lon — now uses CLUSTER lat/lon
# ----------------------------------------------------------

lat0 = np.radians(CLUSTER_LAT)   # cluster center latitude
lon0 = np.radians(CLUSTER_LON)   # cluster center longitude

east_hat = np.array([
    -np.sin(lon0),
     np.cos(lon0),
     0.0
])

north_hat = np.array([
    -np.sin(lat0) * np.cos(lon0),
    -np.sin(lat0) * np.sin(lon0),
     np.cos(lat0)
])

up_hat = np.array([
     np.cos(lat0) * np.cos(lon0),
     np.cos(lat0) * np.sin(lon0),
     np.sin(lat0)
])

# ----------------------------------------------------------
# GEOMETRY GENERATOR
# ----------------------------------------------------------

def generate_geometry(name, s):

    if name == "Linear":
        A = (-s / 2,  0.0,  0.0)
        B = ( 0.0,    0.0,  0.0)
        C = ( s / 2,  0.0,  0.0)

    elif name == "Equilateral":
        h = np.sqrt(3) / 2 * s
        A = ( 0.0,   2 * h / 3, 0.0)
        B = (-s / 2, -h / 3,    0.0)
        C = ( s / 2, -h / 3,    0.0)

    elif name == "Isosceles":
        h = 1.2 * s
        A = ( 0.0,   2 * h / 3, 0.0)
        B = (-s / 2, -h / 3,    0.0)
        C = ( s / 2, -h / 3,    0.0)

    elif name == "Right Triangle":
        A = (-s / 3,  -s / 3,  0.0)
        B = ( 2*s / 3,-s / 3,  0.0)
        C = (-s / 3,   2*s / 3, 0.0)

    else:
        raise ValueError("Unknown Geometry")

    return {"A": A, "B": B, "C": C}

# ----------------------------------------------------------
# GEOMETRIES TO TEST
# ----------------------------------------------------------

geometry_list = [
    "Linear",
    "Equilateral",
    "Isosceles",
    "Right Triangle"
]

# ----------------------------------------------------------
# MAIN LOOP
# ----------------------------------------------------------

for geometry in geometry_list:

    for baseline in baseline_values:

        sats = generate_geometry(geometry, baseline)
        sat_xyz = {}

        # ENU -> ECEF  (using corrected cluster-center ENU basis)
        for sat_name, (east, north, up) in sats.items():
            r_sat = (
                rc
                + east  * east_hat
                + north * north_hat
                + up    * up_hat
            )
            sat_xyz[sat_name] = r_sat

        # Side lengths
        AB = np.linalg.norm(sat_xyz["A"] - sat_xyz["B"]) / 1000
        BC = np.linalg.norm(sat_xyz["B"] - sat_xyz["C"]) / 1000
        AC = np.linalg.norm(sat_xyz["A"] - sat_xyz["C"]) / 1000

        for sat_name in ["A", "B", "C"]:
            xyz = sat_xyz[sat_name]
            gp = earth.transform(
                Vector3D(xyz[0], xyz[1], xyz[2]),
                ITRF,
                None
            )
            results.append({
                "Geometry"     : geometry,
                "Baseline_km"  : baseline / 1000,
                "Satellite"    : sat_name,
                "Latitude_deg" : np.degrees(gp.getLatitude()),
                "Longitude_deg": np.degrees(gp.getLongitude()),
                "Altitude_km"  : gp.getAltitude() / 1000,
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
# BLOCK 4 — GDOP COMPUTATION
# ==============================================================

gdop_results = []
grouped = df_states.groupby(["Geometry", "Baseline_km"])

for (geometry, baseline), group in grouped:

    satA = group[group["Satellite"] == "A"].iloc[0]
    satB = group[group["Satellite"] == "B"].iloc[0]
    satC = group[group["Satellite"] == "C"].iloc[0]

    rA = np.array([satA["X_m"], satA["Y_m"], satA["Z_m"]])
    rB = np.array([satB["X_m"], satB["Y_m"], satB["Z_m"]])
    rC = np.array([satC["X_m"], satC["Y_m"], satC["Z_m"]])

    # LOS unit vectors
    uA = (rA - re) / np.linalg.norm(rA - re)
    uB = (rB - re) / np.linalg.norm(rB - re)
    uC = (rC - re) / np.linalg.norm(rC - re)

    # TDOA geometry matrix
    H = np.array([
        uB - uA,
        uC - uA
    ])

    G     = H.T @ H
    G_inv = np.linalg.pinv(G)
    gdop  = np.sqrt(np.trace(G_inv))

    position_error    = gdop * sigma_TDOA
    meets_requirement = (position_error <= REQUIRED_ACCURACY)

    gdop_results.append({
        "Geometry"          : geometry,
        "Baseline_km"       : baseline,
        "GDOP"              : gdop,
        "Position_Error_m"  : position_error,
        "Meets_Requirement" : meets_requirement
    })

df_gdop = pd.DataFrame(gdop_results).sort_values(
    ["Geometry", "Baseline_km"]
)

csv_name = "GDOP_Trade_Study.csv"
df_gdop.to_csv(csv_name, index=False)

print("\n")
print("=" * 80)
print("GDOP TRADE STUDY CSV SAVED")
print("=" * 80)
print(csv_name)


# ==============================================================
# BLOCK 5 — FINAL TRADE STUDY SUMMARY
# ==============================================================

feasible = df_gdop[df_gdop["Meets_Requirement"] == True]
summary_results = []

for geometry in df_gdop["Geometry"].unique():

    geom = feasible[feasible["Geometry"] == geometry]

    if len(geom) == 0:
        summary_results.append({
            "Geometry"         : geometry,
            "Baseline_km"      : np.nan,
            "GDOP"             : np.nan,
            "Position_Error_m" : np.nan
        })
    else:
        best = geom.sort_values("Baseline_km").iloc[0]
        summary_results.append({
            "Geometry"         : geometry,
            "Baseline_km"      : best["Baseline_km"],
            "GDOP"             : best["GDOP"],
            "Position_Error_m" : best["Position_Error_m"]
        })

df_summary = pd.DataFrame(summary_results)

print("\n")
print("=" * 80)
print("SUMMARY TABLE")
print("=" * 80)
print(df_summary)

df_summary.to_csv("Task2_Final_Summary.csv", index=False)
print("\nTask2_Final_Summary.csv saved")

# ----------------------------------------------------------
# SELECTED CONFIGURATION — minimum feasible Right Triangle
# ----------------------------------------------------------

selected_geometry = "Right Triangle"

rt_row = df_summary[df_summary["Geometry"] == selected_geometry]

if rt_row["Baseline_km"].isna().values[0]:
    print("\nRight Triangle does not meet requirement at tested baselines.")
else:
    selected_baseline = rt_row["Baseline_km"].values[0]

    selected_cluster = df_states[
        (df_states["Geometry"] == selected_geometry) &
        (np.isclose(df_states["Baseline_km"], selected_baseline))
    ]

    print("\n")
    print("=" * 80)
    print("SELECTED CONFIGURATION SATELLITE STATES")
    print("=" * 80)
    print(f"Geometry : {selected_geometry}")
    print(f"Baseline : {selected_baseline:.1f} km")
    print()
    print(
        selected_cluster[
            ["Satellite", "Latitude_deg", "Longitude_deg",
             "Altitude_km", "X_m", "Y_m", "Z_m"]
        ].to_string(index=False)
    )

# ----------------------------------------------------------
# BEST OVERALL CONFIGURATION
# ----------------------------------------------------------

if len(feasible) > 0:
    best_overall  = feasible.loc[feasible["Position_Error_m"].idxmin()]
    best_geometry = best_overall["Geometry"]
    best_baseline = best_overall["Baseline_km"]

    best_cluster = df_states[
        (df_states["Geometry"] == best_geometry) &
        (np.isclose(df_states["Baseline_km"], best_baseline))
    ]

    print("\n")
    print("=" * 80)
    print("BEST OVERALL CONFIGURATION  (minimum position error)")
    print("=" * 80)
    print(f"Geometry        : {best_geometry}")
    print(f"Baseline        : {best_baseline:.1f} km")
    print(f"GDOP            : {best_overall['GDOP']:.3f}")
    print(f"Position Error  : {best_overall['Position_Error_m']:.2f} m")
    print()
    print("Satellite States")
    print("-" * 80)
    print(
        best_cluster[
            ["Satellite", "Latitude_deg", "Longitude_deg",
             "Altitude_km", "X_m", "Y_m", "Z_m"]
        ].to_string(index=False)
    )


# ==============================================================
# BLOCK 6 — PLOTS
# Same style as original — y-axis range widened for readability
# ==============================================================

import matplotlib.pyplot as plt

# ----------------------------------------------------------
# PLOT 1 — GDOP vs Baseline
# ----------------------------------------------------------

plt.figure(figsize=(10, 6))

for geometry in df_gdop["Geometry"].unique():
    geom_df = df_gdop[df_gdop["Geometry"] == geometry]
    plt.plot(
        geom_df["Baseline_km"],
        geom_df["GDOP"],
        marker='o',
        linewidth=2,
        label=geometry
    )

plt.axhline(
    y=GDOP_MAX,
    color='black',
    linestyle='--',
    linewidth=2,
    label=f"Requirement (GDOP={GDOP_MAX:.2f})"
)

plt.xlabel("Baseline (km)", fontsize=12)
plt.ylabel("GDOP", fontsize=12)
plt.title(
    "GDOP vs Baseline for Different Satellite Geometries",
    fontsize=14
)
plt.ylim(0, 50)       # widened — shows all curves; clips Linear extremes
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("GDOP_vs_Baseline.png", dpi=300, bbox_inches="tight")
plt.show()

# ----------------------------------------------------------
# PLOT 2 — Position Error vs Baseline
# ----------------------------------------------------------

plt.figure(figsize=(10, 6))

for geometry in df_gdop["Geometry"].unique():
    geom_df = df_gdop[df_gdop["Geometry"] == geometry]
    plt.plot(
        geom_df["Baseline_km"],
        geom_df["Position_Error_m"],
        marker='o',
        linewidth=2,
        label=geometry
    )

plt.axhline(
    y=REQUIRED_ACCURACY,
    color='black',
    linestyle='--',
    linewidth=2,
    label="200 m Requirement"
)

plt.xlabel("Baseline (km)", fontsize=12)
plt.ylabel("Position Error (m)", fontsize=12)
plt.title("Position Error vs Baseline", fontsize=14)
plt.ylim(0, 2000)     # widened — all geometries visible
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("Position_Error_vs_Baseline.png", dpi=300, bbox_inches="tight")
plt.show()


# ==============================================================
# BLOCK 7 — OFFSET TRADE STUDY
# Sweeps cluster offset from 0° to 5° to find the minimum
# baseline needed at each offset (Right Triangle geometry only)
# ==============================================================

print("\n")
print("=" * 70)
print("OFFSET TRADE STUDY")
print("Right Triangle geometry — sweeping cluster offset 0° to 5°")
print("=" * 70)

offset_values_deg = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

offset_baseline_values = [
    1e3,
    5e3,
    10e3,
    20e3,
    50e3,
    100e3,
    150e3,
    200e3,
    250e3,
    300e3,
    400e3,
    500e3
]

offset_results = []

for offset_deg in offset_values_deg:

    # Recompute cluster center for this offset
    c_lat = np.radians(EMITTER_LAT + offset_deg)
    c_lon = np.radians(EMITTER_LON + CLUSTER_LON_OFFSET_DEG)

    c_gp   = GeodeticPoint(c_lat, c_lon, CLUSTER_ALTITUDE)
    c_ecef = earth.transform(c_gp)
    rc_off = np.array([
        c_ecef.getX(),
        c_ecef.getY(),
        c_ecef.getZ()
    ])

    # ENU at this cluster center
    e_hat_off = np.array([
        -np.sin(c_lon),
         np.cos(c_lon),
         0.0
    ])

    n_hat_off = np.array([
        -np.sin(c_lat) * np.cos(c_lon),
        -np.sin(c_lat) * np.sin(c_lon),
         np.cos(c_lat)
    ])

    # Elevation angle from emitter to this cluster
    los_off  = rc_off - re
    los_u    = los_off / np.linalg.norm(los_off)
    elev_off = np.degrees(
        np.arcsin(np.dot(los_u, up_hat_emitter))
    )

    for baseline in offset_baseline_values:

        s = baseline

        # Right Triangle geometry
        sats_enu = {
            "A": np.array([-s / 3,    -s / 3,   0.0]),
            "B": np.array([ 2*s / 3,  -s / 3,   0.0]),
            "C": np.array([-s / 3,     2*s / 3,  0.0])
        }

        sat_xyz_off = {}
        for name, enu in sats_enu.items():
            east, north, _ = enu
            sat_xyz_off[name] = (
                rc_off
                + east  * e_hat_off
                + north * n_hat_off
            )

        rA = sat_xyz_off["A"]
        rB = sat_xyz_off["B"]
        rC = sat_xyz_off["C"]

        uA = (rA - re) / np.linalg.norm(rA - re)
        uB = (rB - re) / np.linalg.norm(rB - re)
        uC = (rC - re) / np.linalg.norm(rC - re)

        H     = np.array([uB - uA, uC - uA])
        G     = H.T @ H
        G_inv = np.linalg.pinv(G)
        gdop  = np.sqrt(np.trace(G_inv))
        pos_err = gdop * sigma_TDOA

        offset_results.append({
            "Offset_deg"       : offset_deg,
            "Offset_km"        : offset_deg * 111.0,
            "Elevation_deg"    : elev_off,
            "Baseline_km"      : baseline / 1000,
            "GDOP"             : gdop,
            "Position_Error_m" : pos_err,
            "Meets_Req"        : pos_err <= REQUIRED_ACCURACY
        })

df_offset = pd.DataFrame(offset_results)

# ----------------------------------------------------------
# OFFSET SUMMARY TABLE
# ----------------------------------------------------------

print()
print(f"{'Offset':>8}  {'Offset km':>10}  {'Elev':>6}  "
      f"{'Min Baseline':>13}  {'GDOP':>7}  {'Pos Err (m)':>12}")
print("-" * 70)

offset_summary_rows = []

for offset_deg in offset_values_deg:

    sub  = df_offset[
        (df_offset["Offset_deg"] == offset_deg) &
        (df_offset["Meets_Req"]  == True)
    ]
    elev = df_offset[
        df_offset["Offset_deg"] == offset_deg
    ]["Elevation_deg"].iloc[0]

    if len(sub) == 0:
        print(f"{offset_deg:>7.0f}°  {offset_deg*111:>9.0f} km"
              f"  {elev:>5.1f}°  {'N/A':>13}  {'N/A':>7}  {'N/A':>12}")
        offset_summary_rows.append({
            "Offset_deg"       : offset_deg,
            "Offset_km"        : offset_deg * 111.0,
            "Elevation_deg"    : elev,
            "Min_Baseline_km"  : np.nan,
            "GDOP"             : np.nan,
            "Position_Error_m" : np.nan
        })
    else:
        best = sub.sort_values("Baseline_km").iloc[0]
        print(f"{offset_deg:>7.0f}°  {offset_deg*111:>9.0f} km"
              f"  {elev:>5.1f}°"
              f"  {best['Baseline_km']:>10.0f} km"
              f"  {best['GDOP']:>7.3f}"
              f"  {best['Position_Error_m']:>10.1f} m")
        offset_summary_rows.append({
            "Offset_deg"       : offset_deg,
            "Offset_km"        : offset_deg * 111.0,
            "Elevation_deg"    : elev,
            "Min_Baseline_km"  : best["Baseline_km"],
            "GDOP"             : best["GDOP"],
            "Position_Error_m" : best["Position_Error_m"]
        })

df_offset_summary = pd.DataFrame(offset_summary_rows)
df_offset_summary.to_csv("Offset_Trade_Summary.csv", index=False)
print("\nOffset_Trade_Summary.csv saved")

# ----------------------------------------------------------
# OFFSET TRADE PLOTS
# ----------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for offset_deg in offset_values_deg:
    sub = df_offset[df_offset["Offset_deg"] == offset_deg]
    lbl = f"{offset_deg:.0f}° offset ({offset_deg * 111:.0f} km)"

    axes[0].plot(
        sub["Baseline_km"],
        sub["GDOP"],
        marker='o',
        linewidth=2,
        label=lbl
    )
    axes[1].plot(
        sub["Baseline_km"],
        sub["Position_Error_m"],
        marker='o',
        linewidth=2,
        label=lbl
    )

axes[0].axhline(
    GDOP_MAX,
    color='black',
    linestyle='--',
    linewidth=2,
    label=f"Requirement (GDOP={GDOP_MAX:.2f})"
)
axes[0].set_xlabel("Baseline (km)", fontsize=12)
axes[0].set_ylabel("GDOP", fontsize=12)
axes[0].set_title("GDOP vs Baseline — Offset Trade Study", fontsize=14)
axes[0].set_ylim(0, 50)
axes[0].grid(True)
axes[0].legend(fontsize=9)

axes[1].axhline(
    REQUIRED_ACCURACY,
    color='black',
    linestyle='--',
    linewidth=2,
    label="200 m Requirement"
)
axes[1].set_xlabel("Baseline (km)", fontsize=12)
axes[1].set_ylabel("Position Error (m)", fontsize=12)
axes[1].set_title("Position Error vs Baseline — Offset Trade Study", fontsize=14)
axes[1].set_ylim(0, 2000)
axes[1].grid(True)
axes[1].legend(fontsize=9)

plt.tight_layout()
plt.savefig("Offset_Trade_Study.png", dpi=300, bbox_inches="tight")
plt.show()


# ==============================================================
# BLOCK 8 — BEST CONFIGURATION
# ==============================================================

print("\n")
print("=" * 80)
print("BEST CONFIGURATION ANALYSIS")
print("=" * 80)

feasible_offset = df_offset[df_offset["Meets_Req"] == True]

if len(feasible_offset) == 0:
    print("No configuration meets the 200 m requirement in the tested range.")

else:

    # Option A — smallest baseline (operationally easiest)
    best_ops = feasible_offset.sort_values(
        ["Baseline_km", "GDOP"]
    ).iloc[0]

    # Option B — best GDOP (most accurate)
    best_acc = feasible_offset.loc[
        feasible_offset["GDOP"].idxmin()
    ]

    print()
    print("OPTION A — OPERATIONALLY OPTIMAL  (smallest baseline)")
    print("-" * 60)
    print(f"  Cluster Offset   : {best_ops['Offset_deg']:.0f}°  "
          f"({best_ops['Offset_km']:.0f} km north of emitter)")
    print(f"  Elevation Angle  : {best_ops['Elevation_deg']:.2f}°")
    print(f"  Baseline         : {best_ops['Baseline_km']:.0f} km")
    print(f"  Geometry         : Right Triangle")
    print(f"  GDOP             : {best_ops['GDOP']:.3f}")
    print(f"  Position Error   : {best_ops['Position_Error_m']:.2f} m")
    print(f"  Meets 200 m Req  : YES")

    print()
    print("OPTION B — ACCURACY OPTIMAL  (minimum GDOP)")
    print("-" * 60)
    print(f"  Cluster Offset   : {best_acc['Offset_deg']:.0f}°  "
          f"({best_acc['Offset_km']:.0f} km north of emitter)")
    print(f"  Elevation Angle  : {best_acc['Elevation_deg']:.2f}°")
    print(f"  Baseline         : {best_acc['Baseline_km']:.0f} km")
    print(f"  Geometry         : Right Triangle")
    print(f"  GDOP             : {best_acc['GDOP']:.3f}")
    print(f"  Position Error   : {best_acc['Position_Error_m']:.2f} m")
    print(f"  Meets 200 m Req  : YES")

    # ----------------------------------------------------------
    # SATELLITE STATES FOR OPERATIONALLY OPTIMAL CONFIGURATION
    # ----------------------------------------------------------

    opt_offset_deg = best_ops["Offset_deg"]
    opt_baseline   = best_ops["Baseline_km"] * 1000

    c_lat_opt = np.radians(EMITTER_LAT + opt_offset_deg)
    c_lon_opt = np.radians(EMITTER_LON + CLUSTER_LON_OFFSET_DEG)

    c_gp_opt   = GeodeticPoint(c_lat_opt, c_lon_opt, CLUSTER_ALTITUDE)
    c_ecef_opt = earth.transform(c_gp_opt)
    rc_opt     = np.array([
        c_ecef_opt.getX(),
        c_ecef_opt.getY(),
        c_ecef_opt.getZ()
    ])

    e_hat_opt = np.array([
        -np.sin(c_lon_opt),
         np.cos(c_lon_opt),
         0.0
    ])
    n_hat_opt = np.array([
        -np.sin(c_lat_opt) * np.cos(c_lon_opt),
        -np.sin(c_lat_opt) * np.sin(c_lon_opt),
         np.cos(c_lat_opt)
    ])

    s = opt_baseline
    sats_opt = {
        "A": np.array([-s / 3,    -s / 3,   0.0]),
        "B": np.array([ 2*s / 3,  -s / 3,   0.0]),
        "C": np.array([-s / 3,     2*s / 3,  0.0])
    }

    print()
    print("SATELLITE STATES — OPERATIONALLY OPTIMAL CONFIGURATION")
    print("-" * 80)
    print(f"{'Sat':<5} {'Lat (deg)':>12} {'Lon (deg)':>12} "
          f"{'Alt (km)':>10} {'X (m)':>16} {'Y (m)':>16} {'Z (m)':>16}")
    print("-" * 80)

    for name, enu in sats_opt.items():
        east, north, _ = enu
        xyz = rc_opt + east * e_hat_opt + north * n_hat_opt
        gp  = earth.transform(
            Vector3D(xyz[0], xyz[1], xyz[2]),
            ITRF,
            None
        )
        print(
            f"{name:<5}"
            f"{np.degrees(gp.getLatitude()):>12.4f}"
            f"{np.degrees(gp.getLongitude()):>12.4f}"
            f"{gp.getAltitude() / 1000:>10.3f}"
            f"{xyz[0]:>16.3f}"
            f"{xyz[1]:>16.3f}"
            f"{xyz[2]:>16.3f}"
        )

    # ----------------------------------------------------------
    # SATELLITE STATES FOR ACCURACY OPTIMAL CONFIGURATION
    # ----------------------------------------------------------

    acc_offset_deg = best_acc["Offset_deg"]
    acc_baseline   = best_acc["Baseline_km"] * 1000

    c_lat_acc = np.radians(EMITTER_LAT + acc_offset_deg)
    c_lon_acc = np.radians(EMITTER_LON + CLUSTER_LON_OFFSET_DEG)

    c_gp_acc   = GeodeticPoint(c_lat_acc, c_lon_acc, CLUSTER_ALTITUDE)
    c_ecef_acc = earth.transform(c_gp_acc)
    rc_acc     = np.array([
        c_ecef_acc.getX(),
        c_ecef_acc.getY(),
        c_ecef_acc.getZ()
    ])

    e_hat_acc = np.array([
        -np.sin(c_lon_acc),
         np.cos(c_lon_acc),
         0.0
    ])
    n_hat_acc = np.array([
        -np.sin(c_lat_acc) * np.cos(c_lon_acc),
        -np.sin(c_lat_acc) * np.sin(c_lon_acc),
         np.cos(c_lat_acc)
    ])

    s = acc_baseline
    sats_acc = {
        "A": np.array([-s / 3,    -s / 3,   0.0]),
        "B": np.array([ 2*s / 3,  -s / 3,   0.0]),
        "C": np.array([-s / 3,     2*s / 3,  0.0])
    }

    print()
    print("SATELLITE STATES — ACCURACY OPTIMAL CONFIGURATION")
    print("-" * 80)
    print(f"{'Sat':<5} {'Lat (deg)':>12} {'Lon (deg)':>12} "
          f"{'Alt (km)':>10} {'X (m)':>16} {'Y (m)':>16} {'Z (m)':>16}")
    print("-" * 80)

    for name, enu in sats_acc.items():
        east, north, _ = enu
        xyz = rc_acc + east * e_hat_acc + north * n_hat_acc
        gp  = earth.transform(
            Vector3D(xyz[0], xyz[1], xyz[2]),
            ITRF,
            None
        )
        print(
            f"{name:<5}"
            f"{np.degrees(gp.getLatitude()):>12.4f}"
            f"{np.degrees(gp.getLongitude()):>12.4f}"
            f"{gp.getAltitude() / 1000:>10.3f}"
            f"{xyz[0]:>16.3f}"
            f"{xyz[1]:>16.3f}"
            f"{xyz[2]:>16.3f}"
        )

print("\n")
print("=" * 80)
print("END OF TASK 2")
print("=" * 80)

# ==============================================================
# END
# ==============================================================