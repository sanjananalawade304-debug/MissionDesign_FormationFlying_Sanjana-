# ==============================================================
# TASK 1 — CONSTELLATION REVISIT ANALYSIS for 4 planes and 1 cluster each
# ==============================================================

import orekit_jpype as orekit
orekit.initVM()

from org.orekit.data import DataContext, DirectoryCrawler
from java.io import File
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.time import AbsoluteDate
from org.orekit.time import TimeScalesFactory
from org.orekit.frames import FramesFactory
from org.orekit.frames import TopocentricFrame
from org.orekit.utils import Constants
from org.orekit.utils import IERSConventions
from org.orekit.orbits import KeplerianOrbit
from org.orekit.orbits import PositionAngleType
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.bodies import GeodeticPoint
from org.orekit.orbits import OrbitType

# ==============================================================
# LOAD OREKIT DATA
# ==============================================================

# CHANGE THIS PATH

data_path = "C:/Users/Sanjana/Desktop/MissionDesign/orekit-data-main"
manager = DataContext.getDefault().getDataProvidersManager()
manager.addProvider(DirectoryCrawler(File(data_path)))



# ==============================================================
# CONSTANTS
# ==============================================================

MU = Constants.WGS84_EARTH_MU

EARTH_RADIUS = Constants.WGS84_EARTH_EQUATORIAL_RADIUS

ALTITUDE = 550e3

INCLINATION = np.radians(45.0)

ECCENTRICITY = 0.0

ARG_PERIGEE = 0.0

MIN_ELEVATION = np.radians(10.0)

SIM_DURATION = 7 * 24 * 3600


STEP = 60

SAT_MASS = 10.0


# ==============================================================
# TIME
# ==============================================================

utc = TimeScalesFactory.getUTC()

initial_date = AbsoluteDate(
    2026,
    1,
    1,
    0,
    0,
    0.0,
    utc
)


# ==============================================================
# FRAMES
# ==============================================================

EME2000 = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(
    IERSConventions.IERS_2010,
    True
)


# ==============================================================
# EARTH MODEL
# ==============================================================

earth = OneAxisEllipsoid(
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    Constants.WGS84_EARTH_FLATTENING,
    ITRF
)


# ==============================================================
# CREATE NUMERICAL PROPAGATOR
# ==============================================================

def create_propagator(orbit):

    min_step = 0.1
    max_step = 300.0

    position_tolerance = 10.0

    tolerances = NumericalPropagator.tolerances(
        position_tolerance,
        orbit,
        OrbitType.CARTESIAN
    )

    integrator = DormandPrince853Integrator(
        min_step,
        max_step,
        tolerances[0],
        tolerances[1]
    )

    propagator = NumericalPropagator(
        integrator
    )

    propagator.setOrbitType(
        OrbitType.CARTESIAN
    )

    propagator.setInitialState(
        SpacecraftState(
            orbit,
            SAT_MASS
        )
    )

    # ----------------------------------------------------------
    # J2 PERTURBATION
    # ----------------------------------------------------------

    gravity_provider = GravityFieldFactory.getNormalizedProvider(
        2,
        0
    )

    gravity_model = HolmesFeatherstoneAttractionModel(
        ITRF,
        gravity_provider
    )

    propagator.addForceModel(
        gravity_model
    )

    return propagator


# ==============================================================
# CREATE SATELLITE ORBIT
# ==============================================================

def create_satellite_orbit(
    raan_deg,
    anomaly_deg
):

    sma = EARTH_RADIUS + ALTITUDE

    orbit = KeplerianOrbit(
        sma,
        ECCENTRICITY,
        INCLINATION,
        np.radians(ARG_PERIGEE),
        np.radians(raan_deg),
        np.radians(anomaly_deg),
        PositionAngleType.TRUE,
        EME2000,
        initial_date,
        MU
    )

    return orbit


# ==============================================================
# CREATE CONSTELLATION
# ==============================================================

def create_constellation(
    num_planes,
    sats_per_plane
):

    constellation = []

    # ----------------------------------------------------------
    # UNIFORM RAAN SPACING
    # ----------------------------------------------------------

    raan_spacing = 360 / num_planes

    # ----------------------------------------------------------
    # UNIFORM ANOMALY SPACING
    # ----------------------------------------------------------

    anomaly_spacing = 360 / sats_per_plane

    for p in range(num_planes):

        raan = p * raan_spacing

        for s in range(sats_per_plane):

            anomaly = s * anomaly_spacing

            orbit = create_satellite_orbit(
                raan,
                anomaly
            )

            propagator = create_propagator(
                orbit
            )

            constellation.append(
                propagator
            )

    return constellation


# ==============================================================
# TARGET GRID
# ==============================================================

latitudes = [20, 25, 30, 35, 40]

longitudes = [60, 70, 80, 90, 100]

targets = []

for lat in latitudes:

    for lon in longitudes:

        gp = GeodeticPoint(
            np.radians(lat),
            np.radians(lon),
            0.0
        )

        topo = TopocentricFrame(
            earth,
            gp,
            f"T_{lat}_{lon}"
        )

        targets.append({

            "lat": lat,
            "lon": lon,
            "frame": topo

        })


# ==============================================================
# CHECK ACCESS
# ==============================================================

def has_access(
    propagator,
    current_date,
    topo
):

    state = propagator.propagate(
        current_date
    )

    pv = state.getPVCoordinates(
        ITRF
    )

    sat_position = pv.getPosition()

    elevation = topo.getElevation(
        sat_position,
        ITRF,
        current_date
    )

    return elevation > MIN_ELEVATION


# ==============================================================
# COMPUTE REVISIT
# ==============================================================

def compute_revisit(
    constellation
):

    revisit_results = []

    time_array = np.arange(
        0,
        SIM_DURATION,
        STEP
    )

    for target in tqdm(targets):

        topo = target["frame"]

        access_times = []

        for t in time_array:

            current_date = initial_date.shiftedBy(
                float(t)
            )

            visible = False

            for sat in constellation:

                if has_access(
                    sat,
                    current_date,
                    topo
                ):

                    visible = True
                    break

            if visible:

                access_times.append(
                    t
                )

        # ------------------------------------------------------
        # REVISIT CALCULATION
        # ------------------------------------------------------

        if len(access_times) < 2:

            max_revisit = np.inf

        else:

            access_times = np.array(
                access_times
            )

            gaps = np.diff(
                access_times
            )

            revisit_gaps = gaps[
                gaps > STEP
            ]

            if len(revisit_gaps) == 0:

                max_revisit = 0

            else:

                max_revisit = np.max(
                    revisit_gaps
                ) / 3600.0

        revisit_results.append({

            "lat": target["lat"],
            "lon": target["lon"],
            "max_revisit_hr": max_revisit

        })

    return pd.DataFrame(
        revisit_results
    )

def generate_ground_tracks(
    constellation
):

    tracks = []

    sample_times = np.arange(
        0,
        24*3600,
        300
    )

    for sat_id, propagator in enumerate(constellation):

        lats = []
        lons = []

        for t in sample_times:

            current_date = initial_date.shiftedBy(float(t))

            state = propagator.propagate(current_date)

            pv = state.getPVCoordinates(ITRF)

            pos = pv.getPosition()

            gp = earth.transform(
                pos,
                ITRF,
                current_date
            )

            lats.append(
                np.degrees(gp.getLatitude())
            )

            lons.append(
                np.degrees(gp.getLongitude())
            )

        tracks.append((lats,lons))

    return tracks

def plot_ground_tracks(tracks):

    fig = plt.figure(figsize=(14,8))

    ax = plt.axes(
        projection=ccrs.PlateCarree()
    )

    ax.coastlines()

    ax.add_feature(
        cfeature.BORDERS,
        linewidth=0.5
    )

    ax.gridlines(
        draw_labels=True
    )
    ax.set_extent(
    [40, 120, 0, 60],
    crs=ccrs.PlateCarree()
    )

    # AOI

    rect = Rectangle(
    (60,20),
    40,
    20,
    facecolor='red',
    alpha=0.15,
    edgecolor='red',
    linewidth=2
    )

    ax.add_patch(rect)
    ax.text(
    80,
    30,
    "AOI",
    transform=ccrs.PlateCarree(),
    fontsize=12,
    fontweight='bold',
    color='red',
    ha='center'
    )

    for i,(lats,lons) in enumerate(tracks):

        ax.plot(
            lons,
            lats,
            transform=ccrs.PlateCarree(),
            label=f'Plane {i+1}'
        )

    plt.title(
    'Ground Tracks of Selected Constellation\n'
    '4 Orbital Planes × 1 Cluster per Plane'
    )
    plt.legend()

    plt.savefig(
        'Task1_GroundTracks.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()
#Heatmap plot
def plot_revisit_heatmap(df):

    fig = plt.figure(figsize=(12, 8))

    ax = plt.axes(
        projection=ccrs.PlateCarree()
    )

    ax.coastlines()

    ax.add_feature(
        cfeature.BORDERS,
        linewidth=0.5
    )

    ax.gridlines(
        draw_labels=True
    )

    # ----------------------------------------------------------
    # Zoom into region around AOI
    # ----------------------------------------------------------

    ax.set_extent(
        [40, 120, 0, 60],
        crs=ccrs.PlateCarree()
    )

    # ----------------------------------------------------------
    # Create heatmap data
    # ----------------------------------------------------------

    pivot = df.pivot(
        index='lat',
        columns='lon',
        values='max_revisit_hr'
    )

    lon_grid, lat_grid = np.meshgrid(
        pivot.columns,
        pivot.index
    )

    # ----------------------------------------------------------
    # Heatmap
    # ----------------------------------------------------------

    cf = ax.contourf(
        lon_grid,
        lat_grid,
        pivot.values,
        levels=20,
        cmap='viridis',
        transform=ccrs.PlateCarree()
    )

    # ----------------------------------------------------------
    # Display revisit values on grid points
    # ----------------------------------------------------------

    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):

            ax.text(
                pivot.columns[j],
                pivot.index[i],
                f"{pivot.values[i, j]:.2f}",
                fontsize=8,
                ha='center',
                va='center',
                color='black',
                fontweight='bold',
                transform=ccrs.PlateCarree()
            )

    # ----------------------------------------------------------
    # AOI rectangle
    # ----------------------------------------------------------

    rect = Rectangle(
        (60, 20),
        40,
        20,
        facecolor='red',
        alpha=0.15,
        edgecolor='red',
        linewidth=2
    )

    ax.add_patch(rect)

    # ----------------------------------------------------------
    # AOI label
    # ----------------------------------------------------------

    ax.text(
        80,
        30,
        "AOI",
        fontsize=12,
        fontweight='bold',
        color='red',
        ha='center',
        transform=ccrs.PlateCarree()
    )

    # ----------------------------------------------------------
    # Colorbar
    # ----------------------------------------------------------

    plt.colorbar(
        cf,
        label='Worst Revisit Time [hr]'
    )

    # ----------------------------------------------------------
    # Title
    # ----------------------------------------------------------

    plt.title(
        'Worst Revisit Time Across Area of Interest\n'
        '4 Orbital Planes × 1 Cluster per Plane'
    )

    # ----------------------------------------------------------
    # Save Figure
    # ----------------------------------------------------------

    plt.savefig(
        'Task1_RevisitHeatmap.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()

# ==============================================================
# CONSTELLATION ARCHITECTURE
# ==============================================================

from mpl_toolkits.mplot3d import Axes3D

def plot_constellation_architecture():

    fig = plt.figure(figsize=(10,10))

    ax = fig.add_subplot(
        111,
        projection='3d'
    )

    # ----------------------------------------------------------
    # EARTH
    # ----------------------------------------------------------

    Re = EARTH_RADIUS / 1000

    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = Re * np.outer(
        np.cos(u),
        np.sin(v)
    )

    y = Re * np.outer(
        np.sin(u),
        np.sin(v)
    )

    z = Re * np.outer(
        np.ones(np.size(u)),
        np.cos(v)
    )

    ax.plot_surface(
        x, y, z,
        color='lightblue',
        alpha=0.4,
        linewidth=0
    )

    # ----------------------------------------------------------
    # ORBIT PARAMETERS
    # ----------------------------------------------------------

    r = (EARTH_RADIUS + ALTITUDE)/1000

    nu = np.linspace(
        0,
        2*np.pi,
        500
    )

    RAANS = [0,90,180,270]

    # ----------------------------------------------------------
    # ORBITAL PLANES
    # ----------------------------------------------------------

    for i, raan_deg in enumerate(RAANS):

        raan = np.radians(
            raan_deg
        )

        inc = INCLINATION

        x_orb = r*np.cos(nu)
        y_orb = r*np.sin(nu)
        z_orb = np.zeros_like(nu)

        R3 = np.array([
            [ np.cos(raan),-np.sin(raan),0],
            [ np.sin(raan), np.cos(raan),0],
            [ 0,0,1]
        ])

        R1 = np.array([
            [1,0,0],
            [0,np.cos(inc),-np.sin(inc)],
            [0,np.sin(inc), np.cos(inc)]
        ])

        XYZ = R3 @ R1 @ np.vstack(
            [x_orb,y_orb,z_orb]
        )

        ax.plot(
            XYZ[0],
            XYZ[1],
            XYZ[2],
            linewidth=2,
            label=f'Plane {i+1}'
        )

        # ------------------------------------------------------
        # CLUSTER LOCATION
        # ------------------------------------------------------

        sat = XYZ[:,0]

        ax.scatter(
            sat[0],
            sat[1],
            sat[2],
            s=100,
            marker='o'
        )

        ax.text(
            sat[0],
            sat[1],
            sat[2],
            f'C{i+1}'
        )

    # ----------------------------------------------------------
    # FORMATTING
    # ----------------------------------------------------------

    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    ax.set_zlabel('Z [km]')

    ax.set_title(
        'Selected Constellation Architecture\n'
        '4 Orbital Planes × 1 Cluster per Plane'
    )

    ax.set_box_aspect(
        [1,1,1]
    )

    ax.legend()

    plt.tight_layout()

    plt.savefig(
        'Task1_ConstellationArchitecture.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()
# ==============================================================
# TEST CASES
# ==============================================================

cases = [(4,1)]
summary = []


# ==============================================================
# RUN ANALYSIS
# ==============================================================

for planes, sats_per_plane in cases:

    print("\n================================================")

    print(
        f"RUNNING {planes} PLANES × {sats_per_plane} CLUSTERS"
    )

    print("================================================")

    constellation = create_constellation(
        planes,
        sats_per_plane
    )

    plot_constellation_architecture()

    results = compute_revisit(
        constellation
    )
    tracks = generate_ground_tracks(
    constellation
    )

    plot_ground_tracks(
    tracks
    )

    plot_revisit_heatmap(
    results
    )
    results.to_csv(
    "Task1_RevisitResults.csv",
    index=False
    )

    

    print(results)

    worst_revisit = results[
        "max_revisit_hr"
    ].max()

    total_clusters = planes * sats_per_plane

    total_satellites = total_clusters * 3

    summary.append({

        "Planes": planes,

        "Clusters/Plane": sats_per_plane,

        "Clusters": total_clusters,

        "Total Satellites": total_satellites,

        "Worst Revisit [hr]": worst_revisit

    })


# ==============================================================
# FINAL SUMMARY
# ==============================================================

summary_df = pd.DataFrame(summary)

print("\n================================================")
print("FINAL CONSTELLATION TRADE STUDY")
print("================================================")

print(summary_df)


# ==============================================================
# END
# ==============================================================