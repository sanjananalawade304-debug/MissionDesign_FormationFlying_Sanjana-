# ==============================================================
# TASK 1 — CONSTELLATION REVISIT ANALYSIS
# ==============================================================

import orekit_jpype as orekit
orekit.initVM()

from org.orekit.data import DataContext, DirectoryCrawler
from java.io import File
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
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


# ==============================================================
# TEST CASES
# ==============================================================

cases = [

    (3,1),
    (4,1),
    (5,1),
    (6,1),
    (7,1),
    (8,1)

  
]
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

    results = compute_revisit(
        constellation
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
# PLOT RESULTS
# ==============================================================

plt.figure(figsize=(10,6))

plt.plot(

    summary_df["Clusters"],
    summary_df["Worst Revisit [hr]"],
    marker='o',
    linewidth=2

)

plt.axhline(

    3,
    color='red',
    linestyle='--',
    label='3 hr Requirement'

)

plt.xlabel("Number of Clusters")

plt.ylabel("Worst Revisit [hr]")

plt.title(
    "Constellation Revisit Trade Study"
)

plt.grid(True)

plt.legend()
plt.savefig( 'Constellation Revisit Trade Study.png', dpi=300)
plt.show()


# ==============================================================
# END
# ==============================================================