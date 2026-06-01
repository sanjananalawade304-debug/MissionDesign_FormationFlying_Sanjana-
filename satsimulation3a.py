# ============================================================
# BLOCK 1
# OREKIT SETUP + CHIEF ORBIT
# ============================================================
#Analysis 2 with drag
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import orekit_jpype as orekit
orekit.initVM()

from java.io import File

from org.orekit.data import (
    DataContext,
    DirectoryCrawler
)

# ============================================================
# OREKIT DATA
# ============================================================

data_path = r"C:/Users/Sanjana/Desktop/MissionDesign/orekit-data-main"

manager = (
    DataContext
    .getDefault()
    .getDataProvidersManager()
)

manager.addProvider(
    DirectoryCrawler(
        File(data_path)
    )
)

# ============================================================
# IMPORTS
# ============================================================

from org.orekit.frames import FramesFactory
from org.orekit.bodies import OneAxisEllipsoid

from org.orekit.forces.drag import (DragForce,IsotropicDrag)

from org.orekit.models.earth.atmosphere import (NRLMSISE00)
from org.orekit.bodies import CelestialBodyFactory

from org.orekit.models.earth.atmosphere.data import (MarshallSolarActivityFutureEstimation)

from org.orekit.models.earth.atmosphere.data.MarshallSolarActivityFutureEstimation import (StrengthLevel)
from org.orekit.time import (AbsoluteDate,TimeScalesFactory)

from org.orekit.utils import (Constants,PVCoordinates)

from org.orekit.orbits import (KeplerianOrbit,CartesianOrbit,PositionAngleType)
from org.hipparchus.geometry.euclidean.threed import (Vector3D)

# ============================================================
# CONSTANTS
# ============================================================

MU = Constants.WGS84_EARTH_MU

R_E = Constants.WGS84_EARTH_EQUATORIAL_RADIUS

ALTITUDE = 550e3

a = R_E + ALTITUDE

INC = np.radians(45.0)

# ============================================================
# FRAME
# ============================================================

EME2000 = FramesFactory.getEME2000()

# ============================================================
# EPOCH
# ============================================================

utc = TimeScalesFactory.getUTC()

epoch = AbsoluteDate(
    2026,
    1,
    1,
    0,
    0,
    0.0,
    utc
)

# ============================================================
# CHIEF ORBIT (SATELLITE B)
# ============================================================

orbit_B = KeplerianOrbit(
    a,
    0.0,
    INC,
    0.0,
    0.0,
    0.0,
    PositionAngleType.MEAN,
    EME2000,
    epoch,
    MU
)

# ============================================================
# CHIEF STATE
# ============================================================

pvB = orbit_B.getPVCoordinates()

rB = np.array([
    pvB.getPosition().getX(),
    pvB.getPosition().getY(),
    pvB.getPosition().getZ()
])

vB = np.array([
    pvB.getVelocity().getX(),
    pvB.getVelocity().getY(),
    pvB.getVelocity().getZ()
])

# ============================================================
# RTN FRAME
# ============================================================

R_hat = rB / np.linalg.norm(rB)

H_vec = np.cross(
    rB,
    vB
)

N_hat = (
    H_vec /
    np.linalg.norm(H_vec)
)

T_hat = np.cross(
    N_hat,
    R_hat
)

# ============================================================
# CHIEF ORBIT INFO
# ============================================================

n = np.sqrt(
    MU / a**3
)

T_orbit = (
    2*np.pi / n
)

print("="*60)
print("TASK 3 - ANALYSIS 1")
print("="*60)

print(f"\nAltitude       : {ALTITUDE/1000:.1f} km")
print(f"Inclination    : {np.degrees(INC):.1f} deg")

print(
    f"Mean Motion    : "
    f"{n:.6e} rad/s"
)

print(
    f"Orbit Period   : "
    f"{T_orbit/60:.2f} min"
)

print("\nBlock 1 Complete")
# ============================================================
# BLOCK 2
# ROE -> DEPUTY ORBITS
# ============================================================

baseline = 300e3

# ============================================================
# ALONG TRACK DEPUTY (SAT A)
# ============================================================

delta_lambda = -baseline / a

delta_M = delta_lambda

# ============================================================
# CROSS TRACK DEPUTY (SAT C)
# ============================================================

delta_i_y = -baseline / a

delta_RAAN = (
    delta_i_y /
    np.sin(INC)
)

# ============================================================
# SATELLITE A
# ============================================================

orbit_A = KeplerianOrbit(
    a,
    0.0,
    INC,
    0.0,
    0.0,
    delta_M,
    PositionAngleType.MEAN,
    EME2000,
    epoch,
    MU
)

# ============================================================
# SATELLITE C
# ============================================================

orbit_C = KeplerianOrbit(
    a,
    0.0,
    INC,
    0.0,
    delta_RAAN,
    0.0,
    PositionAngleType.MEAN,
    EME2000,
    epoch,
    MU
)

# ============================================================
# INITIAL STATES
# ============================================================

pvA = orbit_A.getPVCoordinates()

pvB = orbit_B.getPVCoordinates()

pvC = orbit_C.getPVCoordinates()

rA = np.array([
    pvA.getPosition().getX(),
    pvA.getPosition().getY(),
    pvA.getPosition().getZ()
])

rB = np.array([
    pvB.getPosition().getX(),
    pvB.getPosition().getY(),
    pvB.getPosition().getZ()
])

rC = np.array([
    pvC.getPosition().getX(),
    pvC.getPosition().getY(),
    pvC.getPosition().getZ()
])

# ============================================================
# CHECK SMA
# ============================================================

print("\n")
print("="*60)
print("INITIAL ORBIT ELEMENTS")
print("="*60)

print(
    f"A SMA = {orbit_A.getA()/1000:.6f} km"
)

print(
    f"B SMA = {orbit_B.getA()/1000:.6f} km"
)

print(
    f"C SMA = {orbit_C.getA()/1000:.6f} km"
)

# ============================================================
# INITIAL GEOMETRY
# ============================================================

AB = np.linalg.norm(
    rA-rB
)/1000

BC = np.linalg.norm(
    rC-rB
)/1000

AC = np.linalg.norm(
    rA-rC
)/1000

print("\n")
print("="*60)
print("INITIAL FORMATION")
print("="*60)

print(f"AB = {AB:.3f} km")
print(f"BC = {BC:.3f} km")
print(f"AC = {AC:.3f} km")

print("\nBlock 2 Complete")
print("\nStarting 24 hr J2 propagation...")
from tqdm import tqdm

from org.orekit.propagation import (
    SpacecraftState
)

from org.orekit.propagation.numerical import (
    NumericalPropagator
)

from org.hipparchus.ode.nonstiff import (
    DormandPrince853Integrator
)

from org.orekit.forces.gravity import (
    HolmesFeatherstoneAttractionModel
)

from org.orekit.forces.gravity.potential import (
    GravityFieldFactory
)

from org.orekit.frames import FramesFactory

from org.orekit.utils import (
    IERSConventions
)

# ============================================================
# ITRF FRAME
# ============================================================

ITRF = FramesFactory.getITRF(
    IERSConventions.IERS_2010,
    True
)
# ============================================================
# EARTH MODEL
# ============================================================

earth = OneAxisEllipsoid(
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    Constants.WGS84_EARTH_FLATTENING,
    ITRF
)

# ============================================================
# J2 GRAVITY MODEL
# ============================================================

gravity_provider = (
    GravityFieldFactory
    .getNormalizedProvider(
        2,
        0
    )
)

gravity_force = (
    HolmesFeatherstoneAttractionModel(
        ITRF,
        gravity_provider
    )
)


# ============================================================
# SOLAR ACTIVITY MODEL
# ============================================================

msafe = MarshallSolarActivityFutureEstimation(
    MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES,
    StrengthLevel.AVERAGE
)

# ============================================================
# ATMOSPHERE
# ============================================================

sun = CelestialBodyFactory.getSun()

atmosphere = NRLMSISE00(
    msafe,
    sun,
    earth
)

# ============================================================
# PROPAGATOR FACTORY
# ============================================================
# ============================================================
# BALLISTIC COEFFICIENTS
# ============================================================

mass_A = 100.0
mass_B = 100.0
mass_C = 100.0

Cd = 2.2

BC_A = 50.0
BC_B = 45.0
BC_C = 45.0

area_A = mass_A/(Cd*BC_A)
area_B = mass_B/(Cd*BC_B)
area_C = mass_C/(Cd*BC_C)

print("\n")
print("="*60)
print("DRAG MODEL")
print("="*60)

print(f"A Area = {area_A:.4f} m²")
print(f"B Area = {area_B:.4f} m²")
print(f"C Area = {area_C:.4f} m²")

def create_propagator(
    orbit,
    area,
    mass
):

    integrator = DormandPrince853Integrator(
        0.1,
        300.0,
        1e-10,
        1e-10
    )

    propagator = NumericalPropagator(
        integrator
    )

    propagator.setInitialState(
        SpacecraftState(
            orbit,
            mass
        )
    )

    propagator.addForceModel(
        gravity_force
    )

    drag_sensitive = IsotropicDrag(
        area,
        Cd
    )

    drag_force = DragForce(
        atmosphere,
        drag_sensitive
    )

    propagator.addForceModel(
        drag_force
    )

    return propagator
# ============================================================
# CREATE PROPAGATORS
# ============================================================

prop_A = create_propagator(
    orbit_A,
    area_A,
    mass_A
)

prop_B = create_propagator(
    orbit_B,
    area_B,
    mass_B
)

prop_C = create_propagator(
    orbit_C,
    area_C,
    mass_C
)

# ============================================================
# ECI -> RTN
# ============================================================

def eci_to_rtn(
    r_ref,
    v_ref,
    r_rel
):

    R_hat = (
        r_ref /
        np.linalg.norm(r_ref)
    )

    h = np.cross(
        r_ref,
        v_ref
    )

    N_hat = (
        h /
        np.linalg.norm(h)
    )

    T_hat = np.cross(
        N_hat,
        R_hat
    )

    transform = np.column_stack(
        (
            R_hat,
            T_hat,
            N_hat
        )
    )

    return (
        transform.T @ r_rel
    )

# ============================================================
# TIME
# ============================================================

SIM_TIME = 3*24*3600

STEP = 60

times = np.arange(
    0,
    SIM_TIME + STEP,
    STEP
)

# ============================================================
# STORAGE
# ============================================================

A_R = []
A_T = []
A_N = []

C_R = []
C_T = []
C_N = []

AB_sep = []
BC_sep = []
AC_sep = []

A_SMA = []
B_SMA = []
C_SMA = []

print("\nStarting propagation...\n")

# ============================================================
# PROPAGATION LOOP
# ============================================================

for t in tqdm(times):

    current_date = epoch.shiftedBy(
        float(t)
    )

    state_A = prop_A.propagate(
        current_date
    )

    state_B = prop_B.propagate(
        current_date
    )

    state_C = prop_C.propagate(
        current_date
    )

    pvA = state_A.getPVCoordinates(
        EME2000
    )

    pvB = state_B.getPVCoordinates(
        EME2000
    )

    pvC = state_C.getPVCoordinates(
        EME2000
    )

    rA = np.array([
        pvA.getPosition().getX(),
        pvA.getPosition().getY(),
        pvA.getPosition().getZ()
    ])

    rB = np.array([
        pvB.getPosition().getX(),
        pvB.getPosition().getY(),
        pvB.getPosition().getZ()
    ])

    rC = np.array([
        pvC.getPosition().getX(),
        pvC.getPosition().getY(),
        pvC.getPosition().getZ()
    ])

    vB = np.array([
        pvB.getVelocity().getX(),
        pvB.getVelocity().getY(),
        pvB.getVelocity().getZ()
    ])

    rAB = rA - rB
    rCB = rC - rB

    rtn_AB = eci_to_rtn(
        rB,
        vB,
        rAB
    )

    rtn_CB = eci_to_rtn(
        rB,
        vB,
        rCB
    )

    A_R.append(
        rtn_AB[0]/1000
    )

    A_T.append(
        rtn_AB[1]/1000
    )

    A_N.append(
        rtn_AB[2]/1000
    )

    C_R.append(
        rtn_CB[0]/1000
    )

    C_T.append(
        rtn_CB[1]/1000
    )

    C_N.append(
        rtn_CB[2]/1000
    )

    AB_sep.append(
        np.linalg.norm(
            rAB
        )/1000
    )

    BC_sep.append(
        np.linalg.norm(
            rCB
        )/1000
    )

    AC_sep.append(
        np.linalg.norm(
            rA-rC
        )/1000
    )

    A_SMA.append(
        state_A.getOrbit().getA()/1000
    )

    B_SMA.append(
        state_B.getOrbit().getA()/1000
    )

    C_SMA.append(
        state_C.getOrbit().getA()/1000
    )

print(AB_sep[0], AB_sep[-1])

print(BC_sep[0], BC_sep[-1])

print(AC_sep[0], AC_sep[-1])

print("\n")
print(A_SMA[0], A_SMA[-1])

print(B_SMA[0], B_SMA[-1])

print(C_SMA[0], C_SMA[-1])

print("\nPropagation Complete")

print("\nBlock 3 Complete")
# ============================================================
# SAVE CSV
# ============================================================

results = pd.DataFrame({

    "Time_hr" : times/3600,

    # Satellite A RTN
    "A_R_km" : A_R,
    "A_T_km" : A_T,
    "A_N_km" : A_N,

    # Satellite C RTN
    "C_R_km" : C_R,
    "C_T_km" : C_T,
    "C_N_km" : C_N,

    # Formation Geometry
    "AB_sep_km" : AB_sep,
    "BC_sep_km" : BC_sep,
    "AC_sep_km" : AC_sep,

    # SMA History
    "A_SMA_km" : A_SMA,
    "B_SMA_km" : B_SMA,
    "C_SMA_km" : C_SMA

})

results.to_csv(
    "Task3_Analysis2_J2_Drag.csv",
    index=False
)

print("\nCSV Saved")

print(
    f"Rows = {len(results)}"
)

print(
    "File = Task3_Analysis2_J2_Drag.csv"
)
print("\nFINAL SEPARATIONS")

print("AB =", AB_sep[-1])
print("BC =", BC_sep[-1])
print("AC =", AC_sep[-1])

print("\nFINAL SMA")

print("A =", A_SMA[-1])
print("B =", B_SMA[-1])
print("C =", C_SMA[-1])

# ============================================================
# PLOT 1
# RELATIVE POSITION VS TIME
# ============================================================

time_hr = times / 3600

plt.figure(figsize=(10,6))

plt.plot(
    time_hr,
    A_T,
    label='Satellite A Along-Track'
)

plt.plot(
    time_hr,
    C_N,
    label='Satellite C Cross-Track'
)

plt.xlabel('Time [hr]')
plt.ylabel('Relative Position [km]')

plt.title(
    'Relative Position Evolution'
)

plt.grid(True)
plt.legend()

plt.tight_layout()

plt.savefig(
    'Task3_RelativePosition_2.png',
    dpi=300
)

plt.show()

# ============================================================
# PLOT 2
# FORMATION SEPARATION
# ============================================================

plt.figure(figsize=(10,6))

plt.plot(
    time_hr,
    AB_sep,
    label='AB'
)

plt.plot(
    time_hr,
    BC_sep,
    label='BC'
)

plt.plot(
    time_hr,
    AC_sep,
    label='AC'
)

plt.xlabel('Time [hr]')
plt.ylabel('Separation [km]')

plt.title(
    'Formation Separation Evolution'
)

plt.grid(True)
plt.legend()

plt.tight_layout()

plt.savefig(
    'Task3_FormationSeparation_2.png',
    dpi=300
)

plt.show()

# ============================================================
# PLOT 3
# SEMI MAJOR AXIS
# ============================================================

plt.figure(figsize=(10,6))

plt.plot(
    time_hr,
    A_SMA,
    label='Satellite A'
)

plt.plot(
    time_hr,
    B_SMA,
    label='Satellite B'
)

plt.plot(
    time_hr,
    C_SMA,
    label='Satellite C'
)

plt.xlabel('Time [hr]')
plt.ylabel('Semi-Major Axis [km]')

plt.title(
    'Semi-Major Axis Evolution'
)

plt.grid(True)
plt.legend()

plt.tight_layout()

plt.savefig(
    'Task3_SMA_2.png',
    dpi=300
)

plt.show()

# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n")
print("="*60)
print("ANALYSIS 2 COMPLETE")
print("="*60)

print("\nFinal Separations")

print(f"AB = {AB_sep[-1]:.3f} km")
print(f"BC = {BC_sep[-1]:.3f} km")
print(f"AC = {AC_sep[-1]:.3f} km")

print("\nFinal SMA")

print(f"A = {A_SMA[-1]:.3f} km")
print(f"B = {B_SMA[-1]:.3f} km")
print(f"C = {C_SMA[-1]:.3f} km")

print("\nFiles Generated")

print("Task3_Analysis1_J2_Drag.csv")
print("Task3_RelativePosition_2.png")
print("Task3_FormationSeparation_2.png")
print("Task3_SMA_2.png")