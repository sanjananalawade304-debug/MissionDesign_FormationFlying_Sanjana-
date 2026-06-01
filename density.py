# ============================================================
# NRLMSISE-00 ATMOSPHERIC DENSITY CALCULATOR
# ============================================================

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

import numpy as np

from org.orekit.frames import FramesFactory
from org.orekit.utils import Constants, IERSConventions
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint

from org.orekit.time import AbsoluteDate
from org.orekit.time import TimeScalesFactory

from org.orekit.models.earth.atmosphere import NRLMSISE00

from org.orekit.models.earth.atmosphere.data import (
    MarshallSolarActivityFutureEstimation
)

from org.orekit.models.earth.atmosphere.data.MarshallSolarActivityFutureEstimation import (
    StrengthLevel
)

# ============================================================
# EARTH MODEL
# ============================================================

itrf = FramesFactory.getITRF(
    IERSConventions.IERS_2010,
    True
)

earth = OneAxisEllipsoid(
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    Constants.WGS84_EARTH_FLATTENING,
    itrf
)

# ============================================================
# SOLAR ACTIVITY MODEL
# ============================================================

msafe = MarshallSolarActivityFutureEstimation(
    MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES,
    StrengthLevel.AVERAGE
)

from org.orekit.bodies import CelestialBodyFactory

sun = CelestialBodyFactory.getSun()

atmosphere = NRLMSISE00(
    msafe,
    sun,
    earth
)

# ============================================================
# DATE
# ============================================================

date = AbsoluteDate(
    2026,
    1,
    1,
    12,
    0,
    0.0,
    TimeScalesFactory.getUTC()
)

# ============================================================
# DENSITY FUNCTION
# ============================================================

def get_density(latitude_deg,
                longitude_deg,
                altitude_km):

    point = GeodeticPoint(
        np.radians(latitude_deg),
        np.radians(longitude_deg),
        altitude_km * 1000.0
    )

    position = earth.transform(point)

    density = atmosphere.getDensity(
        date,
        position,
        earth.getBodyFrame()
    )

    return density


# ============================================================
# COMPUTE DENSITIES
# ============================================================

print("\n")
print("===================================================")
print("NRLMSISE-00 ATMOSPHERIC DENSITY")
print("===================================================")

lat = 35.0
lon = 75.0

for altitude in [200, 550, 800]:

    rho = get_density(
        lat,
        lon,
        altitude
    )

    print(
        f"Altitude = {altitude:4.0f} km"
        f"   Density = {rho:.3e} kg/m³"
    )

print("===================================================")