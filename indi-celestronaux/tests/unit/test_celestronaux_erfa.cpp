/*
    Unit tests for the ErfaTelescopeMixin coordinate transformation pipeline.
    Covers erfaEqToHoriz, erfaHorizToEq, erfaJ2000ToTopocentric, localSiderealTime,
    atmospheric refraction magnitude and monotonicity, and solar tracking accuracy.

    Copyright (C) 2026 Christian Kemper

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <gtest/gtest.h>

#include "erfa_telescope_mixin.h"
#include "tracking_math.h"

using namespace tracking;

#include <libnova/julian_day.h>
#include <libnova/solar.h>
#include <libnova/lunar.h>
#include <libnova/transform.h>
#include <libnova/sidereal_time.h>
#include <libnova/utility.h>

#include <inditelescope.h>   // INDI::J2000toObserved

#include <cmath>

// ---------------------------------------------------------------------------
// Test subclass
// ---------------------------------------------------------------------------

class TestMixin : public ErfaTelescopeMixin
{
public:
    INDI::IGeographicCoordinates getObserverLocation() const override
    {
        // 47.5°N, 122.3°W, 100 m — Seattle area
        return { .latitude = 47.5, .longitude = -122.3, .elevation = 100.0 };
    }

    // Expose protected methods as public for tests.
    using ErfaTelescopeMixin::erfaEqToHoriz;
    using ErfaTelescopeMixin::erfaHorizToEq;
    using ErfaTelescopeMixin::erfaJ2000ToTopocentric;
    using ErfaTelescopeMixin::localSiderealTime;

    void setRefraction(bool on)
    {
        m_RefractionEnabled = on;
    }

    void setAtmosphere(double temp_C, double pressure_hPa, double humidity_pct, double wl_um)
    {
        m_Atmosphere.temp_C       = temp_C;
        m_Atmosphere.pressure_hPa = pressure_hPa;
        m_Atmosphere.humidity     = humidity_pct;
        m_Atmosphere.wavelength_um = wl_um;
    }
};

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------

class ErfaMixinTest : public ::testing::Test
{
protected:
    TestMixin mixin;
    // JD for 2024-06-21 06:00 UTC (summer solstice — high Dec for northern stars)
    const double JD_TEST = 2460482.75;

    void SetUp() override
    {
        mixin.setRefraction(false);
        mixin.setAtmosphere(15.0, 1013.25, 50.0, 0.55);
    }
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Angular separation in arcsec between two Alt/Az positions.
static double sepArcsec(double alt1, double az1, double alt2, double az2)
{
    double da  = (alt2 - alt1) * M_PI / 180.0;
    double daz = (az2  - az1)  * M_PI / 180.0;
    double cosAlt = std::cos((alt1 + alt2) / 2.0 * M_PI / 180.0);
    return std::sqrt(da * da + (daz * cosAlt) * (daz * cosAlt)) * 180.0 / M_PI * 3600.0;
}

// Angular separation in arcsec between two RA/Dec positions.
static double raDecSepArcsec(double ra1_h, double dec1_d, double ra2_h, double dec2_d)
{
    double dra  = (ra2_h  - ra1_h)  * 15.0;   // degrees
    double ddec =  dec2_d - dec1_d;
    double cosDec = std::cos((dec1_d + dec2_d) / 2.0 * M_PI / 180.0);
    return std::sqrt((dra * cosDec) * (dra * cosDec) + ddec * ddec) * 3600.0;
}

// ---------------------------------------------------------------------------
// T6 — Round-trip at multiple sky positions (no refraction)
// ---------------------------------------------------------------------------

struct SkyPos { double ra_h; double dec_d; const char *label; };

class RoundTripTest : public ErfaMixinTest,
                      public ::testing::WithParamInterface<SkyPos>
{};

TEST_P(RoundTripTest, EqHorizEq)
{
    auto pos = GetParam();
    // Derive RA from LST and the desired HA so that we can observe the star.
    double lst_h = mixin.localSiderealTime(JD_TEST, -122.3);
    double ra_h  = pos.ra_h; // pass raw RA; some positions may be below horizon but that's OK for math

    double altOut, azOut;
    mixin.erfaEqToHoriz(ra_h, pos.dec_d, JD_TEST, &altOut, &azOut);

    double raBack_h, decBack_d;
    mixin.erfaHorizToEq(altOut, azOut, JD_TEST, &raBack_h, &decBack_d);

    double sep = raDecSepArcsec(ra_h, pos.dec_d, raBack_h, decBack_d);
    EXPECT_LT(sep, 0.001) << "Round-trip error for " << pos.label << ": " << sep << " arcsec";

    (void)lst_h; // used only for documentation
}

INSTANTIATE_TEST_SUITE_P(
    CanonicalPositions, RoundTripTest,
    ::testing::Values(
        SkyPos{ 12.0,   0.0,  "Meridian equator" },
        SkyPos{ 12.0,  60.0,  "Meridian Dec+60" },
        SkyPos{ 12.0,  85.0,  "Near north pole" },
        SkyPos{  6.0,  30.0,  "East HA=-6h" },
        SkyPos{ 18.0,  30.0,  "West HA=+6h" },
        SkyPos{ 15.0, -30.0,  "Southern star" }
    )
);

// ---------------------------------------------------------------------------
// T7 — Round-trip near horizon with refraction enabled
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, RoundTripHorizonRefraction)
{
    mixin.setRefraction(true);
    mixin.setAtmosphere(15.0, 1013.25, 50.0, 0.55);

    // Use a star that will be near the horizon for our observer.
    // Dec ≈ lat - 90 + 5 ≈ -37.5°  at transit puts star ~5° above horizon.
    double ra_h  = 12.0;
    double dec_d = -37.5;

    double altOut, azOut;
    mixin.erfaEqToHoriz(ra_h, dec_d, JD_TEST, &altOut, &azOut);

    double raBack_h, decBack_d;
    mixin.erfaHorizToEq(altOut, azOut, JD_TEST, &raBack_h, &decBack_d);

    double sep = raDecSepArcsec(ra_h, dec_d, raBack_h, decBack_d);
    // Near-horizon refraction is ~30 arcmin; the forward/inverse ERFA models are not
    // perfectly symmetric at low altitude — 5 arcsec round-trip error is acceptable.
    EXPECT_LT(sep, 5.0) << "Near-horizon round-trip with refraction: " << sep << " arcsec";
}

// ---------------------------------------------------------------------------
// T8 — Refraction magnitude at canonical altitudes
// ---------------------------------------------------------------------------

// Returns the refraction (arcsec) for a star placed at apparent altitude
// approxAlt_deg by finding an appropriate RA/Dec and comparing ERFA with/without
// refraction.
static double refractionAt(TestMixin &mixin, double jd, double approxAlt_deg)
{
    // Transit: RA = LST, Dec = lat - (90 - approxAlt) = approxAlt + lat - 90
    double lst   = mixin.localSiderealTime(jd, -122.3);
    double lat   = 47.5;
    double dec_d = approxAlt_deg + lat - 90.0;

    mixin.setRefraction(false);
    double altNoRef, azNoRef;
    mixin.erfaEqToHoriz(lst, dec_d, jd, &altNoRef, &azNoRef);

    mixin.setRefraction(true);
    double altRef, azRef;
    mixin.erfaEqToHoriz(lst, dec_d, jd, &altRef, &azRef);

    return (altRef - altNoRef) * 3600.0; // arcsec
}

TEST_F(ErfaMixinTest, RefractionMagnitude)
{
    mixin.setAtmosphere(15.0, 1013.25, 50.0, 0.55);

    // At ~30° altitude: USNO formula gives ≈ 1.7 arcmin = 102 arcsec.
    double ref30 = refractionAt(mixin, JD_TEST, 30.0);
    EXPECT_GT(ref30, 80.0)  << "Refraction at 30° too small: " << ref30 << " arcsec";
    EXPECT_LT(ref30, 130.0) << "Refraction at 30° too large: " << ref30 << " arcsec";

    // At ~45° altitude: ≈ 1 arcmin = 60 arcsec.
    double ref45 = refractionAt(mixin, JD_TEST, 45.0);
    EXPECT_GT(ref45, 40.0)  << "Refraction at 45° too small: " << ref45 << " arcsec";
    EXPECT_LT(ref45, 80.0)  << "Refraction at 45° too large: " << ref45 << " arcsec";

    // At ~89° altitude: ERFA's IAU2006 refraction at near-zenith is typically ~1 arcsec
    // (the standard atmosphere gives ~0.7 arcsec; numerical precision adds a little more).
    double ref89 = refractionAt(mixin, JD_TEST, 89.0);
    EXPECT_LT(ref89, 3.0) << "Refraction near zenith: " << ref89 << " arcsec";
}

// ---------------------------------------------------------------------------
// T9 — Refraction always increases apparent altitude
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, RefractionAlwaysIncreasesAlt)
{
    mixin.setAtmosphere(15.0, 1013.25, 50.0, 0.55);
    const double testAlts[] = { 5.0, 20.0, 45.0, 70.0 };

    for (double approxAlt : testAlts)
    {
        double ref = refractionAt(mixin, JD_TEST, approxAlt);
        EXPECT_GT(ref, 0.0)
            << "Refraction should be positive (upward) at alt=" << approxAlt
            << "°, got " << ref << " arcsec";
    }
}

// ---------------------------------------------------------------------------
// T10 — localSiderealTime vs libnova at multiple epochs
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, LocalSiderealTimeVsLibnova)
{
    const double lon = -122.3;
    // Four JDs spanning ~6 months.
    const double jds[] = { 2460389.0, 2460430.0, 2460482.75, 2460540.0 };

    for (double jd : jds)
    {
        double erfa_lst_h = mixin.localSiderealTime(jd, lon);

        // libnova: ln_get_apparent_sidereal_time returns GAST in hours.
        double gast_h = ln_get_apparent_sidereal_time(jd);
        double nova_lst_h = gast_h + lon / 15.0;
        // Normalize to [0,24)
        while (nova_lst_h <  0.0) nova_lst_h += 24.0;
        while (nova_lst_h >= 24.0) nova_lst_h -= 24.0;

        double diffSec = std::abs(erfa_lst_h - nova_lst_h) * 3600.0;
        // The two may differ by up to ±24h due to wrap — take the minimum.
        if (diffSec > 43200.0) diffSec = 86400.0 - diffSec;

        EXPECT_LT(diffSec, 4.0)  // 4 arcsec (≈ 0.27 ms) — ERFA vs older libnova theory
            << "LST difference at JD=" << jd << ": " << diffSec << " arcsec";
    }
}

// ---------------------------------------------------------------------------
// T11 — erfaEqToHoriz vs libnova ln_get_hrz_from_equ (no refraction)
// ---------------------------------------------------------------------------

struct StarPos { double ra_deg; double dec_deg; const char *name; };

class VsLibnovaTest : public ErfaMixinTest,
                      public ::testing::WithParamInterface<StarPos>
{};

TEST_P(VsLibnovaTest, AltAzVsLibnova)
{
    auto star = GetParam();
    mixin.setRefraction(false);

    double altErfa, azErfa;
    mixin.erfaEqToHoriz(star.ra_deg / 15.0, star.dec_deg, JD_TEST, &altErfa, &azErfa);

    // libnova needs apparent RA/Dec (JNow) and observer.
    ln_equ_posn eqPos { star.ra_deg, star.dec_deg };
    ln_lnlat_posn obs { -122.3, 47.5 };   // W longitude is negative in libnova
    ln_hrz_posn  hrz;
    ln_get_hrz_from_equ(&eqPos, &obs, JD_TEST, &hrz);
    // libnova azimuth is south-based (0=S); ERFA is north-based (0=N).
    // Convert libnova az to north-based: az_N = (az_S + 180) mod 360.
    double azNova = std::fmod(hrz.az + 180.0, 360.0);
    double altNova = hrz.alt;

    double sep = sepArcsec(altErfa, azErfa, altNova, azNova);
    EXPECT_LT(sep, 30.0)
        << "ERFA vs libnova for " << star.name << ": " << sep << " arcsec";
}

INSTANTIATE_TEST_SUITE_P(
    KnownStars, VsLibnovaTest,
    ::testing::Values(
        StarPos{ 279.23,  38.78,  "Vega"    },   // RA 18h36m, Dec +38°47'
        StarPos{  79.17,  45.99,  "Capella" },   // RA 5h17m,  Dec +45°59'
        StarPos{ 101.29, -16.72,  "Sirius"  },   // RA 6h45m,  Dec -16°43'
        StarPos{ 213.92,  19.18,  "Arcturus"}    // RA 14h16m, Dec +19°11'
    )
);

// ---------------------------------------------------------------------------
// T12 — erfaJ2000ToTopocentric includes annual aberration
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, J2000TopocentricIncludesAberration)
{
    // Use Vega (J2000): RA=18h36m23s, Dec=+38°47'01"
    double ra_j2000 = 279.234 / 15.0; // hours
    double dec_j2000 = 38.784;        // degrees

    // Two JDs six months apart (aberration changes sign).
    double JD_A = 2460389.0; // March equinox
    double JD_B = 2460569.0; // September equinox (~6 months later)

    double ra_A, dec_A, ra_B, dec_B;
    mixin.erfaJ2000ToTopocentric(ra_j2000, dec_j2000, JD_A, &ra_A, &dec_A);
    mixin.erfaJ2000ToTopocentric(ra_j2000, dec_j2000, JD_B, &ra_B, &dec_B);

    double sep = raDecSepArcsec(ra_A, dec_A, ra_B, dec_B);

    // Annual aberration semi-major axis ≈ 20.5 arcsec; full ellipse diameter ≈ 41 arcsec.
    EXPECT_GT(sep, 10.0) << "Aberration too small between equinoxes: " << sep << " arcsec";
    EXPECT_LT(sep, 60.0) << "Aberration too large between equinoxes: " << sep << " arcsec";
}

// ---------------------------------------------------------------------------
// T13 — erfaJ2000ToTopocentric vs INDI::J2000toObserved
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, J2000TopocentricVsINDI)
{
    const StarPos stars[] = {
        { 279.23,  38.78, "Vega"    },
        {  79.17,  45.99, "Capella" },
        { 213.92,  19.18, "Arcturus"},
    };

    for (auto &star : stars)
    {
        double ra_j2000_h = star.ra_deg / 15.0;
        double erfa_ra, erfa_dec;
        mixin.erfaJ2000ToTopocentric(ra_j2000_h, star.dec_deg, JD_TEST, &erfa_ra, &erfa_dec);

        INDI::IEquatorialCoordinates j2000Eq { ra_j2000_h, star.dec_deg };
        INDI::IEquatorialCoordinates observed;
        INDI::J2000toObserved(&j2000Eq, JD_TEST, &observed);

        double sep = raDecSepArcsec(erfa_ra, erfa_dec,
                                    observed.rightascension, observed.declination);
        EXPECT_LT(sep, 10.0)
            << "ERFA vs INDI::J2000toObserved for " << star.name << ": " << sep << " arcsec";
    }
}

// ---------------------------------------------------------------------------
// T14 — Solar tracking simulation with ERFA conversions (10 minutes)
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, SolarTrackingSimulation10Min)
{
    const int    ticks  = 600;   // 10 minutes at 1 s/tick
    const double dt_day = 1.0 / 86400.0;

    double JD     = JD_TEST;
    bool   primed = false;
    double lastRA = 0, lastDec = 0;
    double targetRA, targetDec;

    {
        ln_equ_posn p0;
        ln_get_solar_equ_coords(JD, &p0);
        targetRA  = p0.ra / 15.0;
        targetDec = p0.dec;
    }

    for (int i = 0; i < ticks; ++i)
    {
        JD += dt_day;
        ln_equ_posn pos;
        ln_get_solar_equ_coords(JD, &pos);
        double newRA  = pos.ra / 15.0;
        double newDec = pos.dec;

        if (primed)
        {
            targetRA += wrapHA(newRA - lastRA);
            while (targetRA <  0.0) targetRA += 24.0;
            while (targetRA >= 24.0) targetRA -= 24.0;
            targetDec += newDec - lastDec;
        }
        lastRA  = newRA;
        lastDec = newDec;
        primed  = true;

        // Convert both the accumulated target and the direct solar position to Alt/Az.
        double altTarget, azTarget;
        mixin.erfaEqToHoriz(targetRA, targetDec, JD, &altTarget, &azTarget);

        double altTruth, azTruth;
        mixin.erfaEqToHoriz(newRA, newDec, JD, &altTruth, &azTruth);

        double sep = sepArcsec(altTarget, azTarget, altTruth, azTruth);
        EXPECT_LT(sep, 5.0)
            << "Solar ERFA tracking error at tick " << i + 1 << ": " << sep << " arcsec";
    }
}

// ---------------------------------------------------------------------------
// T15 — Solar feature tracking with ERFA: offset maintenance (10 minutes)
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, SolarFeatureOffsetMaintained)
{
    const int    ticks  = 600;
    const double dt_day = 1.0 / 86400.0;
    // Feature offset: 10 arcmin north and 10 arcmin west of solar center.
    const double OFFSET_ARCMIN = 10.0;

    double JD    = JD_TEST;
    bool   primed = false;
    double lastRA = 0, lastDec = 0;
    double targetRA, targetDec;

    {
        ln_equ_posn p0;
        ln_get_solar_equ_coords(JD, &p0);
        double sunRA  = p0.ra / 15.0;
        double sunDec = p0.dec;
        double cosDec = std::cos(sunDec * M_PI / 180.0);
        targetRA  = sunRA  - (OFFSET_ARCMIN / 60.0 / 15.0 / cosDec);
        targetDec = sunDec + (OFFSET_ARCMIN / 60.0);
        lastRA  = sunRA;
        lastDec = sunDec;
        primed  = false;
    }

    for (int i = 0; i < ticks; ++i)
    {
        JD += dt_day;
        ln_equ_posn pos;
        ln_get_solar_equ_coords(JD, &pos);
        double newRA  = pos.ra / 15.0;
        double newDec = pos.dec;

        if (primed)
        {
            targetRA += wrapHA(newRA - lastRA);
            while (targetRA <  0.0) targetRA += 24.0;
            while (targetRA >= 24.0) targetRA -= 24.0;
            targetDec += newDec - lastDec;
        }
        lastRA  = newRA;
        lastDec = newDec;
        primed  = true;

        // Current angular offset from the solar center (arcmin).
        double cosDec = std::cos(newDec * M_PI / 180.0);
        double dRAdeg = (targetRA  - newRA)  * 15.0;
        double dDecdeg = targetDec - newDec;
        double raOffsetArcmin  = dRAdeg  * cosDec * 60.0;
        double decOffsetArcmin = dDecdeg * 60.0;

        // The RA offset shrinks slightly as Dec changes (converging meridians).
        // Over 10 minutes the geometric drift is < 0.01 arcmin — use 0.05 arcmin tolerance.
        EXPECT_NEAR(raOffsetArcmin,  -OFFSET_ARCMIN, 0.05)
            << "RA feature offset drifted at tick " << i + 1;
        EXPECT_NEAR(decOffsetArcmin,  OFFSET_ARCMIN, 0.05)
            << "Dec feature offset drifted at tick " << i + 1;
    }
}

// ---------------------------------------------------------------------------
// T16 — Atmospheric parameters wired correctly to eraApco13
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, AtmosphericPressureMonotonic)
{
    // For a star near 10° altitude, refraction should increase with pressure.
    double lst   = mixin.localSiderealTime(JD_TEST, -122.3);
    double lat   = 47.5;
    double dec_d = 10.0 + lat - 90.0; // star at ~10° alt at transit

    const double pressures[] = { 800.0, 900.0, 1013.25 };
    double prevAlt = -1e9;
    for (double p : pressures)
    {
        mixin.setRefraction(true);
        mixin.setAtmosphere(15.0, p, 50.0, 0.55);
        double alt, az;
        mixin.erfaEqToHoriz(lst, dec_d, JD_TEST, &alt, &az);
        EXPECT_GT(alt, prevAlt)
            << "Altitude did not increase with pressure " << p << " hPa";
        prevAlt = alt;
    }
}

TEST_F(ErfaMixinTest, AtmosphericTemperatureMonotonic)
{
    // Refraction decreases with temperature (hotter air is less dense).
    double lst   = mixin.localSiderealTime(JD_TEST, -122.3);
    double lat   = 47.5;
    double dec_d = 10.0 + lat - 90.0;

    const double temps[] = { -20.0, 15.0, 40.0 };
    double prevAlt = 1e9;
    for (double t : temps)
    {
        mixin.setRefraction(true);
        mixin.setAtmosphere(t, 1013.25, 50.0, 0.55);
        double alt, az;
        mixin.erfaEqToHoriz(lst, dec_d, JD_TEST, &alt, &az);
        EXPECT_LT(alt, prevAlt)
            << "Altitude did not decrease with temperature " << t << " °C";
        prevAlt = alt;
    }
}

// ---------------------------------------------------------------------------
// T17 — Disabling refraction overrides stored pressure
// ---------------------------------------------------------------------------
TEST_F(ErfaMixinTest, RefractionDisabledIgnoresPressure)
{
    double ra_h  = 12.0;
    double dec_d = 30.0;

    // With refraction off but high pressure stored.
    mixin.setAtmosphere(15.0, 1013.25, 50.0, 0.55);
    mixin.setRefraction(false);
    double altOff, azOff;
    mixin.erfaEqToHoriz(ra_h, dec_d, JD_TEST, &altOff, &azOff);

    // With refraction on and zero pressure (manual no-refraction reference).
    mixin.setAtmosphere(15.0, 0.0, 50.0, 0.55);
    mixin.setRefraction(true);
    double altZeroPressure, azZeroPressure;
    mixin.erfaEqToHoriz(ra_h, dec_d, JD_TEST, &altZeroPressure, &azZeroPressure);

    // Both should give the same result (no refraction in both cases).
    EXPECT_NEAR(altOff, altZeroPressure, 0.001)
        << "Refraction=off should match pressure=0: diff "
        << (altOff - altZeroPressure) * 3600.0 << " arcsec";
    EXPECT_NEAR(azOff, azZeroPressure, 0.001)
        << "Refraction=off az mismatch";
}

// ---------------------------------------------------------------------------

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
