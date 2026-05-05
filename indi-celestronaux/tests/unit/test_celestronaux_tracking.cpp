/*
    Unit tests for the tracking math primitives in tracking_math.cpp.
    Covers ephemeris delta accumulation, EQ motor rate computation, and
    parabolic AltAz steering.  Links only against tracking_math.cpp and libnova;
    no INDI infrastructure required.

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

#include "tracking_math.h"

using namespace tracking;

#include <libnova/julian_day.h>
#include <libnova/solar.h>
#include <libnova/lunar.h>

#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static double arcsecError(double raErrH, double decErrDeg)
{
    double raDeg  = raErrH * 15.0;
    double errArc = std::sqrt(raDeg * raDeg + decErrDeg * decErrDeg) * 3600.0;
    return errArc;
}

// JD for 2024-03-20 12:00 UTC (near vernal equinox — low solar Dec rate)
static const double JD_EQUINOX = 2460389.0;

// ---------------------------------------------------------------------------
// T1 — AltAz solar delta accumulation over 1 hour
// ---------------------------------------------------------------------------
TEST(EphemDelta, SolarOneHour)
{
    const int    ticks  = 3600;
    const double dt_s   = 1.0;
    const double dt_day = dt_s / 86400.0;

    double JD  = JD_EQUINOX;
    bool   primed  = false;
    double lastRA  = 0, lastDec = 0;
    double targetRA, targetDec;

    // Seed the target with the initial solar position.
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
    }

    // Compare accumulated target with direct libnova query at the end JD.
    ln_equ_posn truth;
    ln_get_solar_equ_coords(JD, &truth);
    double truthRA  = truth.ra / 15.0;
    double truthDec = truth.dec;

    double errRA  = truthRA  - targetRA;
    // Wrap RA error to (-12, 12]
    while (errRA <= -12.0) errRA += 24.0;
    while (errRA >   12.0) errRA -= 24.0;
    double errDec = truthDec - targetDec;
    double errArc = arcsecError(errRA, errDec);

    EXPECT_LT(errArc, 1.0) << "Solar delta drift after 1 hour: " << errArc << " arcsec";
}

// ---------------------------------------------------------------------------
// T2 — AltAz lunar delta accumulation over 1 hour
// ---------------------------------------------------------------------------
TEST(EphemDelta, LunarOneHour)
{
    const int    ticks  = 3600;
    const double dt_s   = 1.0;
    const double dt_day = dt_s / 86400.0;

    double JD = JD_EQUINOX;
    bool   primed  = false;
    double lastRA  = 0, lastDec = 0;
    double targetRA, targetDec;

    {
        ln_equ_posn p0;
        ln_get_lunar_equ_coords(JD, &p0);
        targetRA  = p0.ra / 15.0;
        targetDec = p0.dec;
    }

    for (int i = 0; i < ticks; ++i)
    {
        JD += dt_day;
        ln_equ_posn pos;
        ln_get_lunar_equ_coords(JD, &pos);
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
    }

    ln_equ_posn truth;
    ln_get_lunar_equ_coords(JD, &truth);
    double truthRA  = truth.ra / 15.0;
    double truthDec = truth.dec;

    double errRA  = truthRA  - targetRA;
    while (errRA <= -12.0) errRA += 24.0;
    while (errRA >   12.0) errRA -= 24.0;
    double errDec = truthDec - targetDec;
    double errArc = arcsecError(errRA, errDec);

    // Moon moves ~0.5 deg/hr; per-tick rounding accumulates more than solar.
    EXPECT_LT(errArc, 5.0) << "Lunar delta drift after 1 hour: " << errArc << " arcsec";
}

// ---------------------------------------------------------------------------
// T3 — EQ solar rate sign and magnitude near vernal equinox
// ---------------------------------------------------------------------------
TEST(EqTrackRates, SolarSignAndMagnitude)
{
    const double dt = 1.0; // seconds
    double JD0 = JD_EQUINOX;
    double JD1 = JD0 + dt / 86400.0;

    ln_equ_posn pos0, pos1;
    ln_get_solar_equ_coords(JD0, &pos0);
    ln_get_solar_equ_coords(JD1, &pos1);

    double raRate, decRate;
    computeEqTrackRates(pos0.ra, pos0.dec, pos1.ra, pos1.dec, dt, &raRate, &decRate);

    // Solar RA drifts eastward ~0.04 arcsec/s => raRate ≈ sidereal − 0.04
    double expectedRA = SIDEREAL_RATE - (pos1.ra - pos0.ra) * 3600.0 / dt;
    EXPECT_NEAR(raRate, expectedRA, 0.001)
        << "Solar RA rate: " << raRate << " arcsec/s (expected ~" << expectedRA << ")";

    // Near equinox the Sun's Dec barely changes.
    EXPECT_LT(std::abs(decRate), 0.02)
        << "Solar Dec rate at equinox: " << decRate << " arcsec/s";

    // raRate must be less than sidereal (Sun drifts east, motor runs slower).
    EXPECT_LT(raRate, SIDEREAL_RATE)
        << "Solar RA rate should be below sidereal";
}

// ---------------------------------------------------------------------------
// T4 — EQ lunar rate sign and magnitude
// ---------------------------------------------------------------------------
TEST(EqTrackRates, LunarSignAndMagnitude)
{
    const double dt = 1.0;
    double JD0 = JD_EQUINOX;
    double JD1 = JD0 + dt / 86400.0;

    ln_equ_posn pos0, pos1;
    ln_get_lunar_equ_coords(JD0, &pos0);
    ln_get_lunar_equ_coords(JD1, &pos1);

    double raRate, decRate;
    computeEqTrackRates(pos0.ra, pos0.dec, pos1.ra, pos1.dec, dt, &raRate, &decRate);

    // Moon drifts east ~0.55 arcsec/s => raRate ≈ sidereal − 0.55
    double expectedRA = SIDEREAL_RATE - (pos1.ra - pos0.ra) * 3600.0 / dt;
    EXPECT_NEAR(raRate, expectedRA, 0.001)
        << "Lunar RA rate: " << raRate << " arcsec/s";

    // Verify the lunar rate is meaningfully lower than sidereal (Moon moves east).
    EXPECT_LT(raRate, SIDEREAL_RATE - 0.3)
        << "Lunar RA rate should be well below sidereal";
}

// ---------------------------------------------------------------------------
// T5 — QuadraticInterpolator on synthetic circumpolar path (1 hour)
// ---------------------------------------------------------------------------

// Synthetic circumpolar truth: 5-deg-radius circle at 45° altitude.
static const double CIRC_R     = 5.0;
static const double CIRC_H0    = 45.0;
static const double CIRC_OMEGA = (15.0 / 3600.0) * (M_PI / 180.0); // rad/s

static std::pair<double,double> circumpolarTruth(double t_s)
{
    return { CIRC_R * std::cos(CIRC_OMEGA * t_s),
             CIRC_H0 + CIRC_R * std::sin(CIRC_OMEGA * t_s) };
}

static double range180_local(double a)
{
    while (a <= -180.0) a += 360.0;
    while (a >   180.0) a -= 360.0;
    return a;
}

TEST(QuadraticInterpolator, SyntheticCircumpolar)
{
    const double dt_s   = 1.0;
    const double dt_day = dt_s / 86400.0;
    const double JD0    = 2460389.0;

    // SampleFn: JD → (Az, Alt) from the synthetic circumpolar truth.
    auto sampleFn = [&](double JD) -> std::pair<double, double> {
        double t = (JD - JD0) * 86400.0;
        return circumpolarTruth(t);
    };

    // Linear model: backward-diff rate applied to accumulated position.
    auto start  = circumpolarTruth(0.0);
    double linAz  = start.first,  linAlt  = start.second;
    double maxErrLinear    = 0.0;
    double maxErrParabolic = 0.0;

    QuadraticInterpolator win(sampleFn, dt_day, JD0);

    // Parabolic model: steering from accumulated mount position to predicted next.
    double paraAz = start.first, paraAlt = start.second;

    for (int tick = 0; tick < 3600; ++tick)
    {
        double JDnow  = JD0 + tick * dt_day;
        double JDnext = JDnow + dt_day;

        // Linear model: backward-diff velocity from truth.
        auto pPrev = circumpolarTruth((tick - 1) * dt_s);
        auto pNow  = circumpolarTruth(tick * dt_s);
        double linVAz  = range180_local(pNow.first  - pPrev.first)  / dt_s;
        double linVAlt = (pNow.second - pPrev.second) / dt_s;
        linAz  += linVAz  * dt_s;
        linAlt += linVAlt * dt_s;

        // Parabolic model: steer accumulated position towards predicted next tick.
        auto [azNext, altNext] = win.valueAt(JDnext);
        double vAz  = wrapAngle(azNext - paraAz)  / dt_s;
        double vAlt = (altNext - paraAlt) / dt_s;
        paraAz  += vAz  * dt_s;
        paraAlt += vAlt * dt_s;

        // Ground truth one step ahead.
        auto truth  = circumpolarTruth((tick + 1) * dt_s);
        double cosAlt = std::cos(truth.second * M_PI / 180.0);

        double linErrAz  = range180_local(linAz  - truth.first)  * cosAlt;
        double linErrAlt = linAlt - truth.second;
        maxErrLinear     = std::max(maxErrLinear,
            std::sqrt(linErrAz * linErrAz + linErrAlt * linErrAlt) * 3600.0);

        double paraErrAz  = wrapAngle(paraAz  - truth.first)  * cosAlt;
        double paraErrAlt = paraAlt - truth.second;
        maxErrParabolic   = std::max(maxErrParabolic,
            std::sqrt(paraErrAz * paraErrAz + paraErrAlt * paraErrAlt) * 3600.0);
    }

    EXPECT_LT(maxErrParabolic, 1.0)
        << "Parabolic max error: " << maxErrParabolic << " arcsec";
    EXPECT_GT(maxErrLinear, 0.1)
        << "Linear model should show >0.1 arcsec error (discriminating power check)";
}

// ---------------------------------------------------------------------------

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
