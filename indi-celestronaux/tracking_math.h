/*
    Pure-arithmetic primitives for parabolic interpolation and motor rate computation.
    No INDI or ephemeris dependency — callers supply all sky positions.

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

#pragma once

#include <functional>
#include <utility>

namespace tracking {

// Sidereal tracking rate in arcsec/s — same value as INDI's TRACKRATE_SIDEREAL macro.
// Named SIDEREAL_RATE here to avoid colliding with that macro in files that include indicom.h.
constexpr double SIDEREAL_RATE = 15.041067;

// Normalize angle to (-180, +180].
double wrapAngle(double angle);

// Normalize hour angle to (-12, +12].
double wrapHA(double ha);

// ---------------------------------------------------------------------------
// QuadraticInterpolator — 3-point sliding parabolic interpolator
//
// Maps a caller-supplied function (JD → pair<x, y>) onto a sliding window of
// three evenly-spaced JD samples.  x is treated as angular (wrapped to
// (-180, 180] relative to the centre sample); y is linear.
//
// Typical use: pass a lambda that converts a tracking target RA/Dec to
// (Az_deg, Alt_deg) at the given JD.  Construct with fn, stepJD, and the
// current JD; then call valueAt/rateAt each tick.
// ---------------------------------------------------------------------------
class QuadraticInterpolator
{
public:
    using SampleFn = std::function<std::pair<double, double>(double JD)>;

    QuadraticInterpolator() = default;
    QuadraticInterpolator(SampleFn fn, double stepJD, double JDnow);

    // Parabolically-interpolated (x, y) at JD.
    // Advances the window lazily; re-primes if JD is more than 3 steps ahead.
    std::pair<double, double> valueAt(double JD);

    // First derivative (dx/dt, dy/dt) at JD in units per second.
    std::pair<double, double> rateAt(double JD);

    bool   isReady() const { return m_ready; }
    void   reset()         { m_ready = false; }
    double step()    const { return m_step; }

private:
    void prime(double JDnow);
    void advance();
    void advanceTo(double JD);

    SampleFn       m_fn;
    double         m_step  { 0 };
    mutable double m_x[3]  { 0, 0, 0 };
    mutable double m_y[3]  { 0, 0, 0 };
    mutable double m_JD[3] { 0, 0, 0 };
    bool           m_ready { false };
};

// ---------------------------------------------------------------------------
// Motor rate helpers
// ---------------------------------------------------------------------------

// EQ motor rates from a two-point positional derivative.
//   ra{0,1}_deg, dec{0,1}_deg : positions at t and t+dt (degrees, as libnova returns)
//   dt_sec                    : time step in seconds
//   raRate_arcsec_s           : SIDEREAL_RATE − dRA/dt
//   decRate_arcsec_s          : dDec/dt
void computeEqTrackRates(double ra0_deg,  double dec0_deg,
                         double ra1_deg,  double dec1_deg,
                         double dt_sec,
                         double *raRate_arcsec_s, double *decRate_arcsec_s);

} // namespace tracking
