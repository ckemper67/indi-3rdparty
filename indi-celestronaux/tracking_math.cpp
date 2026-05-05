/*
    Pure-arithmetic primitives for parabolic interpolation and motor rate computation.

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

#include "tracking_math.h"
#include <cmath>
#include <tuple>

namespace tracking {

double wrapAngle(double angle)
{
    while (angle <= -180.0) angle += 360.0;
    while (angle >   180.0) angle -= 360.0;
    return angle;
}

double wrapHA(double ha)
{
    while (ha <= -12.0) ha += 24.0;
    while (ha >   12.0) ha -= 24.0;
    return ha;
}

// ---------------------------------------------------------------------------
// QuadraticInterpolator
// ---------------------------------------------------------------------------

QuadraticInterpolator::QuadraticInterpolator(SampleFn fn, double stepJD, double JDnow)
    : m_fn(std::move(fn)), m_step(stepJD)
{
    prime(JDnow);
}

void QuadraticInterpolator::prime(double JDnow)
{
    m_JD[0] = JDnow - m_step;
    m_JD[1] = JDnow;
    m_JD[2] = JDnow + m_step;
    for (int i = 0; i < 3; ++i)
        std::tie(m_x[i], m_y[i]) = m_fn(m_JD[i]);
    m_ready = true;
}

void QuadraticInterpolator::advance()
{
    m_x[0]  = m_x[1];
    m_y[0]  = m_y[1];
    m_JD[0] = m_JD[1];
    m_x[1]  = m_x[2];
    m_y[1]  = m_y[2];
    m_JD[1] = m_JD[2];
    m_JD[2] += m_step;
    std::tie(m_x[2], m_y[2]) = m_fn(m_JD[2]);
}

// Parabolic coefficients for a 3-point window centred at x0/y0.
// Returns (x_b, x_a, y_b, y_a): value(t) = c + b*t + a*t^2, t in seconds from centre.
// x is unwrapped with wrapAngle to handle azimuth discontinuities.
static std::tuple<double, double, double, double> coefficients(
    double xM, double x0, double xP,
    double yM, double y0, double yP,
    double T)
{
    xM = x0 + wrapAngle(xM - x0);
    xP = x0 + wrapAngle(xP - x0);

    double x_b = (xP - xM)           / (2 * T);
    double x_a = (xP + xM - 2 * x0) / (2 * T * T);
    double y_b = (yP - yM)           / (2 * T);
    double y_a = (yP + yM - 2 * y0) / (2 * T * T);
    return std::make_tuple(x_b, x_a, y_b, y_a);
}

void QuadraticInterpolator::advanceTo(double JD)
{
    if ((JD - m_JD[2]) > 3.0 * m_step)
        prime(JD);
    else
        while (JD >= m_JD[2])
            advance();
}

std::pair<double, double> QuadraticInterpolator::valueAt(double JD)
{
    advanceTo(JD);

    double T  = m_step * 86400.0;
    double dt = (JD - m_JD[1]) * 86400.0;
    auto [x_b, x_a, y_b, y_a] = coefficients(
        m_x[0], m_x[1], m_x[2],
        m_y[0], m_y[1], m_y[2], T);

    return std::make_pair(
        wrapAngle(m_x[1] + x_b * dt + x_a * dt * dt),
        m_y[1] + y_b * dt + y_a * dt * dt
    );
}

std::pair<double, double> QuadraticInterpolator::rateAt(double JD)
{
    advanceTo(JD);

    double T  = m_step * 86400.0;
    double dt = (JD - m_JD[1]) * 86400.0;
    auto [x_b, x_a, y_b, y_a] = coefficients(
        m_x[0], m_x[1], m_x[2],
        m_y[0], m_y[1], m_y[2], T);

    return std::make_pair(
        x_b + 2 * x_a * dt,
        y_b + 2 * y_a * dt
    );
}

void computeEqTrackRates(double ra0_deg,  double dec0_deg,
                         double ra1_deg,  double dec1_deg,
                         double dt_sec,
                         double *raRate_arcsec_s, double *decRate_arcsec_s)
{
    double dRA_arcsec_s  = (ra1_deg  - ra0_deg)  * 3600.0 / dt_sec;
    double dDec_arcsec_s = (dec1_deg - dec0_deg) * 3600.0 / dt_sec;
    *raRate_arcsec_s  = SIDEREAL_RATE - dRA_arcsec_s;
    *decRate_arcsec_s = dDec_arcsec_s;
}

} // namespace tracking
