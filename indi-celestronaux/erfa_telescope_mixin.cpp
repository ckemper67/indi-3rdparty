/*
    ErfaTelescopeMixin — implementation.

    Copyright (C) 2026 Christian Kemper

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.
*/

#include "erfa_telescope_mixin.h"

#include <indicom.h>
#include <defaultdevice.h>

#include <chrono>
#include <cmath>

// ---------------------------------------------------------------------------
// Time helpers
// ---------------------------------------------------------------------------

double ErfaTelescopeMixin::currentJD() const
{
    using namespace std::chrono;
    auto now = system_clock::now();
    auto epoch_ms = duration_cast<milliseconds>(now.time_since_epoch()).count();
    // Unix epoch = 2440587.5 JD
    return 2440587.5 + epoch_ms / 86400000.0;
}

double ErfaTelescopeMixin::localSiderealTime(double jd, double longitude_deg) const
{
    double utc1 = std::floor(jd - 0.5) + 0.5;
    double utc2 = jd - utc1;
    double tai1, tai2, tt1, tt2;
    eraUtctai(utc1, utc2, &tai1, &tai2);
    eraTaitt(tai1, tai2, &tt1, &tt2);
    double gst = eraGst06a(utc1, utc2, tt1, tt2);   // GAST in radians
    double lst = eraAnp(gst + longitude_deg * ERFA_DD2R);
    return lst * ERFA_DR2H;                          // radians → hours
}

// ---------------------------------------------------------------------------
// ERFA context builder
// ---------------------------------------------------------------------------

void ErfaTelescopeMixin::buildErfaContext(double jd, eraASTROM *astrom, double *eo) const
{
    auto loc = getObserverLocation();
    double lon_rad = loc.longitude * ERFA_DD2R;
    double lat_rad = loc.latitude  * ERFA_DD2R;
    double elev_m  = loc.elevation;

    double utc1 = std::floor(jd - 0.5) + 0.5;
    double utc2 = jd - utc1;

    double phpa = m_RefractionEnabled ? m_Atmosphere.pressure_hPa : 0.0;
    double tc   = m_Atmosphere.temp_C;
    double rh   = m_Atmosphere.humidity / 100.0;   // percent → fraction
    double wl   = m_Atmosphere.wavelength_um;

    eraApco13(utc1, utc2, 0.0,
              lon_rad, lat_rad, elev_m,
              0.0, 0.0,
              phpa, tc, rh, wl,
              astrom, eo);
}

// ---------------------------------------------------------------------------
// Coordinate conversions
// ---------------------------------------------------------------------------

void ErfaTelescopeMixin::erfaEqToHoriz(double ra_hours, double dec_deg, double jd,
                                        double *alt_deg, double *az_deg) const
{
    eraASTROM astrom;
    double eo;
    buildErfaContext(jd, &astrom, &eo);

    // Convert apparent RA (hours) → CIRS RA (rad) via equation of the origins
    double ra_cirs = eraAnp(ra_hours * ERFA_DH2R + eo);
    double dec_rad = dec_deg * ERFA_DD2R;

    double aob, zob, hob, dob, rob;
    eraAtioq(ra_cirs, dec_rad, &astrom, &aob, &zob, &hob, &dob, &rob);

    *az_deg  = aob * ERFA_DR2D;
    *alt_deg = (ERFA_DPI / 2.0 - zob) * ERFA_DR2D;
}

void ErfaTelescopeMixin::erfaHorizToEq(double alt_deg, double az_deg, double jd,
                                        double *ra_hours, double *dec_deg) const
{
    eraASTROM astrom;
    double eo;
    buildErfaContext(jd, &astrom, &eo);

    double az_rad = az_deg * ERFA_DD2R;
    double zd_rad = (90.0 - alt_deg) * ERFA_DD2R;

    double ri, di;
    eraAtoiq("A", az_rad, zd_rad, &astrom, &ri, &di);

    // CIRS RA → apparent RA via equation of the origins
    *ra_hours = eraAnp(ri - eo) * ERFA_DR2H;
    *dec_deg  = di * ERFA_DR2D;
}

void ErfaTelescopeMixin::erfaJ2000ToTopocentric(double ra_j2000_hours, double dec_j2000_deg, double jd,
                                                 double *ra_topo_hours, double *dec_topo_deg) const
{
    eraASTROM astrom;
    double eo;
    buildErfaContext(jd, &astrom, &eo);

    double ra_rad  = ra_j2000_hours * ERFA_DH2R;
    double dec_rad = dec_j2000_deg  * ERFA_DD2R;

    double ri, di;
    // eraAtciq applies annual + diurnal aberration (no refraction)
    eraAtciq(ra_rad, dec_rad, 0.0, 0.0, 0.0, 0.0, &astrom, &ri, &di);

    *ra_topo_hours = eraAnp(ri - eo) * ERFA_DR2H;
    *dec_topo_deg  = di * ERFA_DR2D;
}

// ---------------------------------------------------------------------------
// Property lifecycle
// ---------------------------------------------------------------------------

void ErfaTelescopeMixin::initErfaProperties(const char *groupTab)
{
    // downcast to access getDeviceName() — valid because ErfaTelescopeMixin is
    // always mixed into an INDI::DefaultDevice subclass
    auto *dev = dynamic_cast<INDI::DefaultDevice *>(this);
    const char *devName = dev ? dev->getDeviceName() : "";

    AtmosphericNP[ATMOS_TEMPERATURE].fill("TEMPERATURE", "Temperature (C)", "%.1f", -50, 50, 1, 15.0);
    AtmosphericNP[ATMOS_PRESSURE].fill("PRESSURE", "Pressure (hPa)", "%.1f", 500, 1100, 1, 1013.25);
    AtmosphericNP[ATMOS_HUMIDITY].fill("HUMIDITY", "Humidity (%)", "%.0f", 0, 100, 1, 50.0);
    AtmosphericNP[ATMOS_WAVELENGTH].fill("WAVELENGTH", "Wavelength (um)", "%.2f", 0.2, 2.0, 0.01, 0.55);
    AtmosphericNP.fill(devName, "ATMOSPHERIC_CONDITIONS", "Atmosphere", groupTab, IP_RW, 60, IPS_IDLE);

    AtmosphericSourceTP[0].fill("ATMOS_SOURCE", "Source", "MANUAL");
    AtmosphericSourceTP.fill(devName, "ATMOS_SOURCE", "Atmos Source", groupTab, IP_RO, 60, IPS_IDLE);

    RefractionSP[REFRACTION_OFF].fill("REFRACTION_OFF", "Off", ISS_ON);
    RefractionSP[REFRACTION_SOFTWARE].fill("REFRACTION_SOFTWARE", "Software (ERFA)", ISS_OFF);
    RefractionSP.fill(devName, "REFRACTION_CONTROL", "Refraction", groupTab, IP_RW, ISR_1OFMANY, 60, IPS_IDLE);
}

void ErfaTelescopeMixin::updateErfaProperties(bool connected)
{
    auto *dev = dynamic_cast<INDI::DefaultDevice *>(this);
    if (!dev) return;
    if (connected)
    {
        dev->defineProperty(AtmosphericNP);
        dev->defineProperty(AtmosphericSourceTP);
        dev->defineProperty(RefractionSP);
    }
    else
    {
        dev->deleteProperty(AtmosphericNP.getName());
        dev->deleteProperty(AtmosphericSourceTP.getName());
        dev->deleteProperty(RefractionSP.getName());
    }
}

bool ErfaTelescopeMixin::handleErfaNumber(const char *name, double *values, char **names, int n)
{
    if (AtmosphericNP.isNameMatch(name))
    {
        AtmosphericNP.update(values, names, n);
        m_Atmosphere.temp_C        = AtmosphericNP[ATMOS_TEMPERATURE].getValue();
        m_Atmosphere.pressure_hPa  = AtmosphericNP[ATMOS_PRESSURE].getValue();
        m_Atmosphere.humidity      = AtmosphericNP[ATMOS_HUMIDITY].getValue();
        m_Atmosphere.wavelength_um = AtmosphericNP[ATMOS_WAVELENGTH].getValue();
        AtmosphericNP.setState(IPS_OK);
        AtmosphericNP.apply();
        return true;
    }
    return false;
}

bool ErfaTelescopeMixin::handleErfaSwitch(const char *name, ISState *states, char **names, int n)
{
    if (RefractionSP.isNameMatch(name))
    {
        RefractionSP.update(states, names, n);
        m_RefractionEnabled = (RefractionSP[REFRACTION_SOFTWARE].getState() == ISS_ON);
        RefractionSP.setState(IPS_OK);
        RefractionSP.apply();
        return true;
    }
    return false;
}
