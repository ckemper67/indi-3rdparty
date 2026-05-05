/*
    ErfaTelescopeMixin — ERFA coordinate helpers and atmospheric refraction properties.

    All methods and properties in this class are destined for INDI::Telescope once the
    refraction plan is implemented in the base class.  Delete this file and remove the
    mixin from CelestronAUXErfa's inheritance list when that promotion happens.

    Copyright (C) 2026 Christian Kemper

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.
*/

#pragma once

#include <erfa.h>
#include <erfam.h>
#include <indipropertyswitch.h>
#include <indipropertynumber.h>
#include <indipropertytext.h>
#include <inditelescope.h>

// These constants exist in SOFA but not all ERFA releases
#ifndef ERFA_DH2R
#  define ERFA_DH2R (ERFA_DPI / 12.0)
#endif
#ifndef ERFA_DR2H
#  define ERFA_DR2H (12.0 / ERFA_DPI)
#endif

class ErfaTelescopeMixin
{
    protected:
        // -- Atmospheric state --
        struct Atmosphere
        {
            double temp_C        = 15.0;
            double pressure_hPa  = 0.0;   // 0 = refraction disabled
            double humidity      = 50.0;  // percent 0–100 (divided by 100 at eraRefco call site)
            double wavelength_um = 0.55;
        } m_Atmosphere;

        bool m_RefractionEnabled = false;

        // -- INDI properties --
        INDI::PropertyNumber AtmosphericNP {4};
        enum { ATMOS_TEMPERATURE = 0, ATMOS_PRESSURE = 1, ATMOS_HUMIDITY = 2, ATMOS_WAVELENGTH = 3 };

        INDI::PropertyText   AtmosphericSourceTP {1};

        INDI::PropertySwitch RefractionSP {2};
        enum { REFRACTION_OFF = 0, REFRACTION_SOFTWARE = 1 };

        // -- Time helpers (replace libnova calls) --
        double currentJD() const;

        // Returns LST in hours; replaces get_local_sidereal_time(lon)
        double localSiderealTime(double jd, double longitude_deg) const;

        // -- ERFA context builder --
        // Builds eraASTROM at observer for given UTC JD.
        // When m_RefractionEnabled=false, passes pressure_hPa=0 so eraAtioq skips refraction.
        void buildErfaContext(double jd, eraASTROM *astrom, double *eo) const;

        // -- Coordinate conversions (replace INDI::EquatorialToHorizontal / HorizontalToEquatorial) --

        // Apparent RA/Dec (hours/deg) → observed Alt/Az (deg).
        // Refraction applied if m_RefractionEnabled.
        void erfaEqToHoriz(double ra_hours, double dec_deg, double jd,
                           double *alt_deg, double *az_deg) const;

        // Observed Alt/Az (deg) → apparent RA/Dec (hours/deg).
        // Refraction removed if m_RefractionEnabled.
        void erfaHorizToEq(double alt_deg, double az_deg, double jd,
                           double *ra_hours, double *dec_deg) const;

        // J2000 → topocentric CIRS (no refraction); for feeding alignment subsystem.
        void erfaJ2000ToTopocentric(double ra_j2000_hours, double dec_j2000_deg, double jd,
                                    double *ra_topo_hours, double *dec_topo_deg) const;

        // -- Property lifecycle --
        void initErfaProperties(const char *groupTab);
        void updateErfaProperties(bool connected);
        bool handleErfaNumber(const char *name, double *values, char **names, int n);
        bool handleErfaSwitch(const char *name, ISState *states, char **names, int n);

        // -- Observer access (implemented by the driver) --
        virtual INDI::IGeographicCoordinates getObserverLocation() const = 0;
};
