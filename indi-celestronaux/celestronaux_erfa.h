/*
    CelestronAUXErfa — clone of CelestronAUX using ERFA 2000B for all coordinate math.

    This driver is a POC for:
      - replacing libnova/libastro coordinate calls with direct ERFA 2000B calls
      - adding ATMOSPHERIC_CONDITIONS + REFRACTION_CONTROL properties
      - validating the refraction handling design before it is promoted into INDI::Telescope

    When the refraction plan lands in INDI::Telescope:
      1. ErfaTelescopeMixin contents move into INDI::Telescope
      2. Remove ErfaTelescopeMixin from this class's inheritance list
      3. Delete erfa_telescope_mixin.h / .cpp

    Copyright (C) 2020 Paweł T. Jochym
    Copyright (C) 2020 Fabrizio Pollastri
    Copyright (C) 2020-2022 Jasem Mutlaq
    Copyright (C) 2024 Christian Kemper (ERFA migration)

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.
*/

#pragma once

#include <indicom.h>
#include <indiguiderinterface.h>
#include <indifocuserinterface.h>
#include <inditelescope.h>
#include <indielapsedtimer.h>
#include <inditimer.h>
#include <connectionplugins/connectionserial.h>
#include <connectionplugins/connectiontcp.h>
#include <alignment/AlignmentSubsystemForDrivers.h>
#include <indipropertyswitch.h>
#include <indipropertynumber.h>
#include <indipropertytext.h>
#include <pid.h>
#include <termios.h>

#include "auxproto.h"
#include "adaptive_tuner.h"
#include "erfa_telescope_mixin.h"

class CelestronAUXErfa :
    public INDI::Telescope,
    public INDI::GuiderInterface,
    public INDI::FocuserInterface,
    public INDI::AlignmentSubsystem::AlignmentSubsystemForDrivers,
    public ErfaTelescopeMixin
{
    public:
        CelestronAUXErfa();
        ~CelestronAUXErfa() override;

        virtual bool ISNewBLOB(const char *dev, const char *name, int sizes[], int blobsizes[], char *blobs[],
                               char *formats[], char *names[], int n) override;
        virtual bool ISNewNumber(const char *dev, const char *name, double values[], char *names[], int n) override;
        virtual bool ISNewSwitch(const char *dev, const char *name, ISState *states, char *names[], int n) override;
        virtual bool ISNewText(const char *dev, const char *name, char *texts[], char *names[], int n) override;
        virtual bool ISSnoopDevice(XMLEle *root) override;

        // Type defs
        enum ScopeStatus_t
        {
            IDLE,
            SLEWING_FAST,
            APPROACH,
            SLEWING_SLOW,
            SLEWING_MANUAL,
            TRACKING
        };
        ScopeStatus_t ScopeStatus;

        enum AxisStatus
        {
            STOPPED,
            SLEWING
        };

        enum AxisDirection
        {
            FORWARD,
            REVERSE
        };

        enum MountVersion
        {
            GPS_Nexstar       = 0x0001,
            SLT_Nexstar       = 0x0783,
            SE_5_4            = 0x0b83,
            SE_8_6            = 0x0c82,
            CPC_Deluxe        = 0x1189,
            Series_GT         = 0x1283,
            AVX               = 0x1485,
            Evolution_Nexstar = 0x1687,
            CGX               = 0x1788,
            Advanced_GT       = 0x0682
        };

        typedef enum
        {
            PREVIOUS_NS_MOTION_NORTH   = DIRECTION_NORTH,
            PREVIOUS_NS_MOTION_SOUTH   = DIRECTION_SOUTH,
            PREVIOUS_NS_MOTION_UNKNOWN = -1
        } PreviousNSMotion_t;
        typedef enum
        {
            PREVIOUS_WE_MOTION_WEST    = DIRECTION_WEST,
            PREVIOUS_WE_MOTION_EAST    = DIRECTION_EAST,
            PREVIOUS_WE_MOTION_UNKNOWN = -1
        } PreviousWEMotion_t;

        void syncCoordWrapPosition();

    protected:
        virtual void ISGetProperties(const char *dev) override;
        virtual bool initProperties() override;
        virtual bool updateProperties() override;
        virtual bool saveConfigItems(FILE *fp) override;
        virtual bool Handshake() override;
        virtual bool Disconnect() override;

        virtual const char *getDefaultName() override;
        INDI::IHorizontalCoordinates AltAzFromRaDec(double ra, double dec, double ts);

        virtual bool Sync(double ra, double dec) override;
        virtual bool Goto(double ra, double dec) override;
        virtual bool Abort() override;
        virtual bool Park() override;
        virtual bool UnPark() override;

        virtual IPState GuideNorth(uint32_t ms) override;
        virtual IPState GuideSouth(uint32_t ms) override;
        virtual IPState GuideEast(uint32_t ms) override;
        virtual IPState GuideWest(uint32_t ms) override;

        virtual IPState MoveRelFocuser(FocusDirection dir, uint32_t ticks) override;
        virtual IPState MoveAbsFocuser (uint32_t targetTicks) override;
        virtual bool AbortFocuser () override;

        virtual bool MoveNS(INDI_DIR_NS dir, TelescopeMotionCommand command) override;
        virtual bool MoveWE(INDI_DIR_WE dir, TelescopeMotionCommand command) override;

        virtual bool ReadScopeStatus() override;

        void guideAltAzDecomposed(double dNS, double dEW, uint32_t ms);

        virtual void TimerHit() override;
        virtual bool updateLocation(double latitude, double longitude, double elevation) override;

        bool SetCurrentPark() override;
        bool SetDefaultPark() override;

        // ErfaTelescopeMixin observer access
        INDI::IGeographicCoordinates getObserverLocation() const override
        {
            return m_Location;
        }

        /////////////////////////////////////////////////////////////////////////////////////
        /// Motion Control
        /////////////////////////////////////////////////////////////////////////////////////
        bool stopAxis(INDI_HO_AXIS axis);
        bool isSlewing();

        bool slewTo(INDI_HO_AXIS axis, uint32_t steps, bool fast = true);
        bool slewByRate(INDI_HO_AXIS axis, int8_t rate);

        bool goHome(INDI_HO_AXIS axis);
        bool isHomingDone(INDI_HO_AXIS axis);
        bool m_HomingProgress[2] = {false, false};

        bool enforceSlewLimits();

        /////////////////////////////////////////////////////////////////////////////////////
        /// Tracking
        /////////////////////////////////////////////////////////////////////////////////////
        bool SetTrackEnabled(bool enabled) override;
        bool SetTrackMode(uint8_t mode) override;
        bool SetTrackRate(double raRate, double deRate) override;
        void resetTracking();

        bool trackByRate(INDI_HO_AXIS axis, int32_t rate);
        bool trackByMode(INDI_HO_AXIS axis, uint8_t mode);
        bool isTrackingRequested();

        bool getStatus(INDI_HO_AXIS axis);
        bool getEncoder(INDI_HO_AXIS axis);

        /////////////////////////////////////////////////////////////////////////////////////
        /// Coord Wrap
        /////////////////////////////////////////////////////////////////////////////////////
        bool setCordWrapEnabled(bool enable);
        bool getCordWrapEnabled();
        bool setCordWrapPosition(uint32_t steps);
        uint32_t getCordWrapPosition();

        /////////////////////////////////////////////////////////////////////////////////////
        /// Focus
        /////////////////////////////////////////////////////////////////////////////////////
        bool getFocusLimits();
        bool getFocusPosition();
        bool getFocusStatus();
        bool focusTo(uint32_t steps);
        bool focusByRate(int8_t rate);

    private:
        /////////////////////////////////////////////////////////////////////////////////////
        /// Misc
        /////////////////////////////////////////////////////////////////////////////////////
        double getNorthAz();
        bool isNorthHemisphere() const
        {
            return m_Location.latitude >= 0;
        }
        bool startupWithoutHC();
        bool getModel(AUXTargets target);
        bool getVersion(AUXTargets target);
        void getVersions();
        void hex_dump(char *buf, AUXBuffer data, size_t size);

        double AzimuthToDegrees(double degree);
        double DegreesToAzimuth(double degree);

        double EncodersToDegrees(uint32_t steps);
        uint32_t DegreesToEncoders(double degrees);

        double EncodersToHours(uint32_t steps);
        uint32_t HoursToEncoders(double hour);

        double EncodersToDE(uint32_t steps, TelescopePierSide pierSide);
        double DEToEncoders(double de);

        void EncodersToAltAz(INDI::IHorizontalCoordinates &coords);
        void EncodersToRADE(INDI::IEquatorialCoordinates &coords, TelescopePierSide &pierSide);
        void RADEToEncoders(const INDI::IEquatorialCoordinates &coords, uint32_t &haEncoder, uint32_t &deEncoder);

        bool mountToSkyCoords();

        /////////////////////////////////////////////////////////////////////////////////////
        /// Guiding
        /////////////////////////////////////////////////////////////////////////////////////
        bool guidePulse(INDI_EQ_AXIS axis, uint32_t ms, int8_t rate);
        bool guidePulsePhysical(INDI_EQ_AXIS axis, uint32_t ms, int8_t rate);
        bool getGuideRate(AUXTargets target);
        bool setGuideRate(AUXTargets target, uint8_t rate);

    private:
        AxisStatus m_AxisStatus[2] {STOPPED, STOPPED};
        AxisDirection m_AxisDirection[2] {FORWARD, FORWARD};

        double m_TrackRates[2] = {TRACKRATE_SIDEREAL, 0};

        TelescopePierSide m_TargetPierSide {PIER_UNKNOWN};

        INDI::IEquatorialCoordinates m_SkyTrackingTarget { 0, 0 };
        INDI::IEquatorialCoordinates m_SkyGOTOTarget { 0, 0 };
        INDI::IEquatorialCoordinates m_SkyCurrentRADE {0, 0};

        INDI::PropertyNumber EqPENP {2};
        enum { PE_RA = 0, PE_DEC = 1 };

        INDI::IEquatorialCoordinates m_MountCurrentRADE {0, 0};
        INDI::IHorizontalCoordinates m_MountCurrentAltAz {0, 0};

        INDI::ElapsedTimer m_TrackingElapsedTimer;
        INDI::Timer m_GuideRATimer, m_GuideDETimer;

        /////////////////////////////////////////////////////////////////////////////////////
        /// Auxiliary Command Communication
        /////////////////////////////////////////////////////////////////////////////////////
        bool sendAUXCommand(AUXCommand &command);
        void closeConnection();
        void emulateGPS(AUXCommand &m);
        bool serialReadResponse(AUXCommand c);
        bool tcpReadResponse();
        bool readAUXResponse(AUXCommand c);
        bool processResponse(AUXCommand &cmd);
        int sendBuffer(AUXBuffer buf);
        void formatModelString(char *s, int n, uint16_t model);
        void formatVersionString(char *s, int n, uint8_t *verBuf);

        bool m_GPSEmulation {false};

        uint16_t m_ModelVersion {0};
        uint8_t m_MainBoardVersion[4] {0};
        uint8_t m_AltitudeVersion[4] {0};
        uint8_t m_AzimuthVersion[4] {0};
        uint8_t m_HCVersion[4] {0};
        uint8_t m_BATVersion[4] {0};
        uint8_t m_WiFiVersion[4] {0};
        uint8_t m_GPSVersion[4] {0};
        uint8_t m_FocusVersion[4] {0};

        bool m_CordWrapActive {false};
        int32_t m_CordWrapPosition {0};
        uint32_t m_RequestedCordwrapPos;

        bool m_FocusEnabled {false};
        uint32_t m_FocusPosition {0};
        uint32_t m_FocusLimitMax {0};
        uint32_t m_FocusLimitMin {0xffffffff};
        AxisStatus m_FocusStatus {STOPPED};

        bool m_ManualMotionActive { false };

        uint32_t DBG_CAUX {0};
        uint32_t DBG_SERIAL {0};

        ///////////////////////////////////////////////////////////////////////////////
        /// Communication
        ///////////////////////////////////////////////////////////////////////////////
        int m_ModemControl {0};
        void setRTS(bool rts);
        bool waitCTS(float timeout);
        bool detectRTSCTS();
        bool detectHC(char *version, size_t size);
        int response_data_size;
        int aux_tty_read(char *buf, int bufsiz, int timeout, int *n);
        int aux_tty_write (char *buf, int bufsiz, float timeout, int *n);
        bool tty_set_speed(speed_t speed);

        bool m_IsRTSCTS {false};
        bool m_isHandController {false};

        ///////////////////////////////////////////////////////////////////////////////
        /// Celestron AUX Properties
        ///////////////////////////////////////////////////////////////////////////////
        INDI::PropertyText FirmwareTP {9};
        enum {FW_MODEL, FW_HC, FW_MB, FW_AZM, FW_ALT, FW_WiFi, FW_BAT, FW_GPS, FW_FOCUS};

        INDI::PropertySwitch CordWrapToggleSP {2};

        INDI::PropertySwitch CordWrapPositionSP {4};
        enum { CORDWRAP_N, CORDWRAP_E, CORDWRAP_S, CORDWRAP_W };

        INDI::PropertySwitch CordWrapBaseSP {2};
        enum {CW_BASE_ENC, CW_BASE_SKY};

        INDI::PropertySwitch Axis1LimitToggleSP {2};
        INDI::PropertySwitch Axis2LimitToggleSP {2};
        INDI::PropertyNumber SlewLimitPositionNP {4};
        enum { SLEW_LIMIT_AXIS1_MIN, SLEW_LIMIT_AXIS1_MAX, SLEW_LIMIT_AXIS2_MIN, SLEW_LIMIT_AXIS2_MAX };

        INDI::PropertySwitch ApproachDirectionSP {2};
        enum { APPROACH_TRACKING_VECTOR, APPROACH_CONSTANT_OFFSET };

        INDI::PropertySwitch GPSEmuSP {2};
        enum { GPSEMU_OFF, GPSEMU_ON };

        INDI::PropertyNumber HorizontalCoordsNP {2};

        INDI::PropertyNumber GuideRateNP {2};

        INDI::PropertyNumber EncoderNP {2};
        INDI::PropertyNumber AngleNP {2};

        int32_t m_LastTrackRate[2] = {-1, -1};
        double m_TrackStartSteps[2] = {0, 0};
        int32_t m_LastOffset[2] = {0, 0};
        int m_OffsetSwitchSettle[2] = {0, 0};

        INDI::PropertyNumber Axis1PIDNP {3};
        INDI::PropertyNumber Axis2PIDNP {3};

        bool m_IsPipelinePrimed { false };
        INDI::IHorizontalCoordinates m_TrackingWindowCoords[3];
        INDI::IEquatorialCoordinates m_LastTrackingTarget { 0, 0 };
        double m_LastTrackingDt { 0 };

        // Ephemeris tracking state (solar/lunar): last sampled JNow RA/Dec
        double m_lastEphemRA  { 0 };
        double m_lastEphemDec { 0 };
        bool   m_ephemPrimed  { false };

        enum
        {
            Propotional,
            Derivative,
            Integral
        };

        std::unique_ptr<PID> m_Controllers[2];
        std::unique_ptr<AdaptivePIDTuner> m_az_pid_tuner;
        std::unique_ptr<AdaptivePIDTuner> m_al_pid_tuner;

        INDI::PropertySwitch PortTypeSP {2};
        enum
        {
            PORT_AUX_PC,
            PORT_HC_USB,
        };

        int m_ConfigPortType {PORT_AUX_PC};

        INDI::PropertySwitch HomeSP {3};
        enum
        {
            HOME_AXIS1,
            HOME_AXIS2,
            HOME_ALL
        };

        INDI::PropertySwitch AdaptiveTuningAzSP {2};
        INDI::PropertySwitch AdaptiveTuningAlSP {2};
        INDI::PropertyNumber UpdateRateNP {1};

        typedef enum { ALT_AZ, EQ_FORK, EQ_GEM } MountType;
        MountType m_MountType {ALT_AZ};

    private:
        static constexpr int32_t STEPS_PER_REVOLUTION {16777216};
        static constexpr double STEPS_PER_DEGREE {STEPS_PER_REVOLUTION / 360.0};
        static constexpr double STEPS_PER_ARCSEC {STEPS_PER_DEGREE / 3600.0};
        static constexpr double DEGREES_PER_STEP {360.0 / STEPS_PER_REVOLUTION};

        static constexpr double STEPS_PER_HOUR {STEPS_PER_REVOLUTION / 24.0};
        static constexpr double HOURS_PER_STEP {24.0 / STEPS_PER_REVOLUTION};

        static constexpr uint8_t RATE_PER_ARCSEC {4};
        static constexpr uint32_t BUFFER_SIZE {10240};
        static constexpr uint8_t READ_TIMEOUT {1};
        static constexpr uint8_t CTS_TIMEOUT {100};
        static constexpr const char *CORDWRAP_TAB {"Coord Wrap"};
        static constexpr const char *MOUNTINFO_TAB {"Mount Info"};
        static constexpr uint16_t AUX_SIDEREAL {0xffff};
        static constexpr uint16_t AUX_SOLAR {0xfffe};
        static constexpr uint16_t AUX_LUNAR {0xfffd};
        static constexpr uint32_t GEM_HOME {4194304};
};
