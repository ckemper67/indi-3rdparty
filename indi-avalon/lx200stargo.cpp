/*
    Avalon StarGo driver

    Copyright (C) 2018 Christopher Contaxis, Wolfgang Reissenberger,
    Tonino Tasselli and Ken Self

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


#include "lx200stargo.h"

#include "lx200stargofocuser.h"

#include <cmath>
#include <memory>
#include <cstring>
#include <unistd.h>
#ifndef _WIN32
#include <termios.h>
#endif
//#if defined(HAVE_LIBNOVA)
#include <libnova/julian_day.h>
#include <libnova/sidereal_time.h>
//#endif // HAVE_LIBNOVA

// We declare an auto pointer to LX200StarGo
std::unique_ptr<LX200StarGo> telescope;
std::unique_ptr<LX200StarGoFocuser> focuser;

const char *RA_DEC_TAB = "RA / DEC";

void ISInit()
{
    static int isInit = 0;

    if (isInit)
        return;

    isInit = 1;
    if (telescope.get() == 0) 
    {
        LX200StarGo* myScope = new LX200StarGo();
        telescope.reset(myScope);
        focuser.reset(new LX200StarGoFocuser(myScope, "AUX1 Focuser"));
    }
}

void ISGetProperties(const char *dev)
{
    ISInit();
    telescope->ISGetProperties(dev);
//    focuser->ISGetProperties(dev);
}

void ISNewSwitch(const char *dev, const char *name, ISState *states, char *names[], int n)
{
    ISInit();
    telescope->ISNewSwitch(dev, name, states, names, n);
    focuser->ISNewSwitch(dev, name, states, names, n);
}

void ISNewText(const char *dev, const char *name, char *texts[], char *names[], int n)
{
    ISInit();
    telescope->ISNewText(dev, name, texts, names, n);
//    focuser->ISNewText(dev, name, texts, names, n);
}

void ISNewNumber(const char *dev, const char *name, double values[], char *names[], int n)
{
    ISInit();
    telescope->ISNewNumber(dev, name, values, names, n);
    focuser->ISNewNumber(dev, name, values, names, n);
}

void ISNewBLOB(const char *dev, const char *name, int sizes[], int blobsizes[], char *blobs[], char *formats[],
               char *names[], int n)
{
    INDI_UNUSED(dev);
    INDI_UNUSED(name);
    INDI_UNUSED(sizes);
    INDI_UNUSED(blobsizes);
    INDI_UNUSED(blobs);
    INDI_UNUSED(formats);
    INDI_UNUSED(names);
    INDI_UNUSED(n);
}
void ISSnoopDevice(XMLEle *root)
{
    ISInit();
    telescope->ISSnoopDevice(root);
//    focuser->ISSnoopDevice(root);
}

/**************************************************
*** LX200 Generic Implementation
***************************************************/


LX200StarGo::LX200StarGo()
{
    setVersion(0, 5);
    
    DBG_SCOPE = INDI::Logger::DBG_DEBUG;
 
    /* missing capabilities
     * TELESCOPE_HAS_TIME:
     *    missing commands
     *      :GG# (Get UTC offset time)
     *      :GL# (Get Local Time in 24 hour format)
     *
     * TELESCOPE_HAS_LOCATION --Now Implemented
     * reading the location works, setting location not
     *     missing commands
     *       :SgDDD*MM# (Set current site’s longitude)
     *       :StsDD*MM# (Sets the current site latitude)
     *
     * LX200_HAS_ALIGNMENT_TYPE
     *     missing commands
     *        ACK - Alignment Query
     *
     * LX200_HAS_SITES
     *    Makes no sense in combination with KStars?
     *     missing commands
     *        :GM# (Get Site 1 Name)
     *
     * LX200_HAS_TRACKING_FREQ
     *     missing commands
     *        :GT# (Get tracking rate)
     *
     * untested, hence disabled:
     * LX200_HAS_FOCUS
     */

    setLX200Capability(LX200_HAS_PULSE_GUIDING);

    SetTelescopeCapability(TELESCOPE_CAN_PARK | TELESCOPE_CAN_SYNC | TELESCOPE_CAN_GOTO | TELESCOPE_CAN_ABORT |
                           TELESCOPE_HAS_TRACK_MODE | TELESCOPE_HAS_LOCATION | TELESCOPE_CAN_CONTROL_TRACK | TELESCOPE_HAS_PIER_SIDE, 4);
}

/**************************************************************************************
**
***************************************************************************************/
const char *LX200StarGo::getDefaultName()
{
    return (const char *)"Avalon StarGo";
}


/**************************************************************************************
**
***************************************************************************************/
bool LX200StarGo::Handshake()
{
    LOG_DEBUG(__FUNCTION__);
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if(!sendQuery(":GR#", response))
    {
        LOGF_ERROR("%s Error communication with telescope.", getDeviceName());
        return false;
    }
    if (f_scansexa(response, &currentRA))
    {
        LOGF_ERROR("%s: Unable to parse RA response %s", getDeviceName(), response);
        return false;
    }

    return true;
}


/**************************************************************************************
**
***************************************************************************************/
bool LX200StarGo::ISNewSwitch(const char *dev, const char *name, ISState *states, char *names[], int n)
{

    if (dev != nullptr && strcmp(dev, getDeviceName()) == 0)
    {

        // sync home position
        if (!strcmp(name, SyncHomeSP.name))
        {
            return syncHomePosition();
        }

        // goto home position
        if (!strcmp(name, MountGotoHomeSP.name))
        {
            return slewToHome(states, names, n);
        }
        // parking position
        else if (!strcmp(name, MountSetParkSP.name))
        {
            return setParkPosition(states, names, n);
        }
        // tracking mode
        else if (!strcmp(name, TrackModeSP.name))
        {
            if (IUUpdateSwitch(&TrackModeSP, states, names, n) < 0)
                return false;
            int trackMode = IUFindOnSwitchIndex(&TrackModeSP);

            bool result;
            if (trackMode != TRACK_NONE) result = SetTrackMode(trackMode);

            switch (trackMode) {
            case TRACK_SIDEREAL:
                LOG_INFO("Sidereal tracking rate selected.");
                break;
            case TRACK_SOLAR:
                LOG_INFO("Solar tracking rate selected.");
                break;
            case TRACK_LUNAR:
                LOG_INFO("Lunar tracking rate selected");
                break;
            case TRACK_NONE:
                LOG_INFO("Tracking stopped.");
                result = querySetTracking(false);
                break;
            }
            TrackModeSP.s = result ? IPS_OK : IPS_ALERT;

            IDSetSwitch(&TrackModeSP, nullptr);
            return result;
        } 
/*        else if (!strcmp(name, TrackStateSP.name)) 
        {
            LOG_ERROR("Set Track State."); 
            return true;
        }
 * */
        else if (!strcmp(name, ST4StatusSP.name)) 
        {
            bool enabled = (states[0] == ISS_OFF);
            bool result = setST4Enabled(enabled);

            if(result) {
                ST4StatusS[0].s = enabled ? ISS_OFF : ISS_ON;
                ST4StatusS[1].s = enabled ? ISS_ON : ISS_OFF;
                ST4StatusSP.s = IPS_OK;
            } else {
                ST4StatusSP.s = IPS_ALERT;
            }
            IDSetSwitch(&ST4StatusSP, nullptr);
            return result;
        } 
        else if (!strcmp(name, MeridianFlipEnabledSP.name)) 
        {
            return setMeridianFlipEnabled(states[0] == ISS_ON);
        } 
        else if (!strcmp(name, MeridianFlipForcedSP.name)) 
        {
            return setMeridianFlipForced(states[0] == ISS_OFF);
        }

    }

    //  Nobody has claimed this, so pass it to the parent
    return LX200Telescope::ISNewSwitch(dev, name, states, names, n);
}

bool LX200StarGo::ISNewNumber(const char *dev, const char *name, double values[], char *names[], int n) 
{

    if (dev != nullptr && strcmp(dev, getDeviceName()) == 0) 
    {

        // sync home position
        if (!strcmp(name, GuidingSpeedNP.name)) 
        {
            int raSpeed  = round(values[0] * 100);
            int decSpeed = round(values[1] * 100);
            bool result  = setGuidingSpeeds(raSpeed, decSpeed);

            if(result) 
            {
                GuidingSpeedP[0].value = static_cast<double>(raSpeed) / 100.0;
                GuidingSpeedP[1].value = static_cast<double>(decSpeed) / 100.0;
                GuidingSpeedNP.s = IPS_OK;
            } 
            else 
            {
                GuidingSpeedNP.s = IPS_ALERT;
            }
            IDSetNumber(&GuidingSpeedNP, nullptr);
            return result;
        }
    }
        //  Nobody has claimed this, so pass it to the parent
    return LX200Telescope::ISNewNumber(dev, name, values, names, n);
}



/**************************************************************************************
**
***************************************************************************************/
bool LX200StarGo::initProperties()
{
    /* Make sure to init parent properties first */
    if (!LX200Telescope::initProperties()) return false;

    IUFillSwitch(&MountGotoHomeS[0], "MOUNT_GOTO_HOME_VALUE", "Goto Home", ISS_OFF);
    IUFillSwitchVector(&MountGotoHomeSP, MountGotoHomeS, 1, getDeviceName(), "MOUNT_GOTO_HOME", "Goto Home", MAIN_CONTROL_TAB, IP_RW, ISR_ATMOST1, 60, IPS_OK);

    IUFillLight(&MountParkingStatusL[0], "MOUNT_IS_PARKED_VALUE", "Parked", IPS_IDLE);
    IUFillLight(&MountParkingStatusL[1], "MOUNT_IS_UNPARKED_VALUE", "Unparked", IPS_IDLE);
    IUFillLightVector(&MountParkingStatusLP, MountParkingStatusL, 2, getDeviceName(), "PARKING_STATUS", "Parking Status", MAIN_CONTROL_TAB, IPS_IDLE);

    IUFillSwitch(&MountSetParkS[0], "MOUNT_SET_PARK_VALUE", "Set Park", ISS_OFF);
    IUFillSwitchVector(&MountSetParkSP, MountSetParkS, 1, getDeviceName(), "MOUNT_SET_PARK", "Set Park", MAIN_CONTROL_TAB, IP_RW, ISR_ATMOST1, 60, IPS_OK);

    IUFillSwitch(&SyncHomeS[0], "SYNC_HOME", "Sync Home", ISS_OFF);
    IUFillSwitchVector(&SyncHomeSP, SyncHomeS, 1, getDeviceName(), "TELESCOPE_SYNC_HOME", "Home Position", MAIN_CONTROL_TAB,
                       IP_RW, ISR_ATMOST1, 60, IPS_IDLE);

    IUFillText(&MountFirmwareInfoT[0], "MOUNT_FIRMWARE_INFO", "Firmware", "");
    IUFillTextVector(&MountInfoTP, MountFirmwareInfoT, 1, getDeviceName(), "MOUNT_INFO", "Mount Info", INFO_TAB, IP_RO, 60, IPS_OK);

    // Guiding settings
    IUFillNumber(&GuidingSpeedP[0], "GUIDING_SPEED_RA", "RA Speed", "%.2f", 0.0, 2.0, 0.1, 0);
    IUFillNumber(&GuidingSpeedP[1], "GUIDING_SPEED_DEC", "DEC Speed", "%.2f", 0.0, 2.0, 0.1, 0);
    IUFillNumberVector(&GuidingSpeedNP, GuidingSpeedP, 2, getDeviceName(), "GUIDING_SPEED","Autoguiding", RA_DEC_TAB, IP_RW, 60, IPS_IDLE);

    IUFillSwitch(&ST4StatusS[0], "ST4_DISABLED", "disabled", ISS_OFF);
    IUFillSwitch(&ST4StatusS[1], "ST4_ENABLED", "enabled", ISS_OFF);
    IUFillSwitchVector(&ST4StatusSP, ST4StatusS, 2, getDeviceName(), "ST4", "ST4", RA_DEC_TAB, IP_RW, ISR_ATMOST1, 60, IPS_IDLE);

    // meridian flip
    IUFillSwitch(&MeridianFlipEnabledS[0], "MERIDIAN_FLIP_ENABLED", "enabled", ISS_OFF);
    IUFillSwitch(&MeridianFlipEnabledS[1], "MERIDIAN_FLIP_DISABLED", "disabled", ISS_OFF);
    IUFillSwitchVector(&MeridianFlipEnabledSP, MeridianFlipEnabledS, 2, getDeviceName(), "ENABLE_MERIDIAN_FLIP", "Meridian Flip", RA_DEC_TAB, IP_RW, ISR_ATMOST1, 60, IPS_IDLE);


    IUFillSwitch(&MeridianFlipForcedS[0], "MERIDIAN_FLIP_AUTOMATIC", "automatic", ISS_OFF);
    IUFillSwitch(&MeridianFlipForcedS[1], "MERIDIAN_FLIP_FORCED", "forced", ISS_OFF);
    IUFillSwitchVector(&MeridianFlipForcedSP, MeridianFlipForcedS, 2, getDeviceName(), "FORCE_MERIDIAN_FLIP", "Flip Mode", RA_DEC_TAB, IP_RW, ISR_ATMOST1, 60, IPS_IDLE);

    // overwrite the custom tracking mode button
    IUFillSwitch(&TrackModeS[3], "TRACK_NONE", "None", ISS_OFF);

    focuser->initProperties("AUX1 Focuser");

    return true;
}

/**************************************************************************************
**
***************************************************************************************/
bool LX200StarGo::updateProperties()
{
    char firmwareInfo[48] = {0};

    if (! LX200Telescope::updateProperties()) return false;
    if (isConnected())
    {
        defineLight(&MountParkingStatusLP);
        defineSwitch(&SyncHomeSP);
        defineSwitch(&MountGotoHomeSP);
        defineSwitch(&MountSetParkSP);
        defineNumber(&GuidingSpeedNP);
        defineSwitch(&ST4StatusSP);
        defineSwitch(&MeridianFlipEnabledSP);
        defineSwitch(&MeridianFlipForcedSP);

        if (queryFirmwareInfo(firmwareInfo)) 
        {
            MountFirmwareInfoT[0].text = firmwareInfo;
            defineText(&MountInfoTP);
        }
        bool isParked, isSynched;
        if (queryParkSync(&isParked, &isSynched)) 
        {
            SetParked(isParked);
            if (isSynched) 
            {
                SyncHomeS[0].s = ISS_ON;
                SyncHomeSP.s = IPS_OK;
                IDSetSwitch(&SyncHomeSP, nullptr);
            }
        }
        bool isEnabled;
        if (queryGetST4Status(&isEnabled)) 
        {
            ST4StatusS[0].s = isEnabled ? ISS_OFF : ISS_ON;
            ST4StatusS[1].s = isEnabled ? ISS_ON : ISS_OFF;
            ST4StatusSP.s = IPS_OK;
        } 
        else 
        {
            ST4StatusSP.s = IPS_ALERT;
        }
        IDSetSwitch(&ST4StatusSP, nullptr);

        if (queryGetMeridianFlipEnabledStatus(&isEnabled)) 
        {
            MeridianFlipEnabledS[0].s = isEnabled ? ISS_ON  : ISS_OFF;
            MeridianFlipEnabledS[1].s = isEnabled ? ISS_OFF : ISS_ON;
            MeridianFlipEnabledSP.s = IPS_OK;
        } 
        else 
        {
            MeridianFlipEnabledSP.s = IPS_ALERT;
        }
        IDSetSwitch(&MeridianFlipEnabledSP, nullptr);

        int raSpeed, decSpeed;
        if (queryGetGuidingSpeeds(&raSpeed, &decSpeed)) 
        {
            GuidingSpeedP[0].value = static_cast<double>(raSpeed) / 100.0;
            GuidingSpeedP[1].value = static_cast<double>(decSpeed) / 100.0;
            GuidingSpeedNP.s = IPS_OK;
        } 
        else 
        {
            GuidingSpeedNP.s = IPS_ALERT;
        }
        IDSetNumber(&GuidingSpeedNP, nullptr);
    }
    else
    {
        deleteProperty(MountParkingStatusLP.name);
        deleteProperty(MountGotoHomeSP.name);
        deleteProperty(MountSetParkSP.name);
        deleteProperty(SyncHomeSP.name);
        deleteProperty(MountInfoTP.name);
        deleteProperty(GuidingSpeedNP.name);
        deleteProperty(ST4StatusSP.name);
        deleteProperty(MeridianFlipEnabledSP.name);
        deleteProperty(MeridianFlipForcedSP.name);
    }

    focuser->updateProperties();

    return true;
}

/**************************************************************************************
**
***************************************************************************************/
bool LX200StarGo::ReadScopeStatus()
{
    LOG_DEBUG("ReadScopeStatus");
    if (!isConnected())
        return false;

    if (isSimulation())
    {
        mountSim();
        return true;
    }

    if (TrackState == SCOPE_SLEWING)
    {
        // Check if LX200 is done slewing
        if (isSlewComplete()) 
        {
            if (isIdle()) 
            {
                TrackState = SCOPE_IDLE;
                LOG_INFO("Slew is complete. Tracking is off." );
            }  
            else 
            {
                TrackState = SCOPE_TRACKING;
                LOG_INFO("Slew is complete. Tracking...");
            }

            if (MountGotoHomeSP.s == IPS_BUSY) 
            {
                MountGotoHomeSP.s = IPS_OK;
                IDSetSwitch(&MountGotoHomeSP, nullptr);
            }
        }
    }
    else if (TrackState == SCOPE_PARKING)
    {
        if (isSlewComplete()) SetParked(true);
    }

    // update meridian flip status
    bool isEnabled;
    if (queryGetMeridianFlipForcedStatus(&isEnabled)) 
    {
        MeridianFlipForcedS[0].s = isEnabled ? ISS_OFF : ISS_ON;
        MeridianFlipForcedS[1].s = isEnabled ? ISS_ON : ISS_OFF;
        MeridianFlipForcedSP.s = IPS_OK;
    } 
    else 
    {
        MeridianFlipForcedSP.s = IPS_ALERT;
    }
    IDSetSwitch(&MeridianFlipForcedSP, nullptr);

    UpdateMotionStatus();
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if(!sendQuery(":GR#", response))
    {
        LOGF_ERROR("%s: Unable to get RA", getDeviceName());        
    }
    else if (f_scansexa(response, &currentRA))
    {
        LOGF_ERROR("%s: Unable to parse response %s.", getDeviceName(), response);
        EqNP.s = IPS_ALERT;
        IDSetNumber(&EqNP, "Error reading RA.");
    }
    if(!sendQuery(":GD#", response))
    {
        LOGF_ERROR("%s: Unable to get DEC", getDeviceName());                
    }
    else if (f_scansexa(response, &currentDEC))
    {
        LOGF_ERROR("%s: Unable to parse response %s.", getDeviceName(), response);
        EqNP.s = IPS_ALERT;
        IDSetNumber(&EqNP, "Error reading DEC.");
    }

    NewRaDec(currentRA, currentDEC);
    return syncSideOfPier();
}

/**************************************************************************************
**
***************************************************************************************/

bool LX200StarGo::isIdle() 
{
    LOG_DEBUG("isIdle");
    queryMountMotionState();

    return (CurrentTrackMode == TRACK_SIDEREAL);
}

bool LX200StarGo::isSlewComplete()
{
    return (queryIsSlewComplete());
}

/**************************************************************************************
**
***************************************************************************************/

bool LX200StarGo::UpdateMotionStatus() 
{
    LOG_DEBUG("UpdateMotionStatus");
    if(!queryMountMotionState())
        return false;

    switch (CurrentMotorsState) 
    {
    case MOTORS_OFF:
        if (TrackState == SCOPE_PARKING) 
        {
            SetParked(true);
        }
        break;
    default:
        if (TrackState != SCOPE_PARKING && TrackState != SCOPE_SLEWING) 
        {
            TrackState = SCOPE_TRACKING;
        }
        break;
    }
    switch (CurrentTrackMode) 
    {
    case TRACK_NONE:
        if (TrackModeS[TRACK_CUSTOM].s != ISS_ON) 
        {
            IUResetSwitch(&TrackModeSP);
            TrackModeS[TRACK_CUSTOM].s = ISS_ON;
            TrackModeSP.s   = IPS_OK;
            IDSetSwitch(&TrackModeSP, nullptr);
        }
        if (TrackState != SCOPE_PARKED && TrackState != SCOPE_PARKING && TrackState != SCOPE_SLEWING) 
        {
            TrackState = SCOPE_IDLE;
        }
        break;
    case TRACK_LUNAR:
        if (TrackModeS[TRACK_LUNAR].s != ISS_ON) 
        {
            IUResetSwitch(&TrackModeSP);
            TrackModeS[TRACK_LUNAR].s = ISS_ON;
            TrackModeSP.s   = IPS_OK;
            IDSetSwitch(&TrackModeSP, nullptr);
        }
        if (TrackState != SCOPE_PARKED && TrackState != SCOPE_PARKING && TrackState != SCOPE_SLEWING) 
        {
            TrackState = (CurrentMotorsState == MOTORS_OFF) ? SCOPE_IDLE : SCOPE_TRACKING;
        }
        break;
    case TRACK_SOLAR:
        if (TrackModeS[TRACK_SOLAR].s != ISS_ON) 
        {
            IUResetSwitch(&TrackModeSP);
            TrackModeS[TRACK_SOLAR].s = ISS_ON;
            TrackModeSP.s   = IPS_OK;
            IDSetSwitch(&TrackModeSP, nullptr);
        }
        if (TrackState != SCOPE_PARKED && TrackState != SCOPE_PARKING && TrackState != SCOPE_SLEWING) 
        {
            TrackState = (CurrentMotorsState == MOTORS_OFF) ? SCOPE_IDLE : SCOPE_TRACKING;
        }
        break;
    case TRACK_SIDEREAL:
        if (TrackModeS[TRACK_SIDEREAL].s != ISS_ON) 
        {
            IUResetSwitch(&TrackModeSP);
            TrackModeS[TRACK_SIDEREAL].s = ISS_ON;
            TrackModeSP.s   = IPS_OK;
            IDSetSwitch(&TrackModeSP, nullptr);
        }
        if (TrackState != SCOPE_PARKED && TrackState != SCOPE_PARKING && TrackState != SCOPE_SLEWING) 
        {
            TrackState = (CurrentMotorsState == MOTORS_OFF) ? SCOPE_IDLE : SCOPE_TRACKING;
        }
        break;
    }

    switch (CurrentSlewRate) 
    {
    case SLEW_GUIDE:
        if (CurrentSlewRate != SLEW_GUIDE) 
        {
            IUResetSwitch(&SlewRateSP);
            CurrentSlewRate = SLEW_GUIDE;
            SlewRateS[SLEW_GUIDE].s = ISS_ON;
            SlewRateSP.s   = IPS_OK;
            IDSetSwitch(&SlewRateSP, nullptr);
        }
        break;
    case SLEW_CENTERING:
        if (CurrentSlewRate != SLEW_CENTERING) 
        {
            IUResetSwitch(&SlewRateSP);
            CurrentSlewRate = SLEW_CENTERING;
            SlewRateS[SLEW_CENTERING].s = ISS_ON;
            SlewRateSP.s   = IPS_OK;
            IDSetSwitch(&SlewRateSP, nullptr);
        }
        break;
    case SLEW_FIND:
        if (CurrentSlewRate != SLEW_FIND) 
        {
            IUResetSwitch(&SlewRateSP);
            CurrentSlewRate = SLEW_FIND;
            SlewRateS[SLEW_FIND].s = ISS_ON;
            SlewRateSP.s   = IPS_OK;
            IDSetSwitch(&SlewRateSP, nullptr);
        }
        break;
    case SLEW_MAX:
        if (CurrentSlewRate != SLEW_MAX) 
        {
            IUResetSwitch(&SlewRateSP);
            CurrentSlewRate = SLEW_MAX;
            SlewRateS[SLEW_MAX].s = ISS_ON;
            SlewRateSP.s   = IPS_OK;
            IDSetSwitch(&SlewRateSP, nullptr);
        }
        break;
    }
    return true;
}

/**************************************************************************************
**
***************************************************************************************/

bool LX200StarGo::syncHomePosition()
{
    LOG_DEBUG("syncHomePosition");
    char input[AVALON_RESPONSE_BUFFER_LENGTH];
    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    if (!getLST_String(input)) 
    {
        LOGF_WARN("%s Synching home get LST failed.", getDeviceName());
        SyncHomeSP.s = IPS_ALERT;
        return false;
    }

    sprintf(cmd, ":X31%s#", input);
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    
    if (sendQuery(cmd, response))
    {
        LOG_INFO("Synching home position succeeded.");
        SyncHomeSP.s = IPS_OK;
    } 
    else
    {
        LOGF_WARN("%s Synching home position failed.", getDeviceName());
        SyncHomeSP.s = IPS_ALERT;
        return false;
    }
    IDSetSwitch(&SyncHomeSP, nullptr);
    return true; //(result == TTY_OK);
}

/**************************************************************************************
* @author CanisUrsa
***************************************************************************************/

bool LX200StarGo::slewToHome(ISState* states, char* names[], int n) 
{
    LOG_DEBUG("slewToHome");
    IUUpdateSwitch(&MountGotoHomeSP, states, names, n);
    if (querySendMountGotoHome()) 
    {
        MountGotoHomeSP.s = IPS_BUSY;
        TrackState = SCOPE_SLEWING;
    } 
    else 
    {
        MountGotoHomeSP.s = IPS_ALERT;
    }
    MountGotoHomeS[0].s = ISS_OFF;
    IDSetSwitch(&MountGotoHomeSP, nullptr);

    LOG_INFO("Slewing to home position...");
    return true;
}


bool LX200StarGo::setParkPosition(ISState* states, char* names[], int n) 
{
    LOG_DEBUG("setParkPosition");
    IUUpdateSwitch(&MountSetParkSP, states, names, n);
    MountSetParkSP.s = querySendMountSetPark() ? IPS_OK : IPS_ALERT;
    MountSetParkS[0].s = ISS_OFF;
    IDSetSwitch(&MountSetParkSP, nullptr);
    return true;
}

/**************************************************************************************
**
***************************************************************************************/

void LX200StarGo::getBasicData()
{
    LOG_DEBUG("getBasicData");
    if (!isSimulation())
    {
        checkLX200Format();

        if (genericCapability & LX200_HAS_ALIGNMENT_TYPE)
            getAlignment();

        if (genericCapability & LX200_HAS_TRACKING_FREQ)
        {
            if (getTrackFrequency(&TrackFreqN[0].value) < 0)
                LOGF_ERROR("%s Failed to get tracking frequency from device.", getDeviceName());
            else
                IDSetNumber(&TrackingFreqNP, nullptr);
        }
    }

    if (sendLocationOnStartup && (GetTelescopeCapability() & TELESCOPE_HAS_LOCATION))
        sendScopeLocation();
    if (sendTimeOnStartup && (GetTelescopeCapability() & TELESCOPE_HAS_TIME))
        sendScopeTime();
}

/**************************************************************************************
* @author CanisUrsa
***************************************************************************************/
bool LX200StarGo::querySendMountGotoHome() 
{
    // Command  - :X361#
    // Response - pA#
    //            :Z1303#
    //            p0#
    //            :Z1003#
    //            p0#
    LOG_DEBUG("querySendMountGotoHome");
    char response[AVALON_COMMAND_BUFFER_LENGTH] = {0};
    if (!sendQuery(":X361#", response)) 
    {
        LOGF_ERROR("%s: Failed to send mount goto home command.", getDeviceName());
        return false;
    }
    if (strcmp(response, "pA") != 0) 
    {
        LOGF_ERROR("%s: Invalid mount sync goto response '%s'.", getDeviceName(), response);
        return false;
    }
    return true;
}

/**************************************************************************************
**
***************************************************************************************/
bool LX200StarGo::sendScopeLocation()
{
    LOG_DEBUG("sendScopeLocation");
    if (isSimulation())
    {
        return LX200Telescope::sendScopeLocation();
    }

    double siteLat = 0.0, siteLong = 0.0;
    if (!getSiteLatitude(&siteLat))
    {
        LOGF_WARN("%s Failed to get site latitude from device.", getDeviceName());
        return false;
    }
    if (!getSiteLongitude(&siteLong))
    {
        LOGF_WARN("%s Failed to get site longitude from device.", getDeviceName());
        return false;
    }
    LocationNP.np[LOCATION_LATITUDE].value = siteLat;
    LocationNP.np[LOCATION_LONGITUDE].value = siteLong;

    LOGF_DEBUG("Mount Controller Latitude: %g Longitude: %g", LocationN[LOCATION_LATITUDE].value, LocationN[LOCATION_LONGITUDE].value);

    IDSetNumber(&LocationNP, nullptr);

    return true;
}

/**************************************************************************************
**
***************************************************************************************/

bool LX200StarGo::updateLocation(double latitude, double longitude, double elevation)
{
    LOG_DEBUG("updateLocation");
    INDI_UNUSED(elevation);

    if (isSimulation())
        return true;

    LOGF_DEBUG("Setting site longitude '%lf'", longitude);
    if (!isSimulation() && ! querySetSiteLongitude(longitude))
    {
        LOGF_ERROR("%s Error setting site longitude coordinates", getDeviceName());
        return false;
    }

    if (!isSimulation() && ! querySetSiteLatitude(latitude))
    {
        LOGF_ERROR("%s Error setting site latitude coordinates", getDeviceName());
        return false;
    }

    char l[32]={0}, L[32]={0};
    fs_sexa(l, latitude, 3, 3600);
    fs_sexa(L, longitude, 4, 3600);

    LOGF_INFO("Site location updated to Lat %.32s - Long %.32s", l, L);

    return true;
}

/*
 * Determine the site latitude. In contrast to a standard LX200 implementation,
 * StarGo returns the location in arc seconds precision.
 */

bool LX200StarGo::getSiteLatitude(double *siteLat) 
{
    LOG_DEBUG("getSitelatitude");
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":Gt#", response)) {
        LOGF_ERROR("%s: Failed to send query get Site Latitude command.", getDeviceName());
        return false;
    }
    if (f_scansexa(response, siteLat))
    {
        LOGF_ERROR("%s: Unable to parse get Site Latitude response.", getDeviceName());
        return false;
    }
    return true;
}

/*
 * Determine the site longitude. In contrast to a standard LX200 implementation,
 * StarGo returns the location in arc seconds precision.
 */

bool LX200StarGo::getSiteLongitude(double *siteLong) 
{
    LOG_DEBUG("getSiteLongitude");
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":Gg#", response)) 
    {
        LOGF_ERROR("%s: Failed to send query get Site Longitude command.", getDeviceName());
        return false;
    }
    if (f_scansexa(response, siteLong))
    {
        LOGF_ERROR("%s: Unable to parse get Site Longitude response.", getDeviceName());
        return false;
    }
    return true;
}


/**************************************************************************************
**
***************************************************************************************/

bool LX200StarGo::Park() 
{
    LOG_DEBUG("Park");
    // in: :X362#
    // out: "pB#"

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (sendQuery(":X362#", response) && strcmp(response, "pB") == 0) 
    {
        LOG_INFO("Parking scope...");
        // update state
        TrackState = SCOPE_PARKING;
        return true;
    } 
    else 
    {
        LOGF_ERROR("%s: Parking failed.", getDeviceName());
        return false;
    }
}

/**
 * @brief Set parking state to "parked" and reflect the state
 *        in the UI.
 * @param isparked true iff the scope has been parked
 * @return
 */
void LX200StarGo::SetParked(bool isparked) 
{
    LOG_DEBUG("SetParked");
    INDI::Telescope::SetParked(isparked);

    TrackState = isparked ? SCOPE_PARKED : SCOPE_TRACKING;
    ParkS[0].s = isparked ? ISS_ON : ISS_OFF;
    ParkS[1].s = isparked ? ISS_OFF : ISS_ON;
    ParkSP.s   = IPS_OK;
    IDSetSwitch(&ParkSP, nullptr);
    MountParkingStatusL[0].s = isparked ? IPS_OK : IPS_IDLE;
    MountParkingStatusL[1].s = isparked ? IPS_IDLE : IPS_OK;
    IDSetLight(&MountParkingStatusLP, nullptr);
}

bool LX200StarGo::UnPark() 
{
    LOG_DEBUG("UnPark");
    // in: :X370#
    // out: "p0#"

    // set LST to avoid errors
    char input[AVALON_RESPONSE_BUFFER_LENGTH]; 
    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    if (!getLST_String(input))
    {
        LOGF_WARN("%s Obtaining LST before unparking failed.", getDeviceName());
        return false;
    }
    sprintf(cmd, ":X32%s#", input);

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    bool result = sendQuery(cmd, response);

    if (!result || response[0] != '0')
    {
        LOGF_WARN("%s Setting LST before unparking failed.", getDeviceName());
        return false;
    }

    // and now execute unparking
    response[0] = '\0';
    if (sendQuery(":X370#", response) && strcmp(response, "p0") == 0) 
    {
        LOG_INFO("Scope Unparked.");
        TrackState = SCOPE_TRACKING;
        SetParked(false);
        MountParkingStatusL[1].s = IPS_OK;
        MountParkingStatusL[0].s = IPS_IDLE;
        IDSetLight(&MountParkingStatusLP, nullptr);

        return true;
    } 
    else 
    {
        LOGF_ERROR("%s: Unpark failed.", getDeviceName());
        return false;
    }
}

/**
 * @brief Determine the LST with format HHMMSS
 * @return LST value for the current scope locateion
 */

bool LX200StarGo::getLST_String(char* input) 
{
    LOG_DEBUG("getLST_String");
    double siteLong;

    // step one: determine site longitude
    if (!getSiteLongitude(&siteLong))
    {
        LOGF_WARN("%s getLST Failed to get site Longitude from device.", getDeviceName());
        return false;
    }
/*
    double SD = ln_get_apparent_sidereal_time(ln_get_julian_from_sys()) - (360.0 - siteLong) / 15.0;
    double lst =  range24(SD);
    int h=0, m=0, s=0;
*/
/**/
    // determine local sidereal time
    double lst = get_local_sidereal_time(siteLong);
    int h=0, m=0, s=0;
    LOGF_DEBUG("%s Current local sidereal time = %.8f", lst, getDeviceName());
//*/
    // translate into hh:mm:ss
    getSexComponents(lst, &h, &m, &s);

    sprintf(input, "%02d%02d%02d", h, m, s);
    return true;
}


/*********************************************************************************
 * config file
 *********************************************************************************/

bool LX200StarGo::saveConfigItems(FILE *fp)
{
    LOG_DEBUG("saveConfigItems");
    LX200Telescope::saveConfigItems(fp);

    IUSaveConfigText(fp, &SiteNameTP);

    return true;
}

/*********************************************************************************
 * Queries
 *********************************************************************************/

/**
 * @brief Send a LX200 query to the communication port and read the result.
 * @param cmd LX200 query
 * @param response answer
 * @return true if the command succeeded, false otherwise
 */
bool LX200StarGo::sendQuery(const char* cmd, char* response, bool wait) 
{
    LOGF_DEBUG("%s: Sendquery %s %s", getDeviceName(), cmd, wait?"WAIT":"NOWAIT");
//    motorsState = speedState = nrTrackingSpeed = 0;
    response[0] = '\0';
    char lresponse[AVALON_RESPONSE_BUFFER_LENGTH];
    int lbytes=0;
    lresponse [0] = '\0';
    while (receive(lresponse, &lbytes, false)) 
    {
        LOGF_DEBUG("%s: Found unflushed response %s", getDeviceName(), lresponse);
        lbytes=0;
        if(ParseMotionState(lresponse))
        {
            LOGF_DEBUG("%s: Motion state response parsed", getDeviceName());            
        } 
        lresponse [0] = '\0';
    }
    flush();
    if(!transmit(cmd)) 
    {
        LOGF_ERROR("%s: query <%s> failed.", getDeviceName(), cmd);
        return false;
    }
    lresponse[0] = '\0';
    bool lwait = wait;
    while (receive(lresponse, &lbytes, lwait)) 
    {
        LOGF_DEBUG("%s: Found response %s %s", getDeviceName(), wait?"WAIT":"NOWAIT", lresponse);
        lbytes=0;
        if(ParseMotionState(lresponse))
        {
            LOGF_DEBUG("%s: Motion state response parsed", getDeviceName());            
        }
        else // Don't change wait requirement 
        {
            lwait = false;
        }
    }
    flush();
    strcpy(response, lresponse);
    return true;
}

bool LX200StarGo::ParseMotionState(char* state)
{
    LOG_DEBUG("ParseMotionState");
    int lmotor, lmode, lslew;
    if(sscanf(state, ":Z1%01d%01d%01d", &lmotor, &lmode, &lslew)==3)
    {
        LOGF_DEBUG("%s: Mount motion state %s=>%d,%d,%d", getDeviceName(),
                state, lmotor, lmode, lslew); 
    // m = 0 both motors are OFF (no power)
    // m = 1 RA motor OFF DEC motor ON
    // m = 2 RA motor ON DEC motor OFF
    // m = 3 both motors are ON
        switch(lmotor)
        {
            case 0:
                CurrentMotorsState = MOTORS_OFF;
                break;
            case 1:
                CurrentMotorsState = MOTORS_DEC_ONLY;
                break;
            case 2:
                CurrentMotorsState = MOTORS_RA_ONLY;
                break;
            case 3:
                CurrentMotorsState = MOTORS_ON;
                break;
        };
    // Tracking modes
    // t = 0 no tracking at all
    // t = 1 tracking at moon speed
    // t = 2 tracking at sun speed
    // t = 3 tracking at stars speed (sidereal speed)        
        switch(lmode)
        {
            case 0:
                CurrentTrackMode = TRACK_NONE;
                break;
            case 1:
                CurrentTrackMode = TRACK_LUNAR;
                break;
            case 2:
                CurrentTrackMode = TRACK_SOLAR;
                break;
            case 3:
                CurrentTrackMode = TRACK_SIDEREAL;
                break;
        };
    // Slew speed index
    // s = 0 GUIDE speed
    // s = 1 CENTERING speed
    // s = 2 FINDING speed
    // s = 3 MAX speed
        switch(lslew)
        {
            case 0:
                CurrentSlewRate = SLEW_GUIDE;
                break;
            case 1:
                CurrentSlewRate = SLEW_CENTERING;
                break;
            case 2:
                CurrentSlewRate = SLEW_FIND;
                break;
            case 3:
                CurrentSlewRate = SLEW_MAX;
                break;
        };
        return true;
    } 
    else
    {
        return false;
    }
}
bool LX200StarGo::querySendMountSetPark() 
{
    LOG_DEBUG("querySendMountSetPark");
    // Command  - :X352#
    // Response - 0#
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":X352#", response)) 
    {
        LOGF_ERROR("%s: Failed to send mount set park position command.", getDeviceName());
        return false;
    }
    if (response[0] != '0') 
    {
        LOGF_ERROR("%s: Invalid mount set park position response '%s'.", getDeviceName(), response);
        return false;
    }
    return true;
}

/*
 * Determine the site longitude. In contrast to a standard LX200 implementation,
 * StarGo returns the location in arc seconds precision.
 */
bool LX200StarGo::querySetSiteLongitude(double longitude)
{
    LOG_DEBUG("querySetSiteLongitude");
    int d, m, s;
    char command[32];
    if (longitude > 180) longitude = longitude - 360;
    if (longitude < -180) longitude = 360 + longitude;

    getSexComponents(longitude, &d, &m, &s);

    const char* format = ":Sg+%03d*%02d:%02d#";
    if (d < 0 || m < 0 || s < 0) format = ":Sg%04d*%02u:%02u#";

    snprintf(command, sizeof(command), format, d, m, s);

    LOGF_DEBUG("%s: Sending set site longitude request '%s'", getDeviceName(), command);

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    bool result = sendQuery(command, response);

    return (result);
}


/**
 * @brief Set the site latitude
 * @param latitude value
 * @return true iff the command succeeded
 */
bool LX200StarGo::querySetSiteLatitude(double Lat)
{
    LOG_DEBUG("querySetSiteLatitude");
    int d, m, s;
    char command[32];

    getSexComponents(Lat, &d, &m, &s);

    snprintf(command, sizeof(command), ":St%+03d*%02d:%02d#", d, m, s);

    LOGF_DEBUG("%s: Sending set site latitude request '%s'", getDeviceName(), command);

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    return (sendQuery(command, response));
}

/**
 * @brief LX200 query whether slewing is completed
 * @return true iff both motors are back to tracking or stopped
 */
bool LX200StarGo::queryIsSlewComplete()
{
    LOG_DEBUG("queryIsSlewComplete");
    // Command  - :X34#
    // the StarGo replies mxy# where x is the RA / AZ motor status and y
    // the DEC / ALT motor status meaning:
    //    x (y) = 0 motor x (y) stopped or unpowered
    //             (use :X3C# if you want  distinguish if stopped or unpowered)
    //    x (y) = 1 motor x (y) returned in tracking mode
    //    x (y) = 2 motor x (y) acelerating
    //    x (y) = 3 motor x (y) decelerating
    //    x (y) = 4 motor x (y) moving at low speed to refine
    //    x (y) = 5 motor x (y) movig at high speed to target

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":X34#", response)) 
    {
        LOGF_ERROR("%s: Failed to send query slewing state command.", getDeviceName());
        return false;
    }
    int x, y;
    int returnCode = sscanf(response, "m%01d%01d", &x, &y);
    if (returnCode <= 0) 
    {
       LOGF_ERROR("%s: Failed to parse query slewing state response '%s'.", getDeviceName(), response);
       return false;
    }

    LOGF_DEBUG("%s: slewing state motors = (%d, %d)", getDeviceName(), x, y);

    return (x < 2 && y < 2);
}

/**
 * @brief Query the motion state of the mount.
 * @param motorsState - tracking status of RA and DEC motor
 * @param speedState - 0 no tracking at all, 1 tracking at moon speed
 *        2 tracking at sun speed, 3 tracking at stars speed (sidereal speed)
 * @param nrTrackingSpeed
 * @return true if the command succeeded, false otherwise
 */
//bool LX200StarGo::queryMountMotionState(int* motorsState, int* speedState, int* nrTrackingSpeed) 
bool LX200StarGo::queryMountMotionState() 
{
    LOG_DEBUG("queryMountMotionState");
    // Command  - :X3C#
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":X3C#", response, false)) 
    {
        LOGF_ERROR("%s: Failed to send query mount motion state command.", getDeviceName());
        return false;
    }
    return true; 
}

/**
 * @brief Check whether the mount is synched or parked.
 * @param enable if true, tracking is enabled
 * @return true if the command succeeded, false otherwise
 */
bool LX200StarGo::queryParkSync (bool* isParked, bool* isSynched) 
{
    LOG_DEBUG("queryParkSync");
    // Command   - :X38#
    // Answer unparked         - p0
    // Answer at home position - p1
    // Answer parked           - p2

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":X38#", response)) 
    {
        LOGF_ERROR("%s: Failed to send get parking status request.", getDeviceName());
        return false;
    }
    int answer = 0;
    if (! sscanf(response, "p%01d", &answer)) 
    {
        LOGF_ERROR("%s: Unexpected parking status response '%s'.", getDeviceName(), response);
        return false;
    }

    switch (answer) 
    {
    case 0: (*isParked) = false; (*isSynched) = false; break;
    case 1: (*isParked) = false; (*isSynched) = true; break;
    case 2: (*isParked) = true; (*isSynched) = true; break;
    }
    return true;
}

bool LX200StarGo::querySetTracking (bool enable) 
{
    LOG_DEBUG("querySetTracking");
    // Command tracking on  - :X122#
    //         tracking off - :X120#

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (! sendQuery(enable ? ":X122#" : ":X120#", response, false)) 
    {
        LOGF_ERROR("%s: Failed to send query for %s tracking.", getDeviceName(), enable ? "enable" : "disable");
        return false;
    }
    return true;
}

/**
 * @brief Check if the ST4 port is enabled
 * @param isEnabled - true iff the ST4 port is enabled
 * @return
 */
bool LX200StarGo::queryGetST4Status (bool *isEnabled) 
{
    LOG_DEBUG("queryGetST4Status");
    // Command query ST4 status  - :TTGFh#
    //         response enabled  - vh1
    //                  disabled - vh0

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};

    if (!sendQuery(":TTGFh#", response)) 
    {
        LOGF_ERROR("%s: Failed to send query ST4 status request.", getDeviceName());
        return false;
    }
    int answer = 0;
    if (! sscanf(response, "vh%01d", &answer)) 
    {
        LOGF_ERROR("%s: Unexpected ST4 status response '%s'.", getDeviceName(), response);
        return false;
    }

    *isEnabled = (answer == 1);
    return true;
}

/**
 * @brief Determine the guiding speeds for RA and DEC axis
 * @param raSpeed percentage for RA axis
 * @param decSpeed percenage for DEC axis
 * @return
 */
bool LX200StarGo::queryGetGuidingSpeeds (int *raSpeed, int *decSpeed) 
{
    LOG_DEBUG("queryGetGuidingSpeeds");
    // Command query guiding speeds  - :X22#
    //         response              - rrbdd#
    //         rr RA speed percentage, dd DEC speed percentage

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};

    if (!sendQuery(":X22#", response)) 
    {
        LOGF_ERROR("%s: Failed to send query guiding speeds request.", getDeviceName());
        return false;
    }
    if (! sscanf(response, "%02db%2d", raSpeed, decSpeed)) 
    {
        LOGF_ERROR("%s: Unexpected guiding speed response '%s'.", getDeviceName(), response);
        return false;
    }

    return true;
}

/**
 * @brief Set the guiding speeds for RA and DEC axis
 * @param raSpeed percentage for RA axis
 * @param decSpeed percenage for DEC axis
 * @return
 */
bool LX200StarGo::setGuidingSpeeds (int raSpeed, int decSpeed) 
{
    LOG_DEBUG("setGuidingSpeeds");
    // in RA guiding speed  -  :X20rr#
    // in DEC guiding speed - :X21dd#

    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};

    sprintf(cmd, ":X20%2d#", raSpeed);
    if (sendQuery(cmd, response)) 
    {
        LOGF_INFO("Setting RA speed to %2d%%.", raSpeed);
    } 
    else 
    {
        LOGF_ERROR("%s Setting RA speed to %2d %% FAILED", getDeviceName(), raSpeed);
        return false;
    }

    sprintf(cmd, ":X21%2d#", decSpeed);
    if (sendQuery(cmd, response)) 
    {
        LOGF_INFO("Setting DEC speed to %2d%%.", decSpeed);
    } 
    else 
    {
        LOGF_ERROR("%s Setting DEC speed to %2d%% FAILED", getDeviceName(), decSpeed);
        return false;
    }
    return true;
}

/**
 * @brief Enable or disable the ST4 guiding port
 * @param enabled flag whether enable or disable
 * @return
 */

bool LX200StarGo::setST4Enabled(bool enabled) 
{
    LOG_DEBUG("setST4Enabled");

    const char *cmd = enabled ? ":TTSFh#" : ":TTRFh#";
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (sendQuery(cmd, response)) 
    {
        LOG_INFO(enabled ? "ST4 port enabled." : "ST4 port disabled.");
        return true;
    } 
    else 
    {
        LOGF_ERROR("%s Setting ST4 port FAILED", getDeviceName());
        return false;
    }
}

/**
 * @brief Determine whether the meridian flip is enabled
 * @param isEnabled - true iff flip is enabled
 * @return true iff check succeeded
 */
bool LX200StarGo::queryGetMeridianFlipEnabledStatus (bool *isEnabled) 
{
    LOG_DEBUG("queryGetMeridianFlipEnabledStatus");
    // Command query meridian flip enabled status  - :TTGFs#
    //                           response enabled  - vs0
    //                                    disabled - vs1
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":TTGFs#", response)) 
    {
        LOGF_ERROR("%s: Failed to send query meridian flip enabled status request.", getDeviceName());
        return false;
    }
    int answer = 0;
    if (! sscanf(response, "vs%01d", &answer)) 
    {
        LOGF_ERROR("%s: Unexpected meridian flip enabled status response '%s'.", getDeviceName(), response);
        return false;
    }

    *isEnabled = (answer == 0);
    return true;
}

/**
 * @brief Enabling and disabling meridian flip
 * @param enabled - setting enabled iff true
 * @return
 */
bool LX200StarGo::setMeridianFlipEnabled(bool enabled) 
{
    LOG_DEBUG("setMeridianFlipEnabled");

    const char *cmd = enabled ? ":TTRFs#" : ":TTSFs#";
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    bool success = sendQuery(cmd, response);

    if (success)
    {
        if (enabled) 
            LOG_INFO(enabled ? "Meridian flip enabled." : "Meridian flip disabled.");
        else 
            LOG_WARN("Meridian flip disabled. BE CAREFUL, THIS MAY CAUSE DAMAGE TO YOUR MOUNT!");

        MeridianFlipEnabledS[0].s = enabled ? ISS_ON  : ISS_OFF;
        MeridianFlipEnabledS[1].s = enabled ? ISS_OFF : ISS_ON;
        MeridianFlipEnabledSP.s = IPS_OK;
    } 
    else 
    {
        LOGF_ERROR("%s Setting Meridian flip FAILED", getDeviceName());
        MeridianFlipEnabledSP.s = IPS_ALERT;
    }

    IDSetSwitch(&MeridianFlipEnabledSP, nullptr);
    return success;
}

/**
 * @brief Determine whether the forcing meridian flip is enabled
 * @param isEnabled - true iff forcing flip is enabled
 * @return true iff check succeeded
 */
bool LX200StarGo::queryGetMeridianFlipForcedStatus (bool *isEnabled) 
{
    LOG_DEBUG("queryGetMeridianFlipForcedStatus");
    // Command query meridian flip enabled status  - :TTGFd#
    //                           response enabled  - vd1
    //                                    disabled - vd0

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":TTGFd#", response)) 
    {
        LOGF_ERROR("%s: Failed to send query meridian flip forced status request.", getDeviceName());
        return false;
    }
    int answer = 0;
    if (! sscanf(response, "vd%01d", &answer)) 
    {
        LOGF_ERROR("%s: Unexpected meridian flip forced status response '%s'.", getDeviceName(), response);
        return false;
    }

    *isEnabled = (answer == 1);
    return true;
}

/**
 * @brief Enabling and disabling forced meridian flip
 * @param enabled - setting enabled iff true
 * @return
 */
bool LX200StarGo::setMeridianFlipForced(bool enabled) 
{
    LOG_DEBUG("setMeridianFlipForced");
    const char *cmd = enabled ? ":TTSFd#" : ":TTRFd#";
    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    bool success = sendQuery(cmd, response);

    if (success)
    {
        if (enabled) LOG_WARN("Meridian flip forced. BE CAREFUL, THIS MAY CAUSE DAMAGE TO YOUR MOUNT!");
        else LOG_INFO("Meridian flip automatic.");

        MeridianFlipForcedS[0].s = enabled ? ISS_OFF : ISS_ON;
        MeridianFlipForcedS[1].s = enabled ? ISS_ON : ISS_OFF;
        MeridianFlipForcedSP.s = IPS_OK;
    } 
    else 
    {
        LOGF_ERROR("%s Forcing Meridian flip FAILED", getDeviceName() );
        MeridianFlipForcedSP.s = IPS_ALERT;
    }

    IDSetSwitch(&MeridianFlipForcedSP, nullptr);
    return success;
}

/**
 * @brief Retrieve pier side of the mount and sync it back to the client
 * @return true iff synching succeeds
 */
bool LX200StarGo::syncSideOfPier() 
{
    LOG_DEBUG("syncSideOfPier");
    // Command query side of pier - :X39#
    //         side unknown       - PX#
    //         east pointing west - PE#
    //         west pointing east - PW#

    char response[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":X39#", response))
    {
        LOGF_ERROR("%s: Failed to send query pier side.", getDeviceName());
        return false;
    }
    char answer;

    if (! sscanf(response, "P%c", &answer)) 
    {
        LOGF_ERROR("%s: Unexpected query pier side response '%s'.", getDeviceName(), response);
        return false;
    }

    switch (answer) 
    {
    case 'X':
        LOGF_DEBUG("%s: Detected pier side unknown.", getDeviceName());
        setPierSide(INDI::Telescope::PIER_UNKNOWN);
        break;
    case 'W':
        // seems to be vice versa
        LOGF_DEBUG("%s: Detected pier side west.", getDeviceName());
        setPierSide(INDI::Telescope::PIER_WEST);
        break;
    case 'E':
        LOGF_DEBUG("%s: Detected pier side east.", getDeviceName());
        setPierSide(INDI::Telescope::PIER_EAST);
        break;
    default:
        break;
    }

    return true;
}

/**
 * @brief Retrieve the firmware info from the mount
 * @param firmwareInfo - firmware description
 * @return
 */
bool LX200StarGo::queryFirmwareInfo (char* firmwareInfo) 
{
    LOG_DEBUG("queryFirmwareInfo");
    std::string infoStr;
    char manufacturer[AVALON_RESPONSE_BUFFER_LENGTH] = {0};

    // step 1: retrieve manufacturer
    if (!sendQuery(":GVP#", manufacturer))
    {
        LOGF_ERROR("%s: Failed to send get manufacturer request.", getDeviceName());
        return false;
    }
    infoStr.assign(manufacturer); //, bytesReceived);

    // step 2: retrieve firmware version
    char firmwareVersion[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":GVN#", firmwareVersion)) 
    {
        LOGF_ERROR("%s: Failed to send get firmware version request.", getDeviceName());
        return false;
    }
    infoStr.append(" - ").append(firmwareVersion); //, bytesReceived -1);

    // step 3: retrieve firmware date
    char firmwareDate[AVALON_RESPONSE_BUFFER_LENGTH] = {0};
    if (!sendQuery(":GVD#", firmwareDate)) 
    {
        LOGF_ERROR("%s: Failed to send get firmware date request.", getDeviceName());
        return false;
    }
    infoStr.append(" - ").append(firmwareDate); //, 1, bytesReceived);

    strcpy(firmwareInfo, infoStr.c_str());

    return true;
}

/*********************************************************************************
 * Helper functions
 *********************************************************************************/


/**
 * @brief Receive answer from the communication port.
 * @param buffer - buffer holding the answer
 * @param bytes - number of bytes contained in the answer
 * @author CanisUrsa
 * @return true if communication succeeded, false otherwise
 */
bool LX200StarGo::receive(char* buffer, int* bytes, bool wait) 
{
//    LOG_DEBUG("receive");
    int timeout = wait? AVALON_TIMEOUT: 0;
    int returnCode = tty_read_section(PortFD, buffer, '#', timeout, bytes);
    if (returnCode != TTY_OK) 
    {
        char errorString[MAXRBUF];
        tty_error_msg(returnCode, errorString, MAXRBUF);
        if(returnCode==TTY_TIME_OUT && !wait) return false;
        LOGF_WARN("%s Failed to receive full response: %s. (Return code: %d)", getDeviceName(), errorString, returnCode);
        return false;
    }
    if(buffer[*bytes-1]=='#')
        buffer[*bytes - 1] = '\0'; // remove #
    else
        buffer[*bytes] = '\0';
    LOGF_DEBUG("%s Receive response: (%d bytes) '%s'", getDeviceName(), *bytes, buffer);

    return true;
}

/**
 * @brief Flush the communication port.
 * @author CanisUrsa
 */
void LX200StarGo::flush() 
{
//    LOG_DEBUG("flush");
//    tcflush(PortFD, TCIOFLUSH);
}

bool LX200StarGo::transmit(const char* buffer) 
{
//    LOG_DEBUG("transmit");
    int bytesWritten = 0;
    flush();
    int returnCode = tty_write_string(PortFD, buffer, &bytesWritten);
    if (returnCode != TTY_OK) 
    {
        char errorString[MAXRBUF];
        tty_error_msg(returnCode, errorString, MAXRBUF);
        LOGF_WARN("%s: Failed to transmit %s. Wrote %d bytes and got error %s.", getDeviceName(), buffer, bytesWritten, errorString);
        return false;
    }
    LOGF_DEBUG("%s Sending Command '%s'", getDeviceName(), buffer);
    return true;
}

bool LX200StarGo::SetTrackMode(uint8_t mode)
{
    LOG_DEBUG("SetTrackMode");
    if (isSimulation())
        return true;
   
    LOGF_DEBUG("%s: Set Track Mode %d", getDeviceName(), mode);
    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH];
    char s_mode[10]={0};

    switch (mode)
    {
        case TRACK_SIDEREAL:
            strcpy(cmd, ":TQ#");
            strcpy(s_mode, "Sidereal");
            break;
        case TRACK_SOLAR:
            strcpy(cmd, ":TS#");
            strcpy(s_mode, "Solar");
            break;
        case TRACK_LUNAR:
            strcpy(cmd, ":TL#");
            strcpy(s_mode, "Lunar");
            break;
        case TRACK_NONE:
            strcpy(cmd, ":TM#");
            strcpy(s_mode, "None");
            break;
        default:
            return false;
            break;
    }
    if ( !sendQuery(cmd, response, false))  // Dont wait for response - there is none
        return false;
    LOGF_INFO("Tracking mode set to %s", s_mode );

// Only update tracking frequency if it is defined and not deleted by child classes
    if (genericCapability & LX200_HAS_TRACKING_FREQ)
    {
        LOGF_DEBUG("%s: Get Tracking Freq", getDeviceName());
        getTrackFrequency(&TrackFreqN[0].value);
        IDSetNumber(&TrackingFreqNP, nullptr);
    }
    return true;
}
bool LX200StarGo::checkLX200Format()
{
    char response[AVALON_RESPONSE_BUFFER_LENGTH];

    controller_format = LX200_LONG_FORMAT;
//    ::controller_format = LX200_LONG_FORMAT;
    LOGF_DEBUG("%s: checkLX200Format", getDeviceName());

    if (!sendQuery(":GR#", response)) 
    {
        LOGF_ERROR("%s: Failed to get RA for format check", getDeviceName());
        return false;
    }
    /* If it's short format, try to toggle to high precision format */
    if (strlen(response)<= 5 || response[5] == '.')
    {
        LOG_INFO("Detected low precision format, "
            "attempting to switch to high precision.");
        if (!sendQuery(":U#", response, false)) 
        {
            LOGF_ERROR("%s: Failed to switch precision", getDeviceName());
            return false;
        }
        if (!sendQuery(":GR#", response)) 
        {
            LOGF_ERROR("%s: Failed to get high precision RA", getDeviceName());
            return false;
        }
    }
    if (strlen(response)<= 5 || response[5] == '.')
    {
        controller_format = LX200_SHORT_FORMAT;
        LOG_INFO("Coordinate format is low precision.");
        return 0;
    
    }
    else if (strlen(response)> 8 && response[8] == '.')
    {
        controller_format = LX200_LONGER_FORMAT;
        LOG_INFO("Coordinate format is ultra high precision.");
        return 0;
    }
    else
    {
        controller_format = LX200_LONG_FORMAT;
        LOG_INFO("Coordinate format is high precision.");
        return 0;
    }
    return 0;
}

bool LX200StarGo::SetSlewRate(int index)
{
    LOG_DEBUG("SetSlewRate");
    // Convert index to Meade format
    index = 3 - index;

    if (!isSimulation() && !setSlewMode(index))
    {
        SlewRateSP.s = IPS_ALERT;
        IDSetSwitch(&SlewRateSP, "Error setting slew mode.");
        return false;
    }

    SlewRateSP.s = IPS_OK;
    IDSetSwitch(&SlewRateSP, nullptr);
    return true;
}
bool LX200StarGo::setSlewMode(int slewMode)
{
    LOG_DEBUG("setSlewMode");
    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH];

    switch (slewMode)
    {
        case LX200_SLEW_MAX:
            strcpy(cmd, ":RS#");
            break;
        case LX200_SLEW_FIND:
            strcpy(cmd, ":RM#");
            break;
        case LX200_SLEW_CENTER:
            strcpy(cmd, ":RC#");
            break;
        case LX200_SLEW_GUIDE:
            strcpy(cmd, ":RG#");
            break;
        default:
            return false;
            break;
    }
    if (!sendQuery(cmd, response, false)) // Don't wait for response - there isn't one
    {
        return false;
    }
    return true;
}
IPState LX200StarGo::GuideNorth(uint32_t ms)
{
    LOG_DEBUG("GuideNorth");
    if (!usePulseCommand && (MovementNSSP.s == IPS_BUSY || MovementWESP.s == IPS_BUSY))
    {
        LOGF_ERROR("%s Cannot guide while moving.", getDeviceName());
        return IPS_ALERT;
    }

    // If already moving (no pulse command), then stop movement
    if (MovementNSSP.s == IPS_BUSY)
    {
        int dir = IUFindOnSwitchIndex(&MovementNSSP);

        MoveNS(dir == 0 ? DIRECTION_NORTH : DIRECTION_SOUTH, MOTION_STOP);
    }

    if (GuideNSTID)
    {
        IERmTimer(GuideNSTID);
        GuideNSTID = 0;
    }

    if (usePulseCommand)
    {
        SendPulseCmd(LX200_NORTH, ms);
    }
    else
    {
        if (!setSlewMode(LX200_SLEW_GUIDE))
        {
            SlewRateSP.s = IPS_ALERT;
            IDSetSwitch(&SlewRateSP, "Error setting slew mode.");
            return IPS_ALERT;
        }

        MovementNSS[DIRECTION_NORTH].s = ISS_ON;
        MoveNS(DIRECTION_NORTH, MOTION_START);
    }

    // Set slew to guiding
    IUResetSwitch(&SlewRateSP);
    SlewRateS[SLEW_GUIDE].s = ISS_ON;
    IDSetSwitch(&SlewRateSP, nullptr);
    guide_direction_ns = LX200_NORTH;
    GuideNSTID      = IEAddTimer(ms, guideTimeoutHelperNS, this);
    return IPS_BUSY;
}

IPState LX200StarGo::GuideSouth(uint32_t ms)
{
    LOG_DEBUG("GuideSouth");
    if (!usePulseCommand && (MovementNSSP.s == IPS_BUSY || MovementWESP.s == IPS_BUSY))
    {
        LOGF_ERROR("%s Cannot guide while moving.", getDeviceName());
        return IPS_ALERT;
    }

    // If already moving (no pulse command), then stop movement
    if (MovementNSSP.s == IPS_BUSY)
    {
        int dir = IUFindOnSwitchIndex(&MovementNSSP);

        MoveNS(dir == 0 ? DIRECTION_NORTH : DIRECTION_SOUTH, MOTION_STOP);
    }

    if (GuideNSTID)
    {
        IERmTimer(GuideNSTID);
        GuideNSTID = 0;
    }

    if (usePulseCommand)
    {
        SendPulseCmd(LX200_SOUTH, ms);
    }
    else
    {
        if (!setSlewMode(LX200_SLEW_GUIDE))
        {
            SlewRateSP.s = IPS_ALERT;
            IDSetSwitch(&SlewRateSP, "Error setting slew mode.");
            return IPS_ALERT;
        }

        MovementNSS[DIRECTION_SOUTH].s = ISS_ON;
        MoveNS(DIRECTION_SOUTH, MOTION_START);
    }

    // Set slew to guiding
    IUResetSwitch(&SlewRateSP);
    SlewRateS[SLEW_GUIDE].s = ISS_ON;
    IDSetSwitch(&SlewRateSP, nullptr);
    guide_direction_ns = LX200_SOUTH;
    GuideNSTID      = IEAddTimer(ms, guideTimeoutHelperNS, this);
    return IPS_BUSY;
}

IPState LX200StarGo::GuideEast(uint32_t ms)
{
    LOG_DEBUG("GuideEast");
    if (!usePulseCommand && (MovementNSSP.s == IPS_BUSY || MovementWESP.s == IPS_BUSY))
    {
        LOGF_ERROR("%s Cannot guide while moving.", getDeviceName());
        return IPS_ALERT;
    }

    // If already moving (no pulse command), then stop movement
    if (MovementWESP.s == IPS_BUSY)
    {
        int dir = IUFindOnSwitchIndex(&MovementWESP);

        MoveWE(dir == 0 ? DIRECTION_WEST : DIRECTION_EAST, MOTION_STOP);
    }

    if (GuideWETID)
    {
        IERmTimer(GuideWETID);
        GuideWETID = 0;
    }

    if (usePulseCommand)
    {
        SendPulseCmd(LX200_EAST, ms);
    }
    else
    {
        if (!setSlewMode(LX200_SLEW_GUIDE))
        {
            SlewRateSP.s = IPS_ALERT;
            IDSetSwitch(&SlewRateSP, "Error setting slew mode.");
            return IPS_ALERT;
        }

        MovementWES[DIRECTION_EAST].s = ISS_ON;
        MoveWE(DIRECTION_EAST, MOTION_START);
    }

    // Set slew to guiding
    IUResetSwitch(&SlewRateSP);
    SlewRateS[SLEW_GUIDE].s = ISS_ON;
    IDSetSwitch(&SlewRateSP, nullptr);
    guide_direction_we = LX200_EAST;
    GuideWETID      = IEAddTimer(ms, guideTimeoutHelperWE, this);
    return IPS_BUSY;
}

IPState LX200StarGo::GuideWest(uint32_t ms)
{
    LOG_DEBUG("GuideWest");
    if (!usePulseCommand && (MovementNSSP.s == IPS_BUSY || MovementWESP.s == IPS_BUSY))
    {
        LOGF_ERROR("%s Cannot guide while moving.", getDeviceName());
        return IPS_ALERT;
    }

    // If already moving (no pulse command), then stop movement
    if (MovementWESP.s == IPS_BUSY)
    {
        int dir = IUFindOnSwitchIndex(&MovementWESP);

        MoveWE(dir == 0 ? DIRECTION_WEST : DIRECTION_EAST, MOTION_STOP);
    }

    if (GuideWETID)
    {
        IERmTimer(GuideWETID);
        GuideWETID = 0;
    }

    if (usePulseCommand)
    {
        SendPulseCmd(LX200_WEST, ms);
    }
    else
    {
        if (!setSlewMode(LX200_SLEW_GUIDE))
        {
            SlewRateSP.s = IPS_ALERT;
            IDSetSwitch(&SlewRateSP, "Error setting slew mode.");
            return IPS_ALERT;
        }

        MovementWES[DIRECTION_WEST].s = ISS_ON;
        MoveWE(DIRECTION_WEST, MOTION_START);
    }

    // Set slew to guiding
    IUResetSwitch(&SlewRateSP);
    SlewRateS[SLEW_GUIDE].s = ISS_ON;
    IDSetSwitch(&SlewRateSP, nullptr);
    guide_direction_we = LX200_WEST;
    GuideWETID      = IEAddTimer(ms, guideTimeoutHelperWE, this);
    return IPS_BUSY;
}

int LX200StarGo::SendPulseCmd(int8_t direction, uint32_t duration_msec)
{
    LOGF_DEBUG("%s SendPulseCmd dir=%d dur=%d ms", getDeviceName(), 
            direction, duration_msec );
    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH];
    switch (direction)
    {
        case LX200_NORTH:
            sprintf(cmd, ":Mgn%04d#", duration_msec);
            break;
        case LX200_SOUTH:
            sprintf(cmd, ":Mgs%04d#", duration_msec);
            break;
        case LX200_EAST:
            sprintf(cmd, ":Mge%04d#", duration_msec);
            break;
        case LX200_WEST:
            sprintf(cmd, ":Mgw%04d#", duration_msec);
            break;
        default:
            return 1;
    }
    if (!sendQuery(cmd, response, false)) // Don't wait for response - there isn't one
    {
        return false;
    }
    return true;
}

bool LX200StarGo::SetTrackEnabled(bool enabled)
{
    LOG_DEBUG("SetTrackEnabled");
    LOGF_INFO("Tracking being %s", enabled?"enabled":"disabled");
    return querySetTracking(enabled);
}
bool LX200StarGo::SetTrackRate(double raRate, double deRate)
{
    LOG_DEBUG("SetTrackRate");
    INDI_UNUSED(raRate);
    INDI_UNUSED(deRate);
    LOG_WARN("Custom tracking rates is not supported.");
    return false;
}
void LX200StarGo::ISGetProperties(const char *dev)
{
    if (dev != nullptr && strcmp(dev, getDeviceName()) != 0)
        return;

    LX200Telescope::ISGetProperties(dev);
    if (isConnected())
    {
        if (HasTrackMode() && TrackModeS != nullptr)
            defineSwitch(&TrackModeSP);
        if (CanControlTrack())
            defineSwitch(&TrackStateSP);
//        if (HasTrackRate())
//            defineNumber(&TrackRateNP);
    }
/*
    if (isConnected())
    {
        if (genericCapability & LX200_HAS_ALIGNMENT_TYPE)
            defineSwitch(&AlignmentSP);

        if (genericCapability & LX200_HAS_TRACKING_FREQ)
            defineNumber(&TrackingFreqNP);

        if (genericCapability & LX200_HAS_PULSE_GUIDING)
            defineSwitch(&UsePulseCmdSP);

        if (genericCapability & LX200_HAS_SITES)
        {
            defineSwitch(&SiteSP);
            defineText(&SiteNameTP);
        }

        defineNumber(&GuideNSNP);
        defineNumber(&GuideWENP);

        if (genericCapability & LX200_HAS_FOCUS)
        {
            defineSwitch(&FocusMotionSP);
            defineNumber(&FocusTimerNP);
            defineSwitch(&FocusModeSP);
        }
    }
    */
}
bool LX200StarGo::Goto(double ra, double dec)
{
    LOG_DEBUG("Goto");
    const struct timespec timeout = {0, 100000000L};

    targetRA  = ra;
    targetDEC = dec;

//    fs_sexa(RAStr, targetRA, 2, fracbase);
//    fs_sexa(DecStr, targetDEC, 2, fracbase);

    // If moving, let's stop it first.
    if (EqNP.s == IPS_BUSY)
    {
        if (!isSimulation() && !Abort() < 0)
        {
            AbortSP.s = IPS_ALERT;
            IDSetSwitch(&AbortSP, "Abort slew failed.");
            return false;
        }

        AbortSP.s = IPS_OK;
        EqNP.s    = IPS_IDLE;
        IDSetSwitch(&AbortSP, "Slew aborted.");
        IDSetNumber(&EqNP, nullptr);

        if (MovementNSSP.s == IPS_BUSY || MovementWESP.s == IPS_BUSY)
        {
            MovementNSSP.s = MovementWESP.s = IPS_IDLE;
            EqNP.s                          = IPS_IDLE;
            IUResetSwitch(&MovementNSSP);
            IUResetSwitch(&MovementWESP);
            IDSetSwitch(&MovementNSSP, nullptr);
            IDSetSwitch(&MovementWESP, nullptr);
        }

        // sleep for 100 mseconds
        nanosleep(&timeout, nullptr);
    }
    if(!isSimulation() && !setObjectCoords(ra,dec))
    {
         LOG_ERROR("Error setting coords for goto");
         return false;
    }
 //   char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH];

    if (!isSimulation())
    {
        int err = 0;
        if(!sendQuery(":MS#", response))
        /* Slew reads the '0', that is not the end of the slew */
//        if ((err = Slew(PortFD)))
        {
            LOG_ERROR("Error Slewing");
            slewError(err);
            return false;
        }
    }

    TrackState = SCOPE_SLEWING;
    EqNP.s     = IPS_BUSY;

//    LOGF_INFO("Slewing to RA: %s - DEC: %s", RAStr, DecStr);

    return true;
}

bool LX200StarGo::MoveNS(INDI_DIR_NS dir, TelescopeMotionCommand command)
{
    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH];
    LOGF_DEBUG("%s", __FUNCTION__);
    
    sprintf(cmd, ":%s%s#", command==MOTION_START?"M":"Q", dir == DIRECTION_NORTH?"n":"s");
    if (!isSimulation() && !sendQuery(cmd, response, false))
    {
        LOG_ERROR("Error N/S motion direction.");
        return false;
    }

    return true;
}

bool LX200StarGo::MoveWE(INDI_DIR_WE dir, TelescopeMotionCommand command)
{
    char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH];
    LOGF_DEBUG("%s", __FUNCTION__);

    sprintf(cmd, ":%s%s#", command==MOTION_START?"M":"Q", dir == DIRECTION_WEST?"n":"s");

    if (!isSimulation() && !sendQuery(cmd, response, false))
    {
        LOG_ERROR("Error W/E motion direction.");
        return false;
    }

    return true;
}

bool LX200StarGo::Abort()
{
//   char cmd[AVALON_COMMAND_BUFFER_LENGTH];
    char response[AVALON_RESPONSE_BUFFER_LENGTH];
    if (!isSimulation() && !sendQuery(":Q#", response, false))
    {
        LOG_ERROR("Failed to abort slew.");
        return false;
    }

    if (GuideNSNP.s == IPS_BUSY || GuideWENP.s == IPS_BUSY)
    {
        GuideNSNP.s = GuideWENP.s = IPS_IDLE;
        GuideNSN[0].value = GuideNSN[1].value = 0.0;
        GuideWEN[0].value = GuideWEN[1].value = 0.0;

        if (GuideNSTID)
        {
            IERmTimer(GuideNSTID);
            GuideNSTID = 0;
        }

        if (GuideWETID)
        {
            IERmTimer(GuideWETID);
            GuideNSTID = 0;
        }

        LOG_INFO("Guide aborted.");
        IDSetNumber(&GuideNSNP, nullptr);
        IDSetNumber(&GuideWENP, nullptr);

        return true;
    }

    return true;
}
bool LX200StarGo::Sync(double ra, double dec)
{
 //   char syncString[256]={0};
    char response[AVALON_RESPONSE_BUFFER_LENGTH];

    if(!isSimulation() && !setObjectCoords(ra,dec))
    {
         LOG_ERROR("Error setting coords for sync");
         return false;
    }

    if (!isSimulation() && !sendQuery(":CM#", response))
    {
        EqNP.s = IPS_ALERT;
        IDSetNumber(&EqNP, "Synchronization failed.");
        return false;
    }

    currentRA  = ra;
    currentDEC = dec;

    LOG_INFO("Synchronization successful.");

    EqNP.s     = IPS_OK;

    NewRaDec(currentRA, currentDEC);

    return true;
}

bool LX200StarGo::setObjectCoords(double ra, double dec)
{
    LOG_DEBUG("setObjectCoords");

    char RAStr[64]={0}, DecStr[64]={0};
    int h, m, s, d;
//    double d_s;

/*    switch (controller_format) //getLX200Format())
    {
    case LX200_LONGER_FORMAT:
        getSexComponentsIID(ra, &h, &m, &d_s);
        snprintf(RAStr, sizeof(RAStr), ":Sr%02d:%02d:%05.02f#", h, m, d_s);
        getSexComponentsIID(dec, &d, &m, &d_s);
        // case with negative zero
        if (!d && dec < 0)
            snprintf(DecStr, sizeof(DecStr), ":Sd-%02d*%02d:%05.02f#", d, m, d_s);
        else
            snprintf(DecStr, sizeof(DecStr), ":Sd%+03d*%02d:%05.02f#", d, m, d_s);
        break;
    case LX200_LONG_FORMAT:
*/
        getSexComponents(ra, &h, &m, &s);
        snprintf(RAStr, sizeof(RAStr), ":Sr %02d:%02d:%02d#", h, m, s);
        getSexComponents(dec, &d, &m, &s);
        /* case with negative zero */
        if (!d && dec < 0)
            snprintf(DecStr, sizeof(DecStr), ":Sd -%02d*%02d:%02d #", d, m, s);
        else
            snprintf(DecStr, sizeof(DecStr), ":Sd %+03d*%02d:%02d #", d, m, s);
/*
        break;
    case LX200_SHORT_FORMAT:
        double frac_m;
        getSexComponents(ra, &h, &m, &s);
        frac_m = m + (s / 60.0);
        snprintf(RAStr, sizeof(RAStr), ":Sr %02d:%02.1f#", h, frac_m);
        getSexComponents(dec, &d, &m, &s);
        // case with negative zero
        if (!d && dec < 0)
            snprintf(DecStr, sizeof(DecStr), ":Sd-%02d*%02d#", d, m);
        else
            snprintf(DecStr, sizeof(DecStr), ":Sd%+03d*%02d#", d, m);
        break;
    default:
        return false;
        break;
    }
*/
    char response[AVALON_RESPONSE_BUFFER_LENGTH];
    if (!isSimulation())
    {
        if(!sendQuery(RAStr, response, false)  || !sendQuery(DecStr, response, false) )
        {
            EqNP.s = IPS_ALERT;
            IDSetNumber(&EqNP, "Error setting RA/DEC.");
            return false;            
        }
    }
    return true;
}
bool LX200StarGo::setLocalDate(uint8_t days, uint8_t months, uint16_t years)
{
    char cmd[RB_MAX_LEN];
    char response[RB_MAX_LEN];

    int yy = years % 100;

    snprintf(cmd, sizeof(cmd), ":SC %02d/%02d/%02d#", months, days, yy);
    if (!sendQuery(cmd, response))
        return false;

    if (response[0] == '0')
        return false;

    return true;
}

bool LX200StarGo::setLocalTime24(uint8_t hour, uint8_t minute, uint8_t second)
{
    char cmd[RB_MAX_LEN]={0};
    char response[RB_MAX_LEN];

    snprintf(cmd, sizeof(cmd), ":SL %02d:%02d:%02d#", hour, minute, second);

    return (sendQuery(cmd, response, false));
}

bool LX200StarGo::setUTCOffset(double offset)
{
    char cmd[RB_MAX_LEN]={0};
    char response[RB_MAX_LEN];
    int hours = offset * -1.0;

    snprintf(cmd, sizeof(cmd), ":SG %+03d#", hours);

    return (sendQuery(cmd, response, false));
}

bool LX200StarGo::getLocalTime(char *timeString)
{
    if (isSimulation())
    {
        time_t now = time (nullptr);
        strftime(timeString, 32, "%T", localtime(&now));
    }
    else
    {
        double ctime=0;
        int h, m, s;
        char response[RB_MAX_LEN]={0};

        if (!sendQuery(":GL#", response))
            return false;

        if (f_scansexa(response, &ctime))
        {
            LOG_DEBUG("Unable to parse local time response");
            return false;
        }

        getSexComponents(ctime, &h, &m, &s);
        snprintf(timeString, 32, "%02d:%02d:%02d", h, m, s);
    }

    return true;
}

bool LX200StarGo::getLocalDate(char *dateString)
{
    if (isSimulation())
    {
        time_t now = time (nullptr);
        strftime(dateString, 32, "%F", localtime(&now));
    }
    else
    {
        char response[RB_MAX_LEN]={0};
        int dd, mm, yy;
        char mell_prefix[3]={0};
        int vars_read=0;

        if (!sendQuery(":GC#", response))
            return false;
        /* StarGo format is MM/DD/YY */
        vars_read = sscanf(response, "%d%*c%d%*c%d", &mm, &dd, &yy);
        if (vars_read < 3)
            return false;
        /* We consider years 50 or more to be in the last century, anything less in the 21st century.*/
        if (yy > 50)
            strncpy(mell_prefix, "19", 3);
        else
            strncpy(mell_prefix, "20", 3);
        /* We need to have it in YYYY-MM-DD ISO format */
        snprintf(dateString, 32, "%s%02d-%02d-%02d", mell_prefix, yy, mm, dd);
    }
    return true;
}

bool LX200StarGo::getUTFOffset(double *offset)
{
    if (isSimulation())
    {
        *offset = 3;
        return true;
    }

    int lx200_utc_offset = 0;
    char response[RB_MAX_LEN]={0};
    float temp_number;

    if (!sendQuery(":GG#", response))
        return false;

    /* Float */
    if (strchr(response, '.'))
    {
        if (sscanf(response, "%f", &temp_number) != 1)
            return false;
        lx200_utc_offset = static_cast<int>(temp_number);
    }
    /* Int */
    else if (sscanf(response, "%d", &lx200_utc_offset) != 1)
        return false;

    // LX200 TimeT Offset is defined at the number of hours added to LOCAL TIME to get TimeT. This is contrary to the normal definition.
    *offset = lx200_utc_offset * -1;
    return true;
}

bool LX200StarGo::getTrackFrequency(double *value)
{
    float Freq;
    char response[RB_MAX_LEN]={0};

    if (!sendQuery(":GT#", response))
        return false;

    if (sscanf(response, "%f#", &Freq) < 1)
    {
        LOG_ERROR("Unable to parse response");
        return false;
    }

    *value = static_cast<double>(Freq);
    return true;
}
