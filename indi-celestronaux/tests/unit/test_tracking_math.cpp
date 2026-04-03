#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

// Mock for coordinate transforms
struct IHorizontalCoordinates {
    double azimuth;
    double altitude;
};

double range180(double angle) {
    while (angle <= -180) angle += 360;
    while (angle > 180) angle -= 360;
    return angle;
}

// 5-degree radius circumpolar circle
// omega = 15 deg/hour = 15/3600 deg/sec
const double R = 5.0;
const double omega = (15.0 / 3600.0) * (M_PI / 180.0); // rad/sec
const double h0 = 45.0;

IHorizontalCoordinates getTrueCoords(double t) {
    return { R * cos(omega * t), h0 + R * sin(omega * t) };
}

class TrackingModel {
public:
    virtual std::pair<double, double> getRates(double t, double dt, double curAz, double curAlt) = 0;
};

class LinearModel : public TrackingModel {
public:
    std::pair<double, double> getRates(double t, double dt, double curAz, double curAlt) override {
        // Real linear tracking uses the previous rate (lagging)
        auto p_1 = getTrueCoords(t - dt);
        auto p0 = getTrueCoords(t);
        return { range180(p0.azimuth - p_1.azimuth) / dt, (p0.altitude - p_1.altitude) / dt };
    }
};

class ParabolicModel : public TrackingModel {
    IHorizontalCoordinates window[3];
    bool primed = false;
public:
    std::pair<double, double> getRates(double t, double dt, double curAz, double curAlt) override {
        if (!primed) {
            window[0] = getTrueCoords(t - dt);
            window[1] = getTrueCoords(t);
            window[2] = getTrueCoords(t + dt);
            primed = true;
        } else {
            window[0] = window[1];
            window[1] = window[2];
            window[2] = getTrueCoords(t + dt);
        }

        double pAz[3] = { window[0].azimuth, window[1].azimuth, window[2].azimuth };
        double pAlt[3] = { window[0].altitude, window[1].altitude, window[2].altitude };

        pAz[0] = pAz[1] + range180(pAz[0] - pAz[1]);
        pAz[2] = pAz[1] + range180(pAz[2] - pAz[1]);

        double az_a = (pAz[2] + pAz[0] - 2 * pAz[1]) / 2.0;
        double az_b = (pAz[2] - pAz[0]) / 2.0;
        double alt_a = (pAlt[2] + pAlt[0] - 2 * pAlt[1]) / 2.0;
        double alt_b = (pAlt[2] - pAlt[0]) / 2.0;

        double targetAzNext = az_a * dt * dt + az_b * dt + pAz[1];
        double targetAltNext = alt_a * dt * dt + alt_b * dt + pAlt[1];

        return { range180(targetAzNext - curAz) / dt, (targetAltNext - curAlt) / dt };
    }
};

void runSimulation(TrackingModel& model, const std::string& name) {
    double dt = 1.0;
    auto start = getTrueCoords(0);
    double curAz = start.azimuth;
    double curAlt = start.altitude;
    double maxErrorArcsec = 0;

    for (int t = 0; t < 3600; ++t) {
        auto rates = model.getRates(t, dt, curAz, curAlt);
        
        // Update position (simulate mount movement)
        curAz += rates.first * dt;
        curAlt += rates.second * dt;

        // Compare with ground truth at T+1
        auto truth = getTrueCoords(t + 1);
        double errAz = range180(curAz - truth.azimuth) * cos(truth.altitude * M_PI / 180.0);
        double errAlt = curAlt - truth.altitude;
        double errTotal = sqrt(errAz * errAz + errAlt * errAlt) * 3600.0; // in arcsec
        
        maxErrorArcsec = std::max(maxErrorArcsec, errTotal);
    }

    std::cout << name << " - Max Path Deviation: " << std::fixed << std::setprecision(4) << maxErrorArcsec << " arcsec" << std::endl;
}

int main() {
    std::cout << "Starting 1-Hour Circumpolar Tracking Simulation..." << std::endl;
    std::cout << "Radius: 5 degrees, Rate: Sidereal (15 deg/hr)" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    LinearModel linear;
    runSimulation(linear, "LINEAR (1st Order) ");

    ParabolicModel parabolic;
    runSimulation(parabolic, "PARABOLIC (2nd Order)");

    std::cout << "--------------------------------------------------------" << std::endl;
    return 0;
}
