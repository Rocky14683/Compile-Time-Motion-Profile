#include <iostream>
#include <cmath>
#include <map>
#include <fstream>

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Constraints {
    double maxJerk;
    double maxAcceleration;
    double maxVelocity;
};

struct TimeInterval {
    double tj1 = 0;
    double tj2 = 0;
    double ta = 0;//acceleration period
    double tv = 0;//const v stage
    double td = 0;//deceleration period

    double getTotalDuration() const {
        return this->ta + this->tv + this->td;
    }

    bool maxAccelReachable() const {
        return this->ta >= 2. * this->tj1 && this->td >= 2. * this->tj2;
    }
};


struct StartAndEndCondition {
    double p0;
    double pf;
    double v0;
    double vf;

    int dir() const {
        return sgn(this->pf - this->p0);
    }

    double h() const {
        return abs(this->pf - p0);
    }
};

struct SCurveParameters;

struct SCurveInput {
    Constraints constraint;
    StartAndEndCondition condition;

    TimeInterval calcInterval() {
        return this->calcTimeCase1();
    }

    bool isAMaxReached() const {//return true if reached, false if not reached
        return (this->constraint.maxVelocity - this->condition.v0) * this->constraint.maxJerk >=
               this->constraint.maxAcceleration * this->constraint.maxAcceleration;
    }

    bool isDMaxReached() const {
        return (this->constraint.maxVelocity - this->condition.vf) * this->constraint.maxJerk >=
               this->constraint.maxAcceleration * this->constraint.maxAcceleration;
    }


    TimeInterval calcTimeCase1() {
        TimeInterval time;
        auto *newIn = new SCurveInput(*this);
        int dir = this->condition.dir();
        if (!isAMaxReached()) {
            time.tj1 = sqrt(
                    (newIn->constraint.maxVelocity - this->condition.v0) /
                    newIn->constraint.maxJerk);
            time.ta = 2. * time.tj1;
        } else {
            time.tj1 = newIn->constraint.maxAcceleration / newIn->constraint.maxJerk;
            time.ta = time.tj1 +
                      (newIn->constraint.maxVelocity - dir * this->condition.v0) / newIn->constraint.maxAcceleration;
        }

        if (!isDMaxReached()) {
            time.tj2 = sqrt(
                    (newIn->constraint.maxVelocity - this->condition.vf) /
                    newIn->constraint.maxJerk);
            time.td = 2. * time.tj2;
        } else {
            time.tj2 = newIn->constraint.maxAcceleration / newIn->constraint.maxJerk;
            time.td = time.tj2 +
                      (newIn->constraint.maxVelocity - dir * this->condition.vf) / newIn->constraint.maxAcceleration;
        }

        time.tv = this->condition.h() / newIn->constraint.maxVelocity
                  - time.ta / 2.0
                    * (1.0 + dir * this->condition.v0 / newIn->constraint.maxVelocity)
                  - time.td / 2.0
                    * (1.0 + dir * this->condition.vf / newIn->constraint.maxVelocity);

        if (time.tv <= 0) {
            return this->calcTimeCase2(0);
        }

        if (!time.maxAccelReachable()) {
            newIn->constraint.maxAcceleration *= 0.5;
            if (newIn->constraint.maxAcceleration > 0.01) {
                return newIn->calcTimeCase2(0);
            }
            newIn->constraint.maxAcceleration = 0.;
        }
        this->negAccelTimeCase(&time, newIn);
        delete newIn;
        return time;
    }

    TimeInterval calcTimeCase2(int recursionDepth) {
        recursionDepth++;
        TimeInterval time = getTimeCase2();
        auto newIn = new SCurveInput(*this);
        if (!time.maxAccelReachable()) {
            newIn->constraint.maxAcceleration *= 0.5;
            if (newIn->constraint.maxAcceleration > 0.01) {
                return newIn->calcTimeCase2(recursionDepth);
            }
            newIn->constraint.maxAcceleration = 0.;
        }
        this->negAccelTimeCase(&time, newIn);
        if (recursionDepth != 1) {
            newIn->constraint.maxAcceleration *= 2.0;
        }
        return newIn->calcTimeCase2_precise(recursionDepth);
    }

    TimeInterval getTimeCase2() {
        TimeInterval time;
        time.tj1 = this->constraint.maxAcceleration / this->constraint.maxJerk;
        time.tj2 = this->constraint.maxAcceleration / this->constraint.maxJerk;
        time.tv = .0;
        double delta = pow(this->constraint.maxAcceleration, 4) / pow(this->constraint.maxJerk, 2)
                       + 2.0 * (this->condition.v0 * this->condition.v0 + this->condition.vf * this->condition.vf)
                       + this->constraint.maxAcceleration *
                         (4.0 * this->condition.h() - 2.0 * this->constraint.maxAcceleration / this->constraint.maxJerk
                                                      * (this->condition.v0 + this->condition.vf));
        time.ta = (this->constraint.maxAcceleration * this->constraint.maxAcceleration / this->constraint.maxJerk -
                   2.0 * this->condition.v0
                   + sqrt(delta)) / (2.0 * this->constraint.maxAcceleration);
        time.td = (this->constraint.maxAcceleration * this->constraint.maxAcceleration / this->constraint.maxJerk -
                   2.0 * this->condition.vf
                   + sqrt(delta)) / (2.0 * this->constraint.maxAcceleration);

        return time;
    }

    TimeInterval calcTimeCase2_precise(int recursionDepth) {
        recursionDepth++;
        TimeInterval time = this->getTimeCase2();
        auto newIn = new SCurveInput(*this);
        if (!time.maxAccelReachable()) {
            newIn->constraint.maxAcceleration *= 0.99;
            if (newIn->constraint.maxAcceleration > 0.01) {
                return newIn->calcTimeCase2_precise(recursionDepth);
            }
            newIn->constraint.maxAcceleration = 0.;
        }
        this->negAccelTimeCase(&time, newIn);
        return time;
    }

    void negAccelTimeCase(TimeInterval *time, SCurveInput *newIn) {
        if (time->ta < .0) {
            time->tj1 = 0;
            time->ta = 0;
            time->td = this->condition.h() / (this->condition.v0 + this->condition.vf);
            time->tj2 = (newIn->constraint.maxJerk * this->condition.h()
                         - sqrt(newIn->constraint.maxJerk
                                * (newIn->constraint.maxJerk * pow(this->condition.h(), 2)
                                   + pow((this->condition.v0 + this->condition.vf), 2)
                                     * (this->condition.vf - this->condition.v0))
            )) / (newIn->constraint.maxJerk
                  * (this->condition.vf + this->condition.v0));
        }

        if (time->td < .0) {
            time->tj2 = 0;
            time->td = 0;
            time->ta = 2.0 * this->condition.h() / (this->condition.v0 + this->condition.vf);
            time->tj2 = (newIn->constraint.maxJerk * this->condition.h()
                         - sqrt(newIn->constraint.maxJerk
                                * (newIn->constraint.maxJerk * pow(this->condition.h(), 2)
                                   - pow((this->condition.v0 + this->condition.vf), 2)
                                     * (this->condition.vf - this->condition.v0))
            )) / (newIn->constraint.maxJerk
                  * (this->condition.vf + this->condition.v0));
        }
    }

    static double getPosition(const SCurveParameters &p, double t);

    static double getVelocity(const SCurveParameters &p, double t);

    static double getAcceleration(const SCurveParameters &p, double t);

    auto generatePos(SCurveInput &p);
    auto generateVel(SCurveInput &p);
    auto generateAccel(SCurveInput &p);
};


struct SCurveParameters {
    TimeInterval timeIntervals;
    StartAndEndCondition conditions;
    double jMax;
    double jMin;
    double reachedA;//maximum achieved acceleration during the acceleration phase
    double reachedD;//maximum achieved acceleration during the deceleration phase
    double reachedV;

    SCurveParameters(TimeInterval timeInterval, SCurveInput &p) {
        this->timeIntervals = timeInterval;
        this->jMax = p.constraint.maxJerk;
        this->jMin = -p.constraint.maxJerk;
        this->reachedA = p.constraint.maxJerk * timeInterval.tj1;
        this->reachedD = -p.constraint.maxJerk * timeInterval.tj2;
        this->reachedV = p.condition.dir() * p.condition.v0 + (timeInterval.ta - timeInterval.tj1) * reachedA;
        this->conditions = p.condition;
    };
};

double SCurveInput::getPosition(const SCurveParameters &p, double t) {
//    assert(t < 0);
    auto time = p.timeIntervals;
    int dir = p.conditions.dir();
    if (t <= time.tj1) {
        return p.conditions.p0 + p.conditions.v0 * t + dir * p.jMax * pow(t, 3) / 6.;
    } else if (t <= time.ta - time.tj1) {
        return p.conditions.p0
               + p.conditions.v0 * t
               + dir * p.reachedA / 6. * (3. * pow(t, 2) - 3. * time.tj1 * t + pow(time.tj1, 2));
    } else if (t <= time.ta) {
        return p.conditions.p0 + dir * (p.reachedV + dir * p.conditions.v0) * time.ta / 2.
               - dir * p.reachedV * (time.ta - t)
               - dir * p.jMin * pow(time.ta - t, 3) / 6.;
    } else if (t <= time.ta + time.tv) {
        return p.conditions.p0
               + dir * (p.reachedV + dir * p.conditions.v0) * time.ta / 2.
               + dir * p.reachedV * (t - time.ta);
    } else if (t <= time.getTotalDuration() - time.td + time.tj2) {
        return p.conditions.pf - dir * (p.reachedV + dir * p.conditions.vf) * time.td / 2.
               + dir * p.reachedV * (t - time.getTotalDuration() + time.td)
               - dir * p.jMax * pow(t - time.getTotalDuration() + time.td, 3) / 6.;
    } else if (t <= time.getTotalDuration() - time.tj2) {
        return p.conditions.pf - dir * (p.reachedV + dir * p.conditions.vf) * time.td / 2.
               + dir * p.reachedV * (t - time.getTotalDuration() + time.td)
               + dir * p.reachedD / 6.
                 * (3. * pow(t - time.getTotalDuration() + time.td, 2)
                    - 3. * time.tj2 * (t - time.getTotalDuration() + time.td)
                    + pow(time.tj2, 2));
    } else if (t <= time.getTotalDuration()) {
        return p.conditions.pf
               - p.conditions.vf * (time.getTotalDuration() - t)
               - dir * p.jMax * pow(time.getTotalDuration() - t, 3) / 6.;
    } else {
        return p.conditions.pf;
    }

}

double SCurveInput::getVelocity(const SCurveParameters &p, double t) {
//    assert(t > 0);
    auto time = p.timeIntervals;
    int dir = p.conditions.dir();
    if (t <= time.tj1) {
        return p.conditions.v0 + dir * p.jMax * pow(t, 2) / 2.;
    } else if (t <= time.ta - time.tj1) {
        return p.conditions.v0 + dir * p.reachedA * (t - time.tj1 / 2.);
    } else if (t <= time.ta) {
        return dir * p.reachedV + dir * p.jMin * pow(time.ta - t, 2) / 2.;
    } else if (t <= time.ta + time.tv) {
        return dir * p.reachedV;
    } else if (t <= time.getTotalDuration() - time.td + time.tj2) {
        return dir * p.reachedV - dir * p.jMax * pow(t - time.getTotalDuration() + time.td, 2) / 2.;
    } else if (t <= time.getTotalDuration() - time.tj2) {
        return dir * p.reachedV + dir * p.reachedD * (t - time.getTotalDuration() + time.td - time.tj2 / 2.);
    } else if (t <= time.getTotalDuration()) {
        return p.conditions.vf + dir * p.jMax * pow(time.getTotalDuration() - t, 2) / 2.;
    } else {
        return p.conditions.vf;
    }
}

double SCurveInput::getAcceleration(const SCurveParameters &p, double t) {
//    assert(t > 0);
    auto time = p.timeIntervals;
    int dir = p.conditions.dir();
    if (t <= time.tj1) {
        return dir * p.jMax * t;
    } else if (t <= time.ta - time.tj1) {
        return dir * p.reachedA;
    } else if (t <= time.ta) {
        return dir * (-p.jMin) * (time.ta - t);
    } else if (t <= time.ta + time.tv) {
        return .0;
    } else if (t <= time.getTotalDuration() - time.td + time.tj2) {
        return dir * (-p.jMax) * (t - time.getTotalDuration() + time.td);
    } else if (t <= time.getTotalDuration() - time.tj2) {
        return dir * p.reachedD;
    } else if (t <= time.getTotalDuration()) {
        return dir * (-p.jMax) * (time.getTotalDuration() - t);
    } else {
        return .0;
    }
}

auto SCurveInput::generatePos(SCurveInput &p) {
    struct ret {
        SCurveParameters *param;
        std::function<double(const SCurveParameters &, double)> func;//t:v
    };
    auto time = p.calcInterval();
    auto params = SCurveParameters(time, p);
    auto paramsCopy = new SCurveParameters(params);
    std::function<double(const SCurveParameters &, double)> func = this->getPosition;
    return ret{paramsCopy, func};
}

auto SCurveInput::generateVel(SCurveInput &p) {
    struct ret {
        SCurveParameters *param;
        std::function<double(const SCurveParameters &, double)> func;//t:v
    };
    auto time = p.calcInterval();
    auto params = SCurveParameters(time, p);
    auto paramsCopy = new SCurveParameters(params);
    std::function<double(const SCurveParameters &, double)> func = this->getVelocity;
    return ret{paramsCopy, func};
}

auto SCurveInput::generateAccel(SCurveInput &p) {
    struct ret {
        SCurveParameters *param;
        std::function<double(const SCurveParameters &, double)> func;//t:v
    };
    auto time = p.calcInterval();
    auto params = SCurveParameters(time, p);
    auto paramsCopy = new SCurveParameters(params);
    std::function<double(const SCurveParameters &, double)> func = this->getAcceleration;
    return ret{paramsCopy, func};
}

Constraints constraint{
        .maxJerk = 50.0, .maxAcceleration = 30., .maxVelocity = 40.
};

StartAndEndCondition cond{
        .p0 = .0, .pf = 20., .v0 = 5, .vf = 2
};

SCurveInput Scurve{.constraint = constraint, .condition = cond};

int main() {
    auto [paramsP, getP] = Scurve.generatePos(Scurve);
    auto [paramsV, getV] = Scurve.generateVel(Scurve);
    auto [paramsA, getA] = Scurve.generateAccel(Scurve);
    std::ofstream Spos;
    Spos.open("/Users/rockychen/Desktop/ArbitarySCurveProfile/Position.txt");
    std::ofstream Svel;
    Svel.open("/Users/rockychen/Desktop/ArbitarySCurveProfile/velocity.txt");
    std::ofstream Sacc;
    Sacc.open("/Users/rockychen/Desktop/ArbitarySCurveProfile/acceleration.txt");


//    std::cout << paramsV->timeIntervals.getTotalDuration() << std::endl;
//    std::cout << "tj1: " << paramsV->timeIntervals.tj1 << ", tj2: " << paramsV->timeIntervals.tj2 << ", tv: "
//              << paramsV->timeIntervals.tv << ", ta: " << paramsV->timeIntervals.ta << ", td: "
//              << paramsV->timeIntervals.td << std::endl;
    for (double t = 0; t < paramsV->timeIntervals.getTotalDuration(); t += 0.01) {
        Spos << getP(*paramsP, t) << '\n';
        Svel << getV(*paramsV, t) << '\n';
        Sacc << getA(*paramsA, t) << '\n';
    }
//    Spos.close();
    Svel.close();
//    Sacc.close();
    return 0;
}
