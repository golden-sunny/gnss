/**
 * Shared BDS WeekSecond helper.
 */

#ifndef BDS_WEEK_SECOND_HPP
#define BDS_WEEK_SECOND_HPP

#include "TimeStruct.h"
#include "Const.h"

class BDSWeekSecond : public WeekSecond {
public:
    BDSWeekSecond(unsigned int w = 0, double s = 0.0, TimeSystem ts = TimeSystem::BDT)
            : WeekSecond(w, s, ts) { timeSystem = ts; }

    int Nbits(void) const override {
        static const int n = 13;
        return n;
    }

    int bitmask(void) const override {
        static const int bm = 0x1FFF;
        return bm;
    }

    long MJDEpoch(void) const override {
        static const long e = BDS_EPOCH_MJD;
        return e;
    }
};

#endif // BDS_WEEK_SECOND_HPP
