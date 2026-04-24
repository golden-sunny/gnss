/**
 * Exercise 3/4: Create BDS WeekSecond and convert with CommonTime
 */

#include <iostream>
#include <iomanip>
#include "TimeStruct.h"
#include "TimeConvert.h"
#include "BDSWeekSecond.hpp"

int main() {
    // Example: BDT civil time
    CivilTime bdtCivil(2020, 1, 1, 0, 0, 0.0, TimeSystem::BDT);
    CommonTime bdtCT = CivilTime2CommonTime(bdtCivil);

    BDSWeekSecond bdsWS;
    CommonTime2WeekSecond(bdtCT, bdsWS);

    CommonTime backCT;
    WeekSecond2CommonTime(bdsWS, backCT);
    CivilTime backCivil = CommonTime2CivilTime(backCT);

    std::cout << "BDT CivilTime      : " << bdtCivil << std::endl;
    std::cout << "BDS WeekSecond     : week=" << bdsWS.week << " sow=" << std::fixed << std::setprecision(3) << bdsWS.sow << std::endl;
    std::cout << "Back to CommonTime : " << backCT << std::endl;
    std::cout << "Back to CivilTime  : " << backCivil << std::endl;

    return 0;
}
