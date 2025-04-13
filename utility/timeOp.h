#include <cstdint>
#include <chrono>
#ifndef TIMEOP_H
#define TIMEOP_H

using namespace std;

namespace TimeOp {
  // return nansecond from epoch time
  static int64_t getClockNan() {
    return chrono::high_resolution_clock::now().time_since_epoch().count();
  }

  // return time difference in nansecond
  static int64_t diffNan(chrono::_V2::system_clock::time_point start,
                             chrono::_V2::system_clock::time_point end) {
    return chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  }
}  // namespace TimeOp

#endif  // TIMEOP_H