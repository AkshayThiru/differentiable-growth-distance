#ifndef DGD_BENCHMARKS_INTERNAL_HELPERS_TIMER_H_
#define DGD_BENCHMARKS_INTERNAL_HELPERS_TIMER_H_

#include <chrono>

namespace dgd {

namespace internal {

// Timer for benchmarking.
// Adapted from:
// https://github.com/coal-library/coal/blob/devel/include/coal/timings.h
class Timer {
 public:
  explicit Timer(bool start_on_construction = true) : running_(false) {
    if (start_on_construction) Start();
  }

  void Start() {
    if (!running_) {
      running_ = true;
      elapsed_ = 0.0;
      start_ = std::chrono::steady_clock::now();
    }
  }

  void Stop() {
    if (running_) {
      end_ = std::chrono::steady_clock::now();
      running_ = false;
      elapsed_ += static_cast<double>(
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          end_ - start_)
                          .count()) *
                  1e-3;
    }
  }

  void Resume() {
    if (!running_) {
      running_ = true;
      start_ = std::chrono::steady_clock::now();
    }
  }

  double Elapsed() const {
    const auto current = std::chrono::steady_clock::now();

    if (running_) {
      return elapsed_ +
             static_cast<double>(
                 std::chrono::duration_cast<std::chrono::nanoseconds>(current -
                                                                      start_)
                     .count()) *
                 1e-3;
    } else {
      return elapsed_;
    }
  }

 private:
  std::chrono::time_point<std::chrono::steady_clock> start_, end_;
  double elapsed_;
  bool running_;
};

}  // namespace internal

}  // namespace dgd

#endif  // DGD_BENCHMARKS_INTERNAL_HELPERS_TIMER_H_
