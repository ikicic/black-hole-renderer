#include <chrono>
#include <cstdio>
#include <thread>
#include <vector>

void _generate_image_status_thread(
    std::vector<int> &finished_pixels_vec,
    int width,
    int height) {

  auto start_time = std::chrono::system_clock::now();
  int done;
  do {
    done = 0;
    for (int x : finished_pixels_vec)
      done += x;
    std::chrono::duration<double> time_delta =
        std::chrono::system_clock::now() - start_time;
    fprintf(stderr, "%.2lf%% (%.0lfRPS, ETA: %ds)      \r",
        100. * done / (width * height),
        done / time_delta.count(),
        (int)((width * height - done) * time_delta.count() / done)
      );
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  } while (done < width * height);
  fprintf(stderr, "\n");
}


