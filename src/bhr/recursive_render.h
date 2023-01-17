#ifndef RECURSIVE_RENDER_H
#define RECURSIVE_RENDER_H

#include <map>

#include <bhr/geodesic.h>
#include <bhr/utility.h>
#include <bhr/rectangle.h>

#define LOCATION_FACTOR   (1 << 16)


template <typename _FullGeodesicData>
struct ScreenGeodesic {
  _FullGeodesicData full;
  // FullGeodesicData<BasicGeodesicState<_Coord>, GeodesicExtra<_Coord>> full;
  int x, y;  // Screen location multiplied with LOCATION_FACTOR.
};

struct QuadtreeNode {
  int geodesic_ids[4];
  int center_x, center_y;  // Location multiplied with LOCATION_FACTOR.
  int size;  // start size / 2^depth
  int parent;  // ID of the parent node.
  bool split;
};

struct QuadtreeLeaf {
  int geodesic_ids[4];

  QuadtreeLeaf() {}
  QuadtreeLeaf(const QuadtreeNode &node) {
    for (int i = 0; i < 4; ++i)
      geodesic_ids[i] = node.geodesic_ids[i];
  }
};

template <typename _FullGeodesicData>
class SnapshotRecursive : public Snapshot<_FullGeodesicData> {
 public:
  SnapshotRecursive(int width, int height)
      : Snapshot<_FullGeodesicData>(width, height) {}
  virtual bool load(FILE *f) {
    unsigned int n;
    if (fread(&n, sizeof(unsigned int), 1, f) != 1)
      return false;
    geodesics.resize(n);
    if (fread(&geodesics[0], sizeof(geodesics[0]), n, f) != n)
      return false;
    if (fread(&n, sizeof(unsigned int), 1, f) != 1)
      return false;
    leafs.resize(n);
    if (fread(&leafs[0], sizeof(leafs[0]), n, f) != n)
      return false;
    return true;
  }

  virtual bool save(FILE *f) const {
    unsigned int n = (unsigned int)geodesics.size();
    if (fwrite(&n, sizeof(unsigned int), 1, f) != 1)
      return false;
    if (fwrite(&geodesics[0], sizeof(geodesics[0]), n, f) != n)
      return false;
    n = (unsigned int)leafs.size();
    if (fwrite(&n, sizeof(unsigned int), 1, f) != 1)
      return false;
    if (fwrite(&leafs[0], sizeof(leafs[0]), n, f) != n)
      return false;
    return true;
  }

  std::vector<ScreenGeodesic<_FullGeodesicData>> geodesics;
  std::vector<QuadtreeLeaf> leafs;
};

template <typename _FullGeodesicData>
using task_stack_t = std::stack<std::pair<
    typename std::vector<ScreenGeodesic<_FullGeodesicData>>::iterator,
    typename std::vector<ScreenGeodesic<_FullGeodesicData>>::iterator>>;

template <typename _FullGeodesicData,
          typename _Spacetime,
          typename _Field,
          typename DLambdaFunc>
void _integrate_geodesics_worker_thread(
    task_stack_t<_FullGeodesicData> &tasks,
    std::mutex &tasks_mutex,
    const _Spacetime &spacetime,
    const _Field &field,
    const Camera *raytracer_camera,
    const DLambdaFunc &dlambda_func,
    int width,
    int height) {

  for (;;) {
    tasks_mutex.lock();
    if (tasks.empty()) {
      tasks_mutex.unlock();
      break;
    }
    typename std::vector<ScreenGeodesic<_FullGeodesicData>>::iterator begin, end;
    std::tie(begin, end) = tasks.top();
    tasks.pop();
    tasks_mutex.unlock();

    for (; begin != end; ++begin) {
      int dead_reason = integrate_single_geodesic(
          spacetime,
          field,
          raytracer_camera,
          dlambda_func,
          (double)begin->x / LOCATION_FACTOR / width,
          (double)begin->y / LOCATION_FACTOR / height,
          &begin->full);
      if (dead_reason == DEAD_FLAT) {
        typedef decltype(begin->full.basic.position) _Coord;
        typedef typename _Coord::value_type value_type;
        CartesianVector4<value_type> position, direction;
        convert_point_and_diff(
            spacetime.coord_system_parameters(begin->full.basic.position),
            begin->full.basic.position,
            begin->full.basic.direction,
            Null(),
            position,
            direction);
        Vector<value_type, 3> pos3{position[1], position[2], position[3]};
        Vector<value_type, 3> dir3{direction[1], direction[2], direction[3]};
        auto cos_theta = dir3.dot(pos3) / (pos3.length() * dir3.length());
        auto eta = (MAX_r - pos3.length()) / cos_theta / dir3.length();
        mult_add(&begin->full.basic.position, eta, begin->full.basic.direction);
      }
    }
  }
}


template <typename _FullGeodesicData, typename _Spacetime, typename _Field,
         typename _dlambdaFunc>
class ImageRecursiveGenerator {
  typedef typename _FullGeodesicData::basic_type _BasicGeodesicState;
  typedef typename _FullGeodesicData::extra_type _GeodesicExtra;
  typedef ScreenGeodesic<_FullGeodesicData> _ScreenGeodesic;

  const _Spacetime &spacetime;
  const _Field &field;
  const Camera *raytracer_camera;
  const _dlambdaFunc &dlambda_func;
  const int width;
  const int height;
  SnapshotRecursive<_FullGeodesicData> * const output_snapshot;
  const int thread_count;
  const int max_extra_depth;

  std::vector<_ScreenGeodesic> geodesics;
  std::map<std::pair<int, int>, int> geodesics_map;
  std::vector<QuadtreeNode> nodes;
  std::map<std::pair<int, int>, int> nodes_map;
  std::vector<int> leafs;
  std::vector<int> new_leafs;

  size_t last_geodesic;

  const int min_size;
  double split_eps;

 public:
  ImageRecursiveGenerator(
      const _Spacetime &spacetime,
      const _Field &field,
      const Camera *raytracer_camera,
      const _dlambdaFunc &dlambda_func,
      int width,
      int height,
      SnapshotRecursive<_FullGeodesicData> *output_snapshot,
      int thread_count,
      int max_extra_depth) : spacetime(spacetime), field(field),
          raytracer_camera(raytracer_camera),
          dlambda_func(dlambda_func),
          width(width), height(height), output_snapshot(output_snapshot),
          thread_count(thread_count), max_extra_depth(max_extra_depth),
          min_size(LOCATION_FACTOR >> max_extra_depth) {}


  void generate(void) {
    constexpr int BLOCK = 64;
    const int hor = width / BLOCK + 2;  // Number of points (not regions).
    const int ver = height / BLOCK + 2;

    split_eps = 5.;

    last_geodesic = 0;

    geodesics.reserve(hor * ver);
    for (int i = 0; i < ver; ++i)
      for (int j = 0; j < hor; ++j) {
        // These geodesics do not point at the center of pixels.
        // double x = double(width - hor * BLOCK + 2 * j * BLOCK) / (2 * width);
        // double y = double(height - ver * BLOCK + 2 * i * BLOCK) / (2 * height);
        // const real_t x = real_t((2 * j - hor + 1) * BLOCK) / (2 * width);
        // const real_t y = real_t((2 * i - ver + 1) * BLOCK) / (2 * height);
        // We use integers for the screen location and divide by (2 * width) or
        // (2 * height) when passing coordinates to integrate_single_geodesic,
        // and also when rendering the image.
        int x = LOCATION_FACTOR * (2 * j - hor + 1) * BLOCK;
        int y = LOCATION_FACTOR * (2 * i - ver + 1) * BLOCK;;

        _add_geodesic_to_queue(x, y);
      }
    _integrate_geodesics();

    nodes.reserve((hor - 1) * (ver - 1));
    for (int i = 1; i < ver; ++i)
      for (int j = 1; j < hor; ++j) {
        /* 0 1
         * 3 2 */
        int pivot = (i - 1) * hor + j - 1;
        int x = (geodesics[pivot].x + geodesics[pivot + hor + 1].x) / 2;
        int y = (geodesics[pivot].y + geodesics[pivot + hor + 1].y) / 2;
        nodes_map[std::make_pair(x, y)] = (int)nodes.size();
        nodes.push_back(QuadtreeNode{
            {pivot, pivot + 1, pivot + hor + 1, pivot + hor},
            x, y,
            BLOCK * LOCATION_FACTOR,
            -1,    // No parent.
            false  // Not split.
        });
      }

    leafs.resize(nodes.size());
    for (int i = 0; i < (int)leafs.size(); ++i)
      leafs[i] = i;

    fprintf(stderr, "Leafs left...");
    while (!leafs.empty()) {
      fprintf(stderr, " %d", (int)leafs.size());
      new_leafs.clear();
      for (size_t k = 0; k < leafs.size(); ++k) {
        int id = leafs[k];
        if (!_should_leaf_split(nodes[id]))
          continue;

        // Be careful, do not use references to QuadtreeNode, as this might
        // invalidate them.
        _split_node(leafs[k]);

      }
      _integrate_geodesics();
      leafs.swap(new_leafs);
    }
    fprintf(stderr, "  Done!\n");

    output_snapshot->leafs.clear();
    for (const QuadtreeNode &node : nodes)
      if (!node.split)
        output_snapshot->leafs.push_back(QuadtreeLeaf(node));
    output_snapshot->geodesics.swap(geodesics);
  }

 private:

  inline int _add_geodesic_to_queue(int x, int y) {
    geodesics.push_back(_ScreenGeodesic{
        _FullGeodesicData(),
        x, y
      });
    geodesics_map.emplace(std::make_pair(x, y), (int)geodesics.size() - 1);
    return (int)geodesics.size() - 1;
  }

  void _integrate_geodesics(void) {
    std::vector<std::thread> threads;

    task_stack_t<_FullGeodesicData> tasks;
    std::mutex tasks_mutex;

    constexpr int SPLIT = 10;
    size_t diff = geodesics.size() - last_geodesic;
    for (size_t i = 0; i < (diff + SPLIT - 1) / SPLIT; ++i) {
      size_t from = i * SPLIT;
      size_t to = std::min(diff, (i + 1) * SPLIT);
      tasks.emplace(
          geodesics.begin() + last_geodesic + from,
          geodesics.begin() + last_geodesic + to
      );
    }

    for (int i = 0; i < thread_count; ++i) {
      threads.emplace_back(
          _integrate_geodesics_worker_thread<
              _FullGeodesicData, _Spacetime, _Field, _dlambdaFunc>,
          std::ref(tasks),
          std::ref(tasks_mutex),
          std::ref(spacetime),
          std::ref(field),
          raytracer_camera,
          std::ref(dlambda_func),
          width,
          height
      );
    }
    for (auto &thread : threads) {
      if (thread.joinable())
        thread.join();
    }
    last_geodesic = geodesics.size();
  }

  int _split_edge(int a, int b) {
    int x = (geodesics[a].x + geodesics[b].x) / 2;
    int y = (geodesics[a].y + geodesics[b].y) / 2;
    auto it = geodesics_map.find(std::make_pair(x, y));
    if (it != geodesics_map.end())
      return it->second;

    return _add_geodesic_to_queue(x, y);
  }

  bool _should_leaf_split(const QuadtreeNode &leaf) {
    if (leaf.size <= min_size)
      return false;

    // Did they hit the same target?
    int reason_mask = 0;
    int sky_tex_dead_reason = -1;
    bool too_far = true;
    for (int i = 0; i < 4; ++i) {
      const auto &geo = geodesics[leaf.geodesic_ids[i]];
      const int dead_reason = geo.full.extra.dead_reason;
      if (dead_reason >= DEAD_SKY_TEX_OFFSET) {
        if (sky_tex_dead_reason >= 0) {
          if (sky_tex_dead_reason != dead_reason)
            sky_tex_dead_reason = -2;
        } else if (sky_tex_dead_reason == -1) {
          sky_tex_dead_reason = dead_reason;
        }
        reason_mask |= 1 << DEAD_UNUSED;
      } else {
        reason_mask |= 1 << dead_reason;
        if (geo.full.extra.min_r < spacetime.black_hole_radius() * 3)
          too_far = false;
      }
    }

    if (reason_mask & (1 << DEAD_BLACK_HOLE)) {
      // If only black hole, do not split. Otherwise, split.
      return reason_mask == (1 << DEAD_BLACK_HOLE) ? false : true;
    }

    if (sky_tex_dead_reason == -2
        && leaf.size > (LOCATION_FACTOR >> std::max(0, max_extra_depth - 3)))
      return true;
    if (__builtin_popcount(reason_mask) > 1
        && leaf.size > (LOCATION_FACTOR >> std::min(1, max_extra_depth)))
      return true;  // Divide at least once the outer edge of the disk.
    if (too_far
        && leaf.size <= (LOCATION_FACTOR >> std::max(0, max_extra_depth - 3)))
      return false;
    if (reason_mask == (1 << DEAD_DISK) && leaf.size > 4 * LOCATION_FACTOR)
      return true;
    if (__builtin_popcount(reason_mask) > 1)
      return true;

    for (int i = 0; i < 4; ++i) {
      const _ScreenGeodesic &A = geodesics[leaf.geodesic_ids[i]];
      // if (A.dead_reason == DEAD_BLACK_HOLE) return true;
      const _ScreenGeodesic &B = geodesics[leaf.geodesic_ids[(i + 1) % 4]];
      double eps =
          numerical_sqr_distance(A.full.basic.position, B.full.basic.position) +
          numerical_sqr_distance(A.full.basic.direction, B.full.basic.direction);
      if (eps > split_eps)
        return true;
    }
    return false;
  }

  void _split_node(int id) {
    // Copy node as the reference might be invalidated later.
    if (nodes[id].split)
      return;
    nodes[id].split = true;
    QuadtreeNode node = nodes[id];

#define CID(x) node.geodesic_ids[x]
    /* Geodesic indexing:
     * 0 a 1
     * d e b
     * 3 c 2
     */
    const int a = _split_edge(CID(0), CID(1));
    const int b = _split_edge(CID(1), CID(2));
    const int c = _split_edge(CID(2), CID(3));
    const int d = _split_edge(CID(3), CID(0));
    const int e = _split_edge(a, c);

#define ADD_NODE(_a, _b, _c, _d, _geo) {  \
      const int cx = (_geo.x + node.center_x) / 2;  \
      const int cy = (_geo.y + node.center_y) / 2;  \
      new_leafs.push_back((int)nodes.size());  \
      nodes_map[std::make_pair(cx, cy)] = (int)nodes.size();  \
      nodes.push_back(QuadtreeNode{  \
        {_a, _b, _c, _d}, cx, cy, node.size / 2, (int)id, false  \
      });  \
    }

    ADD_NODE(CID(0), a, e, d, geodesics[CID(0)]);
    ADD_NODE(a, CID(1), b, e, geodesics[CID(1)]);
    ADD_NODE(e, b, CID(2), c, geodesics[CID(2)]);
    ADD_NODE(d, e, c, CID(3), geodesics[CID(3)]);
#undef ADD_NODE
#undef CID

    // Split any neighbouring nodes to at least my depth.
    if (node.parent == -1)
      return;  // OK, nothing to split further.

    // Two of the neighbours are already at my level, split only those
    // opposite of the parent's center.
    const int cx = node.center_x;
    const int cy = node.center_y;
    const int dx = cx - nodes[node.parent].center_x;
    const int dy = cy - nodes[node.parent].center_y;

    auto it = nodes_map.find(std::make_pair(cx + 3 * dx, cy - dy));
    if (it != nodes_map.end())
      _split_node(it->second);
    it = nodes_map.find(std::make_pair(cx - dx, cy + 3 * dy));
    if (it != nodes_map.end())
      _split_node(it->second);
  }

};

template <typename _FullGeodesicData, typename _Spacetime, typename _Field,
          typename _dlambdaFunc>
inline void generate_image_recursive(
    const _Spacetime &spacetime,
    const _Field &field,
    const Camera *raytracer_camera,
    const _dlambdaFunc &dlambda_func,
    int width,
    int height,
    SnapshotRecursive<_FullGeodesicData> *output_snapshot,
    int thread_count,
    int max_extra_depth) {
  ImageRecursiveGenerator<_FullGeodesicData, _Spacetime, _Field,_dlambdaFunc>
      generator(
        spacetime, field, raytracer_camera, dlambda_func,
        width, height, output_snapshot, thread_count, max_extra_depth);
  generator.generate();
}


template <typename _Spacetime, typename _FullGeodesicData, typename DiskTex>
void colorize_from_recursive_snapshot(
    const SnapshotRecursive<_FullGeodesicData> &snapshot,
    const _Spacetime &spacetime,
#if SKY_ENABLED
    const Image &/* sky */,
#endif
    const DiskTex &disk_tex,
    RGBd *output) {

  const int width = snapshot.width;
  const int height = snapshot.height;
  const auto &geodesics = snapshot.geodesics;

  for (int k = 0; k < width * height; ++k)
    output[k] = RGBd{0, 0, 0};

  auto to_ij = [width, height](int x, int y) {
    // We introduced integer coordinates to make splitting simpler.
    // Coordinates x, y are now multiplied by 2 * width, 2 * height.
    double j = .5 * (((double)x / LOCATION_FACTOR + width) - 1);
    double i = .5 * (((double)y / LOCATION_FACTOR + height) - 1);
    return std::make_pair(j, i);
  };

  int max_steps = 0;
  for (int k = 0; k < (int)geodesics.size(); ++k)
    max_steps = std::max(max_steps, geodesics[k].full.extra.steps);
  fprintf(stderr, "MAX STEPS: %d\n", max_steps);

  for (int k = 0; k < (int)snapshot.leafs.size(); ++k) {
    const auto &leaf = snapshot.leafs[k];
    const int * const ids = leaf.geodesic_ids;

    const auto &G0 = geodesics[ids[0]];
    const auto &G1 = geodesics[ids[1]];
    const auto &G2 = geodesics[ids[2]];
    const auto &G3 = geodesics[ids[3]];

    double i0, i1, i2, i3;
    double j0, j1, j2, j3;
    std::tie(j0, i0) = to_ij(G0.x, G0.y);
    std::tie(j1, i1) = to_ij(G1.x, G1.y);
    std::tie(j2, i2) = to_ij(G2.x, G2.y);
    std::tie(j3, i3) = to_ij(G3.x, G3.y);

    if (G0.full.extra.dead_reason != G1.full.extra.dead_reason
        || G0.full.extra.dead_reason != G2.full.extra.dead_reason
        || G0.full.extra.dead_reason != G3.full.extra.dead_reason
// #if FAKE_SKY
        || G0.full.extra.dead_reason >= DEAD_SKY_TEX_OFFSET
// #endif
        ) {

#if SKY_ENABLED
      RGBd a = get_single_geodesic_color(spacetime, sky, disk_tex, G0.full);
      RGBd b = get_single_geodesic_color(spacetime, sky, disk_tex, G1.full);
      RGBd c = get_single_geodesic_color(spacetime, sky, disk_tex, G2.full);
      RGBd d = get_single_geodesic_color(spacetime, sky, disk_tex, G3.full);
#else
      RGBd a = get_single_geodesic_color(spacetime, disk_tex, G0.full);
      RGBd b = get_single_geodesic_color(spacetime, disk_tex, G1.full);
      RGBd c = get_single_geodesic_color(spacetime, disk_tex, G2.full);
      RGBd d = get_single_geodesic_color(spacetime, disk_tex, G3.full);
#endif
      render_textured_axes_aligned_rectangle(
          j0, i0, j2, i2,
          a, b, c, d,
          ColorMixTexture(),
          width, height, output);
    } else if (G0.full.extra.dead_reason == DEAD_OBJECT) {
      render_axes_aligned_rectangle(
          j0, i0, j2, i2, OBJECT_COLOR, width, height, output);
    } else if (G0.full.extra.dead_reason == DEAD_BLACK_HOLE) {
      render_axes_aligned_rectangle(
          j0, i0, j2, i2, BLACK_HOLE_COLOR, width, height, output);
#if SKY_ENABLED
#error Not implemented
#elif !FAKE_SKY
    } else if (G0.full.extra.dead_reason == DEAD_FLAT) {
      render_axes_aligned_rectangle(
          j0, i0, j2, i2, SKY_COLOR, width, height, output);
#endif
#if RENDER_DISK
    } else if (G0.full.extra.dead_reason == DEAD_DISK) {
      render_textured_axes_aligned_rectangle(
          j0, i0, j2, i2,
          disk_tex.get_tex_coord(spacetime, G0.full),
          disk_tex.get_tex_coord(spacetime, G1.full),
          disk_tex.get_tex_coord(spacetime, G2.full),
          disk_tex.get_tex_coord(spacetime, G3.full),
          disk_tex, width, height, output);
#endif
    } else {
      fprintf(stderr, "Dead reason %d??", G0.full.extra.dead_reason);
      // Mark with a really bright pink color.
      render_axes_aligned_rectangle(
          j0, i0, j2, i2, RGBd{1e9, 0, 1e9}, width, height, output);
    }

    // RGBd color = RGBd{1, 1, 1};
    // RGBd color = RGBd{-1, -1, -1};
    // render_axes_aligned_rectangle_border(j0, i0, j2, i2, color, 0.3, 0.3, width, height, output);
  }

  fprintf(stderr, "total geodesics: %d\n", (int)snapshot.geodesics.size());
}

#endif
