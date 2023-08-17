#ifndef _MT_MESH_H_
#define _MT_MESH_H_ 1

#include "mt_defs.h"
#include "nanoflann.hpp"

#include <deal.II/fe/mapping.h>

#include <memory>

class TriaFinder {
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, TriaFinder>,
                                              TriaFinder, 2>
      KDTree;

public:
  size_t kdtree_get_point_count() const { return mesh->n_vertices(); }

  double kdtree_distance(const double *p, const size_t idx, size_t) const {
    Point pt;
    double distance;

    if (mesh->vertex_used(idx)) {
      pt = mesh->get_vertices()[idx];
      distance = std::sqrt((p[0] - pt[0]) * (p[0] - pt[0]) + (p[1] - pt[1]) * (p[1] - pt[1]));
    } else {
      distance = std::numeric_limits<double>::max();
    }

    return distance;
  }

  double kdtree_get_pt(const size_t idx, int dim) const { return mesh->get_vertices()[idx][dim]; }

  template <typename BBOX> bool kdtree_get_bbox(BBOX &) const { return false; }

public:
  void reinit(const Triangulation &);
  int find_closest_vertex(const Point &p) const;
  std::pair<Triangulation::active_cell_iterator, Point>
  find_active_cell_around_point(const Point &p) const;

private:
  const Triangulation *mesh;
  std::shared_ptr<KDTree> v_index;
  std::shared_ptr<dealii::Mapping<2>> mapping;
  std::vector<Triangulation::active_cell_iterator> v_to_c;
};

#endif
