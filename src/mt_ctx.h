#ifndef _MT_CTX_H_
#define _MT_CTX_H_

#include "mt_defs.h"
#include "mt_mesh.h"

#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/base/quadrature_lib.h>



#include <memory>
#include <string>
#include <vector>

struct MT2DCtx {
  std::string iprefix, oprefix, adaptive_rsp_oprefix, goal_rsp_oprefix, adaptive_vtk_oprefix, goal_vtk_oprefix;

  std::vector<Point> rx;
  std::vector<double> freqs;
  std::vector<int> fidx, ridx, otype;
  std::vector<Complex> obs, rsp, obserr;

  VectorD rho, rho_lb, rho_ub;

  std::shared_ptr<FE_Q> fe;
  std::shared_ptr<Triangulation> coarse_mesh;
  std::shared_ptr<TriaFinder> tria_finder;

  BoundingBox bbox;
  std::vector<Layer> layered_models[2];

  dealii::ParameterHandler prm;
  // double refine_fraction, max_rx_cell_size;
  double max_rx_cell_size, refine_threshold;
  // int min_degree, max_degree;
  int max_adaptive_refinements, max_dofs, n_global_refinements, n_rx_cell_refinements;
  int refine_mode;
};

struct FEMData {
  TriaFinder tf;
  Triangulation mesh;
  DoFHandler dh;
  AffineConstraintsZ ac;
  dealii::SparsityPattern sp;
  SparseMatrixZ K;
  VectorZ u, s, delta_s;
  VectorD  delta;
};

void declare_parameter(MT2DCtx *);
void read_parameters(MT2DCtx *, const std::string &);
void create_ctx(MT2DCtx *);
void destroy_ctx(MT2DCtx *);

#endif
