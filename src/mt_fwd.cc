#include "mt_ctx.h"
#include "mt_io.h"
#include "mt_pw.h"
#include "reEstimate.h"
#include"mt_refine.h"

#include <Eigen/SparseLU>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_tools.h>

 
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_refinement.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/refinement.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/smoothness_estimator.h>
#include <deal.II/numerics/error_estimator.h>

#include <stdio.h>
#include <fstream>
#include <mutex>
#include <map>
#include <time.h>

void get_vertical_layer(MT2DCtx *ctx, int c) {
  Point p, fp;
  int pidx, i;
  double sigma;
  std::vector<Layer> layers;
  Triangulation::active_cell_iterator cell;

  p[1] = ctx->bbox.first[1];
  switch (c) {
  case Left:
    p[0] = ctx->bbox.first[0];
    break;
  case Right:
    p[0] = ctx->bbox.second[0];
    break;
  default:
    break;
  }

  fp = p;
  while (fp[1] < ctx->bbox.second[1]) {
    pidx = ctx->tria_finder->find_closest_vertex(Point(fp[0], fp[1] + 1.0));
    cell = ctx->tria_finder->find_active_cell_around_point(Point(fp[0], fp[1] + 1.0)).first;

    if (cell.state() != dealii::IteratorState::valid) {
    }

    sigma = 1.0 / ctx->rho[cell->user_index()];
    layers.push_back(std::make_pair(fp[1], sigma));

    for (i = 0; i < (int)dealii::GeometryInfo<2>::vertices_per_cell; ++i) {
      if ((int)cell->vertex_index(i) != pidx && std::abs(cell->vertex(i)[0] - p[0]) < EPS) {
        fp = cell->vertex(i);
        break;
      }
    }
  }
  ctx->layered_models[c] = layers;
}

void calculate_layered_models(MT2DCtx *ctx) {
  int c;
  Point p;

  p[0] = p[1] = -1.0E08;
  ctx->bbox.first = ctx->coarse_mesh->get_vertices()[ctx->tria_finder->find_closest_vertex(p)];
  p[0] = p[1] = 1.0E08;
  ctx->bbox.second = ctx->coarse_mesh->get_vertices()[ctx->tria_finder->find_closest_vertex(p)];

  for (c = 0; c < 2; ++c) {
    get_vertical_layer(ctx, c);
  }
}

void create_delta_vector(MT2DCtx *ctx, FEMData *fd, VectorD &delta)
{
  VectorD per_rx_rhs(fd->dh.n_dofs());
  delta = 0.0;
  for (unsigned int i = 0; i < ctx->rx.size(); i++)
  {
    per_rx_rhs = 0.0;

    dealii::VectorTools::create_point_source_vector(fd->dh, ctx->rx[i], per_rx_rhs);

    delta += per_rx_rhs;
  } 
}

void setup_system(MT2DCtx *ctx, FEMData *fd, int fidx, int mode) {
  BoundaryFunction bf;

  fd->dh.reinit(fd->mesh);
  fd->dh.distribute_dofs(*ctx->fe);
  bf.reinit(ctx->bbox, ctx->layered_models, ctx->freqs[fidx], mode);

  fd->ac.clear();
  dealii::DoFTools::make_hanging_node_constraints(fd->dh, fd->ac);
  dealii::VectorTools::interpolate_boundary_values(fd->dh, 0, bf, fd->ac);
  fd->ac.close();

  dealii::DynamicSparsityPattern dsp(fd->dh.n_dofs(), fd->dh.n_dofs());
  dealii::DoFTools::make_sparsity_pattern(fd->dh, dsp, fd->ac, false);
  fd->sp.copy_from(dsp);

  fd->K.reinit(fd->sp);
  fd->u.reinit(fd->dh.n_dofs());
  fd->s.reinit(fd->dh.n_dofs());
  fd->delta_s.reinit(fd->dh.n_dofs());
  fd->delta.reinit(fd->dh.n_dofs());
  create_delta_vector(ctx, fd, fd->delta);
}

void assemble_linear_system(MT2DCtx *ctx, FEMData *fd, int fidx, int mode) {
  dealii::QGauss<2> quadrature_formula(5);
  dealii::FEValues<2> fe_values(*ctx->fe, quadrature_formula,
                                dealii::update_values | dealii::update_gradients |
                                    dealii::update_quadrature_points | dealii::update_JxW_values);
  const size_t dofs_per_cell = ctx->fe->dofs_per_cell;
  const size_t n_q_points = quadrature_formula.size();

  dealii::FullMatrix<Complex> cell_matrix(dofs_per_cell, dofs_per_cell);
  dealii::Vector<Complex> cell_rhs(dofs_per_cell);

  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

  size_t i, j, q;
  double omega, sigma;
  Complex coeff_a, coeff_b;

  omega = 2 * dealii::numbers::PI * ctx->freqs[fidx];

  fd->K = 0.0;
  fd->s = 0.0;
  for (auto cell = fd->dh.begin_active(); cell != fd->dh.end(); ++cell) {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;
    sigma = 1.0 / ctx->rho[cell->user_index()];

    if (mode == TE) {
      coeff_a = 1.0;
      coeff_b = II * omega * MU * sigma;
    } else {
      coeff_a = 1.0 / sigma;
      coeff_b = II * omega * MU;
    }

    for (i = 0; i < dofs_per_cell; ++i) {
      for (j = 0; j < dofs_per_cell; ++j) {
        for (q = 0; q < n_q_points; ++q) {
          cell_matrix(i, j) += (coeff_a * fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)
              + coeff_b * fe_values.shape_value(i, q) * fe_values.shape_value(j, q)) * fe_values.JxW(q);
        }
      }
    }
    cell->get_dof_indices(local_dof_indices);

    fd->ac.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, fd->K, fd->s);
  }
}

void solve_linear_system(FEMData *fd) {
  int i;
  Eigen::VectorXcd x1(fd->K.n()), b1(fd->K.n());
  Eigen::SparseMatrix<Complex, Eigen::ColMajor, int> A(fd->K.m(), fd->K.n());
  std::vector<Eigen::Triplet<Complex, int> > triplet_list(fd->K.n_nonzero_elements());
  Eigen::SparseLU<Eigen::SparseMatrix<Complex, Eigen::ColMajor, int>, Eigen::COLAMDOrdering<int> > solver;

  for (SparseMatrixZ::const_iterator it = fd->K.begin(); it != fd->K.end(); ++it) {
    triplet_list.push_back(Eigen::Triplet<Complex, int>(it->row(), it->column(), it->value()));
  }
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());
  A.makeCompressed();

  solver.analyzePattern(A);
  solver.factorize(A);

  for (i = 0; i < (int)fd->s.size(); ++i) {
    b1[i] = fd->s[i];
  }
  x1 = solver.solve(b1);
  for (i = 0; i < (int)fd->s.size(); ++i) {
    fd->u[i] = x1[i];
  }

  fd->ac.distribute(fd->u);
}

void solve_delta_fun(FEMData *fd)
{
  int i;
  Eigen::VectorXcd x2(fd->K.n()), b2(fd->K.n());
  Eigen::SparseMatrix<Complex, Eigen::ColMajor, int> A(fd->K.m(), fd->K.n());
  std::vector<Eigen::Triplet<Complex, int> > triplet_list(fd->K.n_nonzero_elements());
  Eigen::SparseLU<Eigen::SparseMatrix<Complex, Eigen::ColMajor, int>, Eigen::COLAMDOrdering<int> > solver;

  for (SparseMatrixZ::const_iterator it = fd->K.begin(); it != fd->K.end(); ++it) {
    triplet_list.push_back(Eigen::Triplet<Complex, int>(it->row(), it->column(), it->value()));
  }
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());
  A.makeCompressed();

  solver.analyzePattern(A);
  solver.factorize(A);

  for (i = 0; i < (int)fd->delta.size(); ++i) {
    b2[i] = fd->delta[i];
  }
  x2 = solver.solve(b2);
  for (i = 0; i < (int)fd->delta.size(); ++i) {
    fd->delta_s[i] = x2[i];
  }
  fd->ac.distribute(fd->delta_s);
}

void compute_interpolation_vector(
    const std::pair<DoFHandler::active_cell_iterator, Point> &cell_point,
    std::vector<dealii::types::global_dof_index> &dof_indices, std::vector<double> dof_coeffs[2]) {
  size_t i, j, dofs_per_cell;
  const dealii::Quadrature<2> quadrature(cell_point.second);
  dealii::FEValues<2> fe_values(dealii::StaticMappingQ1<2>::mapping, cell_point.first->get_fe(),
                                quadrature, dealii::update_values | dealii::update_gradients);

  fe_values.reinit(cell_point.first);

  dofs_per_cell = fe_values.get_fe().dofs_per_cell;

  dof_indices.resize(dofs_per_cell);
  for (j = 0; j < 2; ++j) {
    dof_coeffs[j].resize(dofs_per_cell);
  }

  cell_point.first->get_dof_indices(dof_indices);

  for (i = 0; i < dofs_per_cell; ++i) {
    dof_coeffs[0][i] = fe_values.shape_value_component(i, 0, 0);
    dof_coeffs[1][i] = fe_values.shape_grad_component(i, 0, 0)[1];
  }
}

void interpolate_fields(std::vector<dealii::types::global_dof_index> &dof_indices,
                        std::vector<double> dof_coeffs[2], const VectorZ &u, Complex &e,
                        Complex &h) {
  size_t i;

  e = h = 0.0;
  for (i = 0; i < dof_indices.size(); ++i) {
    e += u[dof_indices[i]] * dof_coeffs[0][i];
    h += u[dof_indices[i]] * dof_coeffs[1][i];
  }
}

void calculate_response_mt(MT2DCtx *ctx, FEMData *fd, int fidx, int mode) {
  int ridx, o;
  double omega;
  Complex e, h;
  std::vector<double> dof_coeffs[2];
  std::vector<dealii::types::global_dof_index> dof_indices;
  std::pair<Triangulation::active_cell_iterator, Point> cell_point;

  omega = 2 * dealii::numbers::PI * ctx->freqs[fidx];

  for (ridx = 0; ridx < (int)ctx->rx.size(); ++ridx) {
    cell_point = fd->tf.find_active_cell_around_point(ctx->rx[ridx]);

    if (cell_point.first.state() != dealii::IteratorState::valid) {
    }

    compute_interpolation_vector(
        std::make_pair(DoFHandler::active_cell_iterator(*(cell_point.first), &fd->dh),
                       cell_point.second),
        dof_indices, dof_coeffs);

    if (mode == TE) {
      interpolate_fields(dof_indices, dof_coeffs, fd->u, e, h);
      h /= (II * omega * MU);
    } else {
      interpolate_fields(dof_indices, dof_coeffs, fd->u, h, e);
      e *= ctx->rho[cell_point.first->user_index()];
    }

    for (o = 0; o < (int)ctx->rsp.size(); ++o) {
      if (ctx->fidx[o] != fidx || TYPE_TO_MODE.at(ctx->otype[o]) != mode || ctx->ridx[o] != ridx) {
        continue;
      }

      switch (ctx->otype[o]) {
      case R_XY_AP:
      case R_YX_AP:
        ctx->rsp[o] = std::log((e / h) * (e / h) / (omega * MU));
        break;
        break;
      default:
        break;
      }
    }
  }
}

void refine_receiving_area(MT2DCtx *ctx, FEMData *fd, int fidx) {
  int r, i;
  std::set<int> rx_indices;
  std::set<int>::iterator it;
  Triangulation::active_cell_iterator cell;

  for (i = 0; i < (int)ctx->obs.size(); ++i) {
    if (ctx->fidx[i] == fidx) {
      rx_indices.insert(ctx->ridx[i]);
    }
  }

  for (r = 0; r < ctx->n_rx_cell_refinements; ++r) {
    fd->tf.reinit(fd->mesh);
    for (it = rx_indices.begin(); it != rx_indices.end(); ++it) {
      cell = fd->tf.find_active_cell_around_point(ctx->rx[*it]).first;

      if (cell.state() != dealii::IteratorState::valid) {
      }

      if (cell->diameter() > ctx->max_rx_cell_size) {
        cell->set_refine_flag();
      }
    }
    fd->mesh.execute_coarsening_and_refinement();
  }
}

void refine_mesh(MT2DCtx *ctx, FEMData *fd, int fidx, int index)
{
  switch (ctx->refine_mode)
  {
  case 0:
    refine_goal(ctx, fd,fidx, index);
    break;
  case 1:
    refine_adaptive(ctx, fd, fidx, index);
    break;
  case 2:
    refine_global(ctx, fd, fidx);
    break;
  default:
    break;
  }
  
}

void mt_forward_per_freq_mode(MT2DCtx *ctx, int fidx, int mode) {
  FEMData fd;
 
  for(int i = 0; i <= ctx->max_adaptive_refinements+1; i++)
  {
    clock_t t_now;
    t_now = clock();

    if(i == 0)
    {
      fd.mesh.copy_triangulation(*ctx->coarse_mesh);
    }
    else
    {
      refine_mesh(ctx, &fd, fidx, i);
      if(i == ctx->max_adaptive_refinements+1)
      {
        continue;
      }  
    }
    fd.tf.reinit(fd.mesh);
    std::cout << "   Number of active cells:       "
                << fd.mesh.n_active_cells() << std::endl;

    setup_system(ctx, &fd, fidx, mode);
    std::cout << "   Number of degrees of freedom: " << fd.dh.n_dofs()
                << std::endl;

    assemble_linear_system(ctx, &fd, fidx, mode);

    solve_linear_system(&fd);
   
    solve_delta_fun(&fd);

    calculate_response_mt(ctx, &fd, fidx, mode);

    clock_t t = clock() - t_now;

    std::cout<< "Time diff:" << float(t)/CLOCKS_PER_SEC << "s" << std::endl;

    if(ctx->refine_mode == 0)
    {
      save_rsp(ctx, (ctx->goal_rsp_oprefix+ "-adap" + std::to_string(fidx) + "-" + std::to_string(i) + ".rsp").c_str());
    }
    else if(ctx->refine_mode == 1)
    {
      save_rsp(ctx, (ctx->adaptive_rsp_oprefix+ "-adap" + std::to_string(fidx) + "-" + std::to_string(i) + ".rsp").c_str());
    }
  }

  
}

void mt_forward(MT2DCtx *ctx) {
  int n_finished = 0;
  std::set<std::pair<int, int> > tasks;

  for (int i = 0; i < (int)ctx->obs.size(); ++i) {
    tasks.insert(std::make_pair(ctx->fidx[i], TYPE_TO_MODE.at(ctx->otype[i])));
  }

  calculate_layered_models(ctx);

  for (const std::pair<int, int>& task : tasks) {
    mt_forward_per_freq_mode(ctx, task.first, task.second);
    printf("Process: %3d of %3d finished.\n", ++n_finished, (int)tasks.size());
  }

 
}
