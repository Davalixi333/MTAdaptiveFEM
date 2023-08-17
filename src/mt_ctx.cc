#include "mt_ctx.h"
void declare_parameter(MT2DCtx *ctx) {
  ctx->prm.declare_entry("iprefix", "", dealii::Patterns::FileName(), "Prefix of the input files");
  ctx->prm.declare_entry("oprefix", "",
                         dealii::Patterns::FileName(dealii::Patterns::FileName::output),
                         "Prefix of the output files");
  ctx->prm.declare_entry("adaptive_rsp_oprefix", "",
                         dealii::Patterns::FileName(dealii::Patterns::FileName::output),
                         "Prefix of the output files");
  ctx->prm.declare_entry("goal_rsp_oprefix", "",
                         dealii::Patterns::FileName(dealii::Patterns::FileName::output),
                         "Prefix of the output files");
  ctx->prm.declare_entry("goal_vtk_oprefix", "",
                         dealii::Patterns::FileName(dealii::Patterns::FileName::output),
                         "Prefix of the output files");
  ctx->prm.declare_entry("adaptive_vtk_oprefix", "",
                         dealii::Patterns::FileName(dealii::Patterns::FileName::output),
                         "Prefix of the output files");
  ctx->prm.declare_entry("max_adaptive_refinements", "1", dealii::Patterns::Integer(0, 20), "Maximum adaptive refinements");
  // ctx->prm.declare_entry("refine_fraction", "0.3", dealii::Patterns::Double(0.0, 1.0), "Refine fraction");
  ctx->prm.declare_entry("n_global_refinements", "0", dealii::Patterns::Integer(0, 5), "Maximum global refinements");
  ctx->prm.declare_entry("n_rx_cell_refinements", "0", dealii::Patterns::Integer(0, 10), "Number of refinements at receiver point");
  ctx->prm.declare_entry("max_rx_cell_size", "-1.0", dealii::Patterns::Double(-1.0, 100), "Maximum cell size at receiver point");
  // ctx->prm.declare_entry("min_degrees", "1", dealii::Patterns::Integer(0, 10), "Minmum degree");
  // ctx->prm.declare_entry("max_degrees", "1", dealii::Patterns::Integer(0, 10), "Maxmum degree");
  ctx->prm.declare_entry("refine_mode", "0", dealii::Patterns::Integer(0,5), "Refine_Mode");
  ctx->prm.declare_entry("refine_threshold", "0.1", dealii::Patterns::Double(0.0, 1.0), "Refine_Threshold");
}

void read_parameters(MT2DCtx *ctx, const std::string &fn) {
  ctx->prm.parse_input(fn);
  ctx->iprefix = ctx->prm.get("iprefix");
  ctx->oprefix = ctx->prm.get("oprefix");
  ctx->adaptive_rsp_oprefix = ctx->prm.get("adaptive_rsp_oprefix");
  ctx->goal_rsp_oprefix = ctx->prm.get("goal_rsp_oprefix");
  ctx->adaptive_vtk_oprefix = ctx->prm.get("adaptive_vtk_oprefix");
  ctx->goal_vtk_oprefix = ctx->prm.get("goal_vtk_oprefix");
  ctx->max_adaptive_refinements = ctx->prm.get_integer("max_adaptive_refinements");
  // ctx->refine_fraction = ctx->prm.get_double("refine_fraction");
  ctx->n_global_refinements = ctx->prm.get_integer("n_global_refinements");
  ctx->n_rx_cell_refinements = ctx->prm.get_integer("n_rx_cell_refinements");
  ctx->max_rx_cell_size = ctx->prm.get_double("max_rx_cell_size");
  // ctx->min_degree = ctx->prm.get_integer("min_degrees");
  // ctx->max_degree = ctx->prm.get_integer("max_degrees");
  ctx->refine_mode = ctx->prm.get_integer("refine_mode");
  ctx->refine_threshold = ctx->prm.get_double("refine_threshold");
}

void create_ctx(MT2DCtx *ctx) {
  ctx->fe.reset(new FE_Q(1));
  ctx->coarse_mesh.reset(new Triangulation);
  ctx->tria_finder.reset(new TriaFinder);
}

void destroy_ctx(MT2DCtx *ctx) {
  // ctx->fe = nullptr;
  ctx->coarse_mesh = nullptr;
}
