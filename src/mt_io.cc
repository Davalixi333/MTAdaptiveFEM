#include "mt_ctx.h"
#include "tethex.h"

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>

static std::string parse_string(const std::string &s) {
  size_t pos;
  std::string l = s;

  pos = l.find("#");
  if (pos != std::string::npos) {
    l.erase(pos, std::string::npos);
  }

  pos = l.find_first_not_of(" \t");
  if (pos != 0 && pos != std::string::npos) {
    l.erase(l.begin(), l.begin() + pos);
  }

  pos = l.find_last_not_of(" \t");
  if (pos != (l.size() - 1) && pos != std::string::npos) {
    l.erase(pos + 1, std::string::npos);
  }

  if (l.size() > 0) {
    l += "\n";
  }

  return l;
}

void read_mdl(MT2DCtx *ctx) {
  std::string l;
  std::stringstream ss;
  std::vector<Point> points;
  std::vector<dealii::CellData<2> > cells;
  size_t n_points, n_cells, i, j, type, idx;

  tethex::Mesh mesh;
  std::vector<size_t> corners;
  std::vector<double> vertices;

  std::ifstream ifs(ctx->iprefix + ".mdl");
  if (!ifs.good()) {
  }

  while (std::getline(ifs, l)) {
    ss << parse_string(l);
  }

  ss >> n_points;
  vertices.resize(n_points * 2);
  for (i = 0; i < n_points; ++i) {
    ss >> vertices[i * 2 + 0] >> vertices[i * 2 + 1];
  }

  ss >> n_cells >> type;
  corners.resize(n_cells * (type + 1));
  for (i = 0; i < n_cells; ++i) {
    for (j = 0; j < type; ++j) {
      ss >> corners[i * (type + 1) + j];
    }
    ss >> idx;
    corners[(type + 1) * i + j] = idx;
  }

  ss >> n_cells;
  ctx->rho.reinit(n_cells);
  ctx->rho_lb.reinit(n_cells);
  ctx->rho_ub.reinit(n_cells);

  for (i = 0; i < n_cells; ++i) {
    ss >> ctx->rho(i) >> ctx->rho_lb(i) >> ctx->rho_ub(i);
  }

  mesh.read(type, vertices, corners);
  mesh.convert();

  points.resize(mesh.get_n_vertices());
  for (i = 0; i < points.size(); ++i) {
    points[i][0] = mesh.get_vertex(i).get_coord(0);
    points[i][1] = mesh.get_vertex(i).get_coord(1);
  }

  cells.resize(mesh.get_n_quadrangles());
  for (i = 0; i < cells.size(); ++i) {
    for (j = 0; j < 4; ++j) {
      cells[i].vertices[j] = mesh.get_quadrangle(i).get_vertex(j);
    }
    if (type == 3) {
      std::swap(cells[i].vertices[2], cells[i].vertices[3]);
    }
    cells[i].manifold_id = mesh.get_quadrangle(i).get_material_id();
  }

  ctx->coarse_mesh->clear();
  ctx->coarse_mesh->create_triangulation(points, cells, dealii::SubCellData());

  for (auto cell = ctx->coarse_mesh->begin_active(); cell != ctx->coarse_mesh->end(); ++cell) {
    cell->set_user_index(cell->manifold_id());
    cell->set_manifold_id(dealii::numbers::flat_manifold_id);
  }

  ctx->tria_finder->reinit(*ctx->coarse_mesh);
}

void read_emd(MT2DCtx *ctx) {
  std::string l;
  std::stringstream ss;
  double re, im, err_re, err_im;
  size_t n_freqs, n_rxes, n_obses, i;

  std::ifstream ifs(ctx->iprefix + ".emd");
  if (!ifs.good()) {
  }

  while (std::getline(ifs, l)) {
    ss << parse_string(l);
  }

  ss >> n_freqs;
  ctx->freqs.resize(n_freqs);
  for (i = 0; i < n_freqs; ++i) {
    ss >> ctx->freqs[i];
  }

  ss >> n_rxes;
  ctx->rx.resize(n_rxes);
  for (i = 0; i < n_rxes; ++i) {
    ss >> ctx->rx[i][0] >> ctx->rx[i][1];
  }

  ss >> n_obses;
  ctx->otype.resize(n_obses);
  ctx->fidx.resize(n_obses);
  ctx->ridx.resize(n_obses);
  ctx->obs.resize(n_obses);
  ctx->rsp.resize(n_obses);
  ctx->obserr.resize(n_obses);

  std::fill(ctx->rsp.begin(), ctx->rsp.end(), Complex(0.0));

  for (i = 0; i < n_obses; ++i) {
    ss >> ctx->otype[i] >> ctx->fidx[i] >> ctx->ridx[i];
    ss >> re >> im >> err_re >> err_im;
    switch (ctx->otype[i]) {
    case R_XY_AP:
    case R_YX_AP:
      err_re = err_re / 100.0;
      err_im = err_im * DTOR;
      ctx->obs[i] = Complex(std::log(re), im * 2 * DTOR);
      ctx->obserr[i] = Complex(err_re * err_re, (err_im * 2) * (err_im * 2));
      break;
    default:
      ctx->obs[i] = Complex(re, im);
      ctx->obserr[i] = Complex(err_re * err_re, err_im * err_im);
      break;
    }
  }
}

void save_rsp(MT2DCtx *ctx, const char *fn) {
  FILE *fp;
  Complex obs, fm;
  int i, nrx, nfreq, nobs;

  fp = fopen(fn, "w");

  nfreq = ctx->freqs.size();
  fprintf(fp, "# frequencies\n");
  fprintf(fp, "%d\n", nfreq);
  for (i = 0; i < nfreq; ++i) {
    fprintf(fp, "%.4E\n", ctx->freqs[i]);
  }

  nrx = ctx->rx.size();
  fprintf(fp, "# recievers\n");
  fprintf(fp, "%d\n", nrx);
  for (i = 0; i < nrx; ++i) {
    fprintf(fp, "% .4E % .4E\n", ctx->rx[i][0], ctx->rx[i][1]);
  }

  nobs = ctx->obs.size();
  fprintf(fp, "# observations\n");
  fprintf(fp, "%d\n", nobs);
  fprintf(fp, "#%6s %6s %7s %13s %13s %13s %13s\n", "otype", "fidx", "ridx", "amp", "phase", "amp`", "phase`");
  for (i = 0; i < nobs; ++i) {
    switch (ctx->otype[i]) {
    case R_XY_AP:
    case R_YX_AP:
      obs = Complex(std::exp(std::real(ctx->obs[i])), std::imag(ctx->obs[i]) * RTOD / 2);
      fm = Complex(std::exp(std::real(ctx->rsp[i])), std::imag(ctx->rsp[i]) * RTOD / 2);
      fprintf(fp, "% 7d % 6d % 7d % 13.6E % 13.3f % 13.6E % 13.3f\n", ctx->otype[i],
              ctx->fidx[i], ctx->ridx[i], std::real(obs), std::imag(obs), std::real(fm), std::imag(fm));
      break;
    default:
      obs = ctx->obs[i];
      fm = ctx->rsp[i];
      fprintf(fp, "% 7d % 6d % 7d % 13.6E % 13.6E % 13.6E % 13.6E\n", ctx->otype[i],
              ctx->fidx[i], ctx->ridx[i], std::real(obs), std::imag(obs), std::real(fm), std::imag(fm));
      break;
    }
  }

  fclose(fp);
}
