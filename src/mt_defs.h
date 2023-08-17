#ifndef _MT_DEFS_H_
#define _MT_DEFS_H_

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/hp/fe_collection.h>

#include <complex>
#include <map>

typedef std::complex<double> Complex;
typedef dealii::Point<2> Point;
typedef dealii::Triangulation<2> Triangulation;
typedef dealii::DoFHandler<2> DoFHandler;
typedef dealii::FE_Q<2> FE_Q;
typedef dealii::FESystem<2> FESystem;
typedef dealii::Vector<double> VectorD;
typedef dealii::Vector<float> VectorF;
typedef dealii::Vector<Complex> VectorZ;
typedef dealii::SparseMatrix<Complex> SparseMatrixZ;
typedef dealii::AffineConstraints<Complex> AffineConstraintsZ;
typedef dealii::Quadrature<1> Quadrature;
// typedef dealii::hp::FECollection<2> FECollection;
// typedef dealii::hp::QCollection<2> QCollection2;
// typedef dealii::hp::QCollection<1> QCollection1;

typedef std::pair<double, double> Layer;
typedef std::pair<Point, Point> BoundingBox;

const double MU = 4 * dealii::numbers::PI * 1E-7;
const double EPS = 8.8541878176E-12;
const Complex II = Complex(0.0, 1.0);

const double RTOD = 180.0 / dealii::numbers::PI;
const double DTOR = dealii::numbers::PI / 180.0;

enum Mode { TE = -1, TM = -2 };

enum Type { R_XY_AP = 112, R_YX_AP = 122, T_ZX_AP = 132 };

enum Corner { Left = 0, Right = 1 };

static const std::map<int, int> TYPE_TO_MODE = { { 112, TE }, { 122, TM } };

#endif
