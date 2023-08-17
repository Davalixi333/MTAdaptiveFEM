#include"mt_ctx.h"
#include"mt_defs.h"

#include <deal.II/base/numbers.h>
#include <deal.II/fe/fe_values.h>



void Re_Estimate(MT2DCtx *ctx, FEMData *fd, VectorZ& solution, VectorF& error);

void integrate_on_face();


void integrate_on_regular_cell();

void integrate_on_irregular_cell();
