#include"mt_ctx.h"
#include"mt_defs.h"
#include"reEstimate.h"

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include<vector>
#include<memory>

class Get_Rho : public dealii::Function<2>
{
public:
    Get_Rho(MT2DCtx *ctx, FEMData *fd);
    
    void value_list(const std::vector<Point> &point_list, std::vector<double> &value_list,
                        unsigned int component) const;

    MT2DCtx *ctx;
    FEMData *fd;
    double tolerance = 1.e-10;
};


void write_output_error(FEMData *fd, int fidx, VectorF &error );

void write_output(FEMData *fd, int fidx, int index);

void refine_global(MT2DCtx *ctx, FEMData *fd, int fidx);


void refine_adaptive(MT2DCtx *ctx, FEMData *fd,int fidx, int index);


void refine_goal(MT2DCtx *ctx, FEMData *fd, int fidx, int index);
