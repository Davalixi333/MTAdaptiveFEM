#include"mt_refine.h"
#include"reEstimate.h"
#include <deal.II/grid/grid_refinement.h>


#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/mapping.h>

#include <fstream>
Get_Rho::Get_Rho(MT2DCtx *ctx, FEMData *fd)
{
    this->ctx = ctx;
    this->fd = fd;
}                  

void Get_Rho::value_list(const std::vector<Point> &point_list, std::vector<double> &value_list,
                        unsigned int component) const
{
    auto cell_hint = Triangulation::active_cell_iterator();
    AssertDimension(value_list.size(), point_list.size());
    for (unsigned int p = 0; p< point_list.size(); ++p)
    {
        std::vector<bool> marked_vertices = {};
        // Get_Rho::vector_value(point_list[p], value_list[p]);
        cell_hint = dealii::GridTools::find_active_cell_around_point(
                        this->fd->mesh, point_list[p],marked_vertices,this->tolerance);
        if(cell_hint != this->fd->mesh.end())
        {
            // cell_hint = cell_and_point.first;
            value_list[p]= 1./this->ctx->rho[cell_hint->user_index()];
        }
    }
    (void)component;

}

void write_output_error(FEMData *fd, int fidx, VectorF &error )
{
    dealii::DataOut<2> data_out;
    data_out.attach_dof_handler(fd->dh);
	data_out.add_data_vector(fd->u, "U");
	data_out.add_data_vector(error, "error");
	data_out.build_patches();

	std::ofstream output("grid-adp" + std::to_string(fidx)+ ".vtk");
	data_out.write_vtk(output);
}

void write_output(FEMData *fd, int fidx, int index)
{
    dealii::DataOut<2> data_out;
    data_out.attach_dof_handler(fd->dh);
	data_out.add_data_vector(fd->u, "U");
	
	data_out.build_patches();

	std::ofstream output("grid-block" + std::to_string(fidx) + std::to_string(index)+ ".vtk");
	data_out.write_vtk(output);
}

void refine_global(MT2DCtx *ctx, FEMData *fd, int fidx)
{
    VectorF estimated_error_per_cell(fd->mesh.n_active_cells());
    Re_Estimate(ctx, fd, fd->u, estimated_error_per_cell);

    write_output_error(fd, fidx, estimated_error_per_cell);

    dealii::GridRefinement::refine_and_coarsen_fixed_number(fd->mesh,
		                                    estimated_error_per_cell,
		                                    ctx->refine_threshold,
		                                    0);
	fd->mesh.execute_coarsening_and_refinement();
    for (auto cell = fd->mesh.begin(0); cell != fd->mesh.end(0); ++cell) {
		cell->recursively_set_user_index(cell->user_index());
	}
}

void refine_adaptive(MT2DCtx *ctx, FEMData *fd,int fidx, int index)
{
    // dealii::ComponentMask mask;
    // Get_Rho gh(ctx, fd);
    // dealii::Function<2> *fun = &gh;
    VectorF estimated_error_per_cell(fd->mesh.n_active_cells());
	dealii::KellyErrorEstimator<2>::estimate(fd->dh,
                                      dealii::QGauss<1>(ctx->fe->degree + 1),
                                      {},
                                      fd->u,
                                      estimated_error_per_cell);
                                    //   mask,
                                    //   fun);


	dealii::DataOut<2> data_out;
	data_out.attach_dof_handler(fd->dh);
	data_out.add_data_vector(fd->u, "U");
	data_out.add_data_vector(estimated_error_per_cell, "error");
	data_out.build_patches();

	std::ofstream output(ctx->adaptive_vtk_oprefix + "-adap-" + std::to_string(index-1) + std::to_string(fidx) + ".vtk");
	data_out.write_vtk(output);

	dealii::GridRefinement::refine_and_coarsen_fixed_number(fd->mesh,
		                                    estimated_error_per_cell,
                                            0.3,
		                                    0);
	fd->mesh.execute_coarsening_and_refinement();

	for (auto cell = fd->mesh.begin(0); cell != fd->mesh.end(0); ++cell) {
		cell->recursively_set_user_index(cell->user_index());
	}
}

void refine_goal(MT2DCtx *ctx, FEMData *fd, int fidx, int index)
{
    VectorF estimated_error_per_cell(fd->mesh.n_active_cells());
    dealii::KellyErrorEstimator<2>::estimate(fd->dh,
                                     dealii::QGauss<1>(ctx->fe->degree + 1),
                                      {},
                                     fd->u,
                                     estimated_error_per_cell);
  
    VectorF delta_error(fd->mesh.n_active_cells());
    dealii::KellyErrorEstimator<2>::estimate(fd->dh,
                                     dealii::QGauss<1>(ctx->fe->degree + 1),
                                      {},
                                     fd->delta_s,
                                     delta_error);

    VectorF total_estimated_error(fd->mesh.n_active_cells());
    for(unsigned int i = 0; i < total_estimated_error.size(); i++)
    {
        total_estimated_error[i] = estimated_error_per_cell[i] + delta_error[i];
    }
    dealii::DataOut<2> data_out;
    data_out.attach_dof_handler(fd->dh);
    data_out.add_data_vector(fd->u, "U");
    data_out.add_data_vector(estimated_error_per_cell, "error");
    data_out.add_data_vector(total_estimated_error, "total_error");
    data_out.build_patches();
    
    std::ofstream output(ctx->goal_vtk_oprefix + "-adap-" + std::to_string(index-1) + std::to_string(fidx) +".vtk");
    data_out.write_vtk(output);

    dealii::GridRefinement::refine_and_coarsen_fixed_number(fd->mesh,
                                                    total_estimated_error,
                                                    0.3,
                                                    0);
    fd->mesh.execute_coarsening_and_refinement();
    
    //refine_receiving_area(ctx, fd, fidx);
    for (auto cell = fd->mesh.begin(0); cell != fd->mesh.end(0); ++cell) {
        cell->recursively_set_user_index(cell->user_index());
    }
}