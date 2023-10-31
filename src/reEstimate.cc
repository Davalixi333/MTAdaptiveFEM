#include"reEstimate.h"

#include<deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include<vector>
#include<map>
#include<list>
#include<ctime>
#include<algorithm>
#include<fstream>

void Re_Estimate(MT2DCtx *ctx, FEMData *fd, VectorZ& solution, VectorF& error)
{
    dealii::QGauss<1> face_quadrature_formula(ctx->fe->degree + 1);

    dealii::FEFaceValues<2> fe_face_values(*(ctx->fe), face_quadrature_formula,
                                dealii::update_values | dealii::update_JxW_values | dealii::update_gradients |
                                dealii::update_normal_vectors);
    dealii::FEFaceValues<2> fe_face_neighbor_values(*(ctx->fe), face_quadrature_formula,
                                dealii::update_values | dealii::update_JxW_values | dealii::update_gradients |
                                dealii::update_normal_vectors);                      
    
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    // unsigned int n_component = fd->dh.get_fe(0).n_components();
    double sigme;
    double sigme_neighbor;
   
    std::vector<dealii::Tensor<1, 2, Complex>> fe_face_gradient(n_face_q_points);
    std::vector<dealii::Tensor<1, 2, Complex>> fe_face_neighbor_gradient(n_face_q_points);
    unsigned int i = 0;
    for(const auto &cell : fd->dh.active_cell_iterators())
    {
     
        float cell_error = 0;
        sigme = 1./ctx->rho[cell->user_index()];
        
        for(const unsigned int face_num : cell->face_indices())
        {   
            float jump = 0;
            fe_face_values.reinit(cell, face_num);

            if(cell->face(face_num)->at_boundary() == false)
            {
                DoFHandler::active_cell_iterator neighbor =
                                        cell->neighbor(face_num);
                sigme_neighbor = 1. / ctx->rho[neighbor->user_index()];

                unsigned int neighbor_neighbor = 
                                        cell->neighbor_of_neighbor(face_num);

                fe_face_neighbor_values.reinit(neighbor, neighbor_neighbor);
                
                fe_face_values.get_function_gradients(solution, fe_face_gradient);
                fe_face_neighbor_values.get_function_gradients(solution, fe_face_neighbor_gradient);

                for(unsigned int q = 0; q < n_face_q_points; ++q)
                {
                    jump += dealii::numbers::NumberTraits<Complex>::abs_square(sigme * fe_face_gradient[q] * fe_face_values.normal_vector(q) +
                            sigme_neighbor * fe_face_neighbor_gradient[q] * fe_face_neighbor_values.normal_vector(q)) *
                             fe_face_values.JxW(q);
                    
                }
                cell_error += jump;
            }
        }
        error[i] = std::sqrt(cell_error);
        i++;
    }
}