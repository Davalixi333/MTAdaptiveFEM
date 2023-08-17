SET(mt2d_SOURCES mt_ctx.cc mt_io.cc mt_mesh.cc reEstimate.cc mt_refine.cc mt_fwd.cc tethex.cc mt2d.cc)

ADD_EXECUTABLE(mt2d ${mt2d_SOURCES})
DEAL_II_SETUP_TARGET(mt2d)
