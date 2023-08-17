#include "mt_ctx.h"
#include "mt_fwd.h"
#include "mt_io.h"

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Parameter file must be given." << std::endl;
    return 1;
  }

  MT2DCtx ctx;

  create_ctx(&ctx);

  declare_parameter(&ctx);
  read_parameters(&ctx, argv[1]);

  read_mdl(&ctx);
  read_emd(&ctx);

  mt_forward(&ctx);

  destroy_ctx(&ctx);

  return 0;
}