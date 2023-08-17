#ifndef _MT_IO_H_
#define _MT_IO_H_

struct MT2DCtx;

void read_mdl(MT2DCtx *);
void read_emd(MT2DCtx *);

void save_rsp(MT2DCtx *, const char *);

#endif