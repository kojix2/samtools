#ifndef SAM_VIEW_MRUBY_H
#define SAM_VIEW_MRUBY_H

#include <htslib/sam.h>
#include <mruby.h>

struct RProc;

mrb_state *init_mruby(void);
void finalize_mruby(mrb_state *mrb);
struct RProc *compile_expression(mrb_state *mrb, const char *expr);
int evaluate_expression(mrb_state *mrb, const struct RProc *proc, const sam_hdr_t *h, bam1_t *b);

#endif
