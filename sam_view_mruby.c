#include <htslib/sam.h>

#include "sam_view_mruby.h"

#include <mruby.h>
#include <mruby/array.h>
#include <mruby/compile.h>
#include <mruby/proc.h>
#include <mruby/string.h>
#include <mruby/variable.h>

#include <stdio.h>

static bam1_t *current_record(mrb_state *mrb)
{
    return (bam1_t *)mrb_cptr(mrb_gv_get(mrb, mrb_intern_cstr(mrb, "$b1")));
}

static const sam_hdr_t *current_header(mrb_state *mrb)
{
    return (const sam_hdr_t *)mrb_cptr(mrb_gv_get(mrb, mrb_intern_cstr(mrb, "$hdr")));
}

static mrb_value endpos_method(mrb_state *mrb, mrb_value self)
{
    return mrb_fixnum_value(bam_endpos(current_record(mrb)));
}

static mrb_value flags_method(mrb_state *mrb, mrb_value self)
{
    return mrb_fixnum_value(current_record(mrb)->core.flag);
}

static mrb_value flag_bit_method(mrb_state *mrb, uint16_t bit)
{
    return mrb_bool_value(current_record(mrb)->core.flag & bit);
}

#define DEFINE_FLAG_METHOD(name, bit)                    \
    static mrb_value name##_method(mrb_state *mrb, mrb_value self) \
    {                                                    \
        return flag_bit_method(mrb, bit);                \
    }

DEFINE_FLAG_METHOD(flag_paired, BAM_FPAIRED)
DEFINE_FLAG_METHOD(flag_proper_pair, BAM_FPROPER_PAIR)
DEFINE_FLAG_METHOD(flag_unmap, BAM_FUNMAP)
DEFINE_FLAG_METHOD(flag_munmap, BAM_FMUNMAP)
DEFINE_FLAG_METHOD(flag_reverse, BAM_FREVERSE)
DEFINE_FLAG_METHOD(flag_mreverse, BAM_FMREVERSE)
DEFINE_FLAG_METHOD(flag_read1, BAM_FREAD1)
DEFINE_FLAG_METHOD(flag_read2, BAM_FREAD2)
DEFINE_FLAG_METHOD(flag_secondary, BAM_FSECONDARY)
DEFINE_FLAG_METHOD(flag_qcfail, BAM_FQCFAIL)
DEFINE_FLAG_METHOD(flag_dup, BAM_FDUP)
DEFINE_FLAG_METHOD(flag_supplementary, BAM_FSUPPLEMENTARY)

static mrb_value hclen_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    uint32_t *cigar = bam_get_cigar(rec);
    int hclen = 0;

    for (uint32_t i = 0; i < rec->core.n_cigar; i++)
        if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)
            hclen += bam_cigar_oplen(cigar[i]);

    return mrb_fixnum_value(hclen);
}

static mrb_value mapq_method(mrb_state *mrb, mrb_value self)
{
    return mrb_fixnum_value(current_record(mrb)->core.qual);
}

static mrb_value mrefid_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    return rec->core.mtid < 0 ? mrb_nil_value() : mrb_fixnum_value(rec->core.mtid);
}

static mrb_value ncigar_method(mrb_state *mrb, mrb_value self)
{
    return mrb_fixnum_value(current_record(mrb)->core.n_cigar);
}

static mrb_value pnext_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    return rec->core.mpos < 0 ? mrb_nil_value() : mrb_fixnum_value(rec->core.mpos + 1);
}

static mrb_value pos_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    return rec->core.pos < 0 ? mrb_nil_value() : mrb_fixnum_value(rec->core.pos + 1);
}

static mrb_value qlen_method(mrb_state *mrb, mrb_value self)
{
    return mrb_fixnum_value(current_record(mrb)->core.l_qseq);
}

static mrb_value qname_method(mrb_state *mrb, mrb_value self)
{
    return mrb_str_new_cstr(mrb, bam_get_qname(current_record(mrb)));
}

static mrb_value qual_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    uint8_t *qual = bam_get_qual(rec);
    mrb_value array = mrb_ary_new_capa(mrb, rec->core.l_qseq);

    for (int32_t i = 0; i < rec->core.l_qseq; i++)
        mrb_ary_push(mrb, array, mrb_fixnum_value(qual[i]));

    return array;
}

static mrb_value refid_method(mrb_state *mrb, mrb_value self)
{
    return mrb_fixnum_value(current_record(mrb)->core.tid);
}

static mrb_value rlen_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    uint32_t *cigar = bam_get_cigar(rec);
    int rlen = 0;

    for (uint32_t i = 0; i < rec->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
            rlen += bam_cigar_oplen(cigar[i]);
    }

    return mrb_fixnum_value(rlen);
}

static mrb_value rname_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    const char *name = sam_hdr_tid2name(current_header(mrb), rec->core.tid);
    return name ? mrb_str_new_cstr(mrb, name) : mrb_nil_value();
}

static mrb_value rnext_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    if (rec->core.mtid < 0) return mrb_nil_value();

    const char *name = sam_hdr_tid2name(current_header(mrb), rec->core.mtid);
    return name ? mrb_str_new_cstr(mrb, name) : mrb_nil_value();
}

static mrb_value sclen_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    uint32_t *cigar = bam_get_cigar(rec);
    int sclen = 0;

    for (uint32_t i = 0; i < rec->core.n_cigar; i++)
        if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP)
            sclen += bam_cigar_oplen(cigar[i]);

    return mrb_fixnum_value(sclen);
}

static mrb_value seq_method(mrb_state *mrb, mrb_value self)
{
    bam1_t *rec = current_record(mrb);
    uint8_t *seq = bam_get_seq(rec);
    mrb_value seq_mrb = mrb_str_new(mrb, NULL, rec->core.l_qseq);
    char *seq_str = RSTRING_PTR(seq_mrb);

    for (int32_t i = 0; i < rec->core.l_qseq; i++)
        seq_str[i] = seq_nt16_str[bam_seqi(seq, i)];

    return seq_mrb;
}

static mrb_value tlen_method(mrb_state *mrb, mrb_value self)
{
    return mrb_fixnum_value(current_record(mrb)->core.isize);
}

static mrb_value bam_tag_method(mrb_state *mrb, mrb_value self)
{
    const char *tag;
    mrb_get_args(mrb, "z", &tag);

    uint8_t *aux = bam_aux_get(current_record(mrb), tag);
    if (!aux) return mrb_nil_value();

    switch (*aux) {
    case 'i':
    case 'I':
    case 'c':
    case 'C':
    case 's':
    case 'S':
        return mrb_fixnum_value(bam_aux2i(aux));
    case 'f':
    case 'd':
        return mrb_float_value(mrb, bam_aux2f(aux));
    case 'Z':
    case 'H':
        return mrb_str_new_cstr(mrb, bam_aux2Z(aux));
    case 'A':
        return mrb_str_new(mrb, (char *)(aux + 1), 1);
    default:
        return mrb_nil_value();
    }
}

static void register_methods(mrb_state *mrb)
{
    struct RClass *kernel = mrb->kernel_module;

    mrb_define_method(mrb, kernel, "endpos", endpos_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "flags", flags_method, MRB_ARGS_NONE());

    mrb_define_method(mrb, kernel, "paired", flag_paired_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "proper_pair", flag_proper_pair_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "unmap", flag_unmap_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "munmap", flag_munmap_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "reverse", flag_reverse_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "mreverse", flag_mreverse_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "read1", flag_read1_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "read2", flag_read2_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "secondary", flag_secondary_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "qcfail", flag_qcfail_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "dup", flag_dup_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "supplementary", flag_supplementary_method, MRB_ARGS_NONE());

    mrb_define_alias(mrb, kernel, "paired?", "paired");
    mrb_define_alias(mrb, kernel, "proper_pair?", "proper_pair");
    mrb_define_alias(mrb, kernel, "unmap?", "unmap");
    mrb_define_alias(mrb, kernel, "munmap?", "munmap");
    mrb_define_alias(mrb, kernel, "reverse?", "reverse");
    mrb_define_alias(mrb, kernel, "mreverse?", "mreverse");
    mrb_define_alias(mrb, kernel, "read1?", "read1");
    mrb_define_alias(mrb, kernel, "read2?", "read2");
    mrb_define_alias(mrb, kernel, "secondary?", "secondary");
    mrb_define_alias(mrb, kernel, "qcfail?", "qcfail");
    mrb_define_alias(mrb, kernel, "dup?", "dup");
    mrb_define_alias(mrb, kernel, "supplementary?", "supplementary");

    mrb_define_method(mrb, kernel, "hclen", hclen_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "mapq", mapq_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "mrefid", mrefid_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "ncigar", ncigar_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "pnext", pnext_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "pos", pos_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "qlen", qlen_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "qname", qname_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "qual", qual_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "refid", refid_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "rlen", rlen_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "rname", rname_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "rnext", rnext_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "sclen", sclen_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "seq", seq_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "tlen", tlen_method, MRB_ARGS_NONE());
    mrb_define_method(mrb, kernel, "tag", bam_tag_method, MRB_ARGS_REQ(1));

    mrb_define_alias(mrb, kernel, "mpos", "pnext");
    mrb_define_alias(mrb, kernel, "mrname", "rnext");
}

mrb_state *init_mruby(void)
{
    mrb_state *mrb = mrb_open();
    if (!mrb) {
        fprintf(stderr, "Could not initialize mruby\n");
        return NULL;
    }

    register_methods(mrb);
    return mrb;
}

void finalize_mruby(mrb_state *mrb)
{
    if (mrb) mrb_close(mrb);
}

struct RProc *compile_expression(mrb_state *mrb, const char *expr)
{
    if (!mrb) {
        fprintf(stderr, "mruby is not initialized\n");
        return NULL;
    }

    mrb_ccontext *cxt = mrb_ccontext_new(mrb);
    if (!cxt) {
        fprintf(stderr, "Could not initialize mruby compiler context\n");
        return NULL;
    }

    cxt->capture_errors = TRUE;
    struct mrb_parser_state *parser = mrb_parse_string(mrb, expr, cxt);
    if (!parser || parser->nerr || !parser->tree) {
        if (mrb->exc) {
            mrb_print_error(mrb);
            mrb->exc = NULL;
        }
        else if (parser && parser->nerr) {
            fprintf(stderr, "mruby syntax error at line %d: %s\n",
                    parser->error_buffer[0].lineno,
                    parser->error_buffer[0].message);
        }
        else {
            fprintf(stderr, "mruby syntax error\n");
        }

        if (parser) mrb_parser_free(parser);
        mrb_ccontext_free(mrb, cxt);
        return NULL;
    }

    struct RProc *proc = mrb_generate_code(mrb, parser);
    mrb_parser_free(parser);
    mrb_ccontext_free(mrb, cxt);

    if (!proc) {
        if (mrb->exc) {
            mrb_print_error(mrb);
            mrb->exc = NULL;
        }
        else {
            fprintf(stderr, "mruby code generation failed\n");
        }
        return NULL;
    }

    MRB_PROC_SET_TARGET_CLASS(proc, mrb->object_class);
    mrb_gv_set(mrb, mrb_intern_cstr(mrb, "$samtools_mruby_proc"), mrb_obj_value(proc));
    return proc;
}

int evaluate_expression(mrb_state *mrb, const struct RProc *proc, const sam_hdr_t *h, bam1_t *b)
{
    if (!mrb || !proc) {
        fprintf(stderr, "mruby expression is not initialized\n");
        return -1;
    }

    mrb_gv_set(mrb, mrb_intern_cstr(mrb, "$b1"), mrb_cptr_value(mrb, b));
    mrb_gv_set(mrb, mrb_intern_cstr(mrb, "$hdr"), mrb_cptr_value(mrb, (void *)h));

    int ai = mrb_gc_arena_save(mrb);
    mrb_value result = mrb_vm_run(mrb, proc, mrb_top_self(mrb), 0);
    mrb_gc_arena_restore(mrb, ai);

    if (mrb->exc) {
        mrb_print_error(mrb);
        mrb->exc = NULL;
        return -1;
    }

    return mrb_bool(result);
}
