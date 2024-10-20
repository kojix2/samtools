#include <mruby.h>
#include <mruby/compile.h>

#include <stdlib.h>
#include <stdio.h>

int main_tanuki(void)
{
    mrb_state *mrb;
    mrb = mrb_open();
    mrb_load_string(mrb, "puts \""
"　　　　__,-─-､__\n"
"　　　（〆-─-ヽ)\n"
"　　　 （ ´・ω・｀ ）\n"
"　　　/ 　,r‐‐‐､ヽ\n"
"　 　 し　ｌ　 x　）J\n"
"　　　_.'､ ヽ　 ノ.人\n"
"　(_((__,ノＵ´U. (酒)\n\"");
    mrb_load_string(mrb, "puts 'Tanuki in mruby'");
    mrb_close(mrb);

    return 0;
}

