/****************************************************************\
*                                                                *
*  Module of splice site and intron modelling                    *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2009.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU General Public License, version 3. See the file COPYING   *
*  or http://www.gnu.org/licenses/gpl.txt for details            *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include "intron.h"
#include "ungapped.h"

Intron_ArgumentSet *Intron_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Intron_ArgumentSet ias;
    if(arg){
        as = ArgumentSet_create("Intron Modelling Options");
        ArgumentSet_add_option(as, 0, "minintron", "length",
            "Minimum intron length", "30",
            Argument_parse_int, &ias.min_intron);
        ArgumentSet_add_option(as, 0, "maxintron", "length",
            "Maximum intron length", "200000",
            Argument_parse_int, &ias.max_intron);
        ArgumentSet_add_option(as, 'i', "intronpenalty", "score",
            "Intron Opening penalty", "-30",
            Argument_parse_int, &ias.intron_open_penalty);
        Argument_absorb_ArgumentSet(arg, as);
        ias.sps = NULL;
        if(ias.intron_open_penalty > 0)
            g_warning("Intron open penalty should be negative [%d]",
                      ias.intron_open_penalty);
    } else {
        if(!ias.sps)
            ias.sps = SplicePredictorSet_create();
            /* FIXME: this should be freed somewhere */
        }
    return &ias;
    }

/**/

static gchar *Intron_ChainData_get_seq_func(gint start,
                                            gint len,
                                            gpointer user_data){
    register Sequence *s = user_data;
    return Sequence_get_substr(s, start, len);
    }

void Intron_ChainData_init_splice_prediction(Intron_ChainData *icd,
        Sequence *s, SplicePredictorSet *sps, SpliceType type,
        gboolean use_single){
    if(!icd->sps)
        icd->sps = g_new0(SplicePrediction_Set, 1);
    icd->sps = SplicePrediction_Set_add(icd->sps, sps, type, s->len,
            Intron_ChainData_get_seq_func, use_single, s);
    return;
    }

static void Intron_ChainData_init(Intron_ChainData *icd){
    icd->curr_intron_start = 0;
    return;
    }

static Intron_ChainData *Intron_ChainData_create(void){
    register Intron_ChainData *icd = g_new0(Intron_ChainData, 1);
    Intron_ChainData_init(icd);
    return icd;
    }

static void Intron_ChainData_destroy(Intron_ChainData *icd){
    if(icd->sps)
        SplicePrediction_Set_destroy(icd->sps);
    g_free(icd);
    return;
    }

/**/

Intron_Data *Intron_Data_create(void){
    register Intron_Data *id = g_new0(Intron_Data, 1);
    id->ias = Intron_ArgumentSet_create(NULL);
    id->query_data = Intron_ChainData_create();
    id->target_data = Intron_ChainData_create();
    return id;
    }

void Intron_Data_destroy(Intron_Data *id){
    Intron_ChainData_destroy(id->query_data);
    Intron_ChainData_destroy(id->target_data);
    g_free(id);
    return;
    }

/**/

static void intron_init_func(Region *region, gpointer user_data){
    register Ungapped_Data *ud = user_data;
    register Intron_Data *id = Ungapped_Data_get_Intron_Data(ud);
    if(!id->curr_region){
        id->curr_region = Region_share(region);
        Intron_ChainData_init(id->query_data);
        Intron_ChainData_init(id->target_data);
        }
    return;
    }

static gchar *intron_init_macro =
    "if(!id->curr_region){\n"
    "    id->curr_region = Region_share(region);\n"
    "    Intron_ChainData_init(id->query_data);\n"
    "    Intron_ChainData_init(id->target_data);\n"
    "    }\n";

static void intron_exit_func(Region *region, gpointer user_data){
    register Ungapped_Data *ud = user_data;
    register Intron_Data *id = Ungapped_Data_get_Intron_Data(ud);
    if(id->curr_region){
        Region_destroy(id->curr_region);
        id->curr_region = NULL;
        }
    return;
    }

static gchar *intron_exit_macro =
    "if(id->curr_region){\n"
    "    Region_destroy(id->curr_region);\n"
    "    id->curr_region = NULL;\n"
    "    }\n";

/**/

#define Intron_CalcFunc(site, chain, is_pre)                       \
    static C4_Score Intron_calc_##site##_##chain(                  \
        gint query_pos, gint target_pos, gpointer user_data){      \
        register Ungapped_Data *ud = user_data;                    \
        register Intron_Data *id                                   \
                = Ungapped_Data_get_Intron_Data(ud);               \
        register C4_Score score = 0;                               \
        register gint intron_length;                               \
        g_assert(id->chain##_data->sps->site);                     \
        g_assert(id->curr_region);                                 \
        if(is_pre){                                                \
            score = id->ias->intron_open_penalty;                  \
        } else {                                                   \
            intron_length = chain##_pos                            \
                          - id->chain##_data->curr_intron_start    \
                          + 2;                                     \
            if((intron_length < id->ias->min_intron)               \
            || (intron_length > id->ias->max_intron))              \
                return C4_IMPOSSIBLY_LOW_SCORE;                    \
            }                                                      \
        score += SplicePrediction_get(id->chain##_data->sps->site, \
                                      chain##_pos);                \
        return score;                                              \
        }

#define Intron_JointCalcFunc(site, is_pre)                          \
    static C4_Score Intron_joint_calc_##site(                       \
        gint query_pos, gint target_pos, gpointer user_data){       \
        register Ungapped_Data *ud = user_data;                     \
        register Intron_Data *id                                    \
                = Ungapped_Data_get_Intron_Data(ud);                \
        register C4_Score score = 0;                                \
        register gint intron_length;                                \
        g_assert(id->query##_data->sps->site);                      \
        g_assert(id->target##_data->sps->site);                     \
        g_assert(id->curr_region);                                  \
        if(is_pre){                                                 \
            score = id->ias->intron_open_penalty;                   \
        } else {                                                    \
            intron_length = query_pos                               \
                          - id->query_data->curr_intron_start       \
                          + 2;                                      \
            if((intron_length < id->ias->min_intron)                \
            || (intron_length > id->ias->max_intron))               \
                return C4_IMPOSSIBLY_LOW_SCORE;                     \
            intron_length = target_pos                              \
                          - id->target_data->curr_intron_start      \
                          + 2;                                      \
            if((intron_length < id->ias->min_intron)                \
            || (intron_length > id->ias->max_intron))               \
                return C4_IMPOSSIBLY_LOW_SCORE;                     \
            }                                                       \
        score += SplicePrediction_get(id->query##_data->sps->site,  \
                                      query_pos);                   \
        score += SplicePrediction_get(id->target##_data->sps->site, \
                                      target_pos);                  \
        return score;                                               \
        }

static gchar *Intron_get_calc_macro_splice(gboolean on_query,
                                           gchar *site){
    register gchar *chain = on_query?"query":"target";
    register gchar *chain_macro = on_query?"QP":"TP";
    return g_strdup_printf("SplicePrediction_get(id->%s_data->sps->%s, %%%s)",
            chain, site, chain_macro);
    }

static gchar *Intron_get_calc_macro_length(gboolean on_query,
                                           gchar *input_macro){
    register gchar *chain = on_query?"query":"target";
    register gchar *chain_macro = on_query?"QP":"TP";
    register gchar *intron_length = g_strdup_printf(
                    "(%%%s - id->%s_data->curr_intron_start + 2)",
                     chain_macro, chain);
    register gchar *macro = g_strdup_printf(
                    "(((%s < id->ias->min_intron)"
                    "||(%s > id->ias->max_intron))"
                    "?C4_IMPOSSIBLY_LOW_SCORE:%s)",
                    intron_length, intron_length, input_macro);
    g_free(intron_length);
    return macro;
    }

static gchar *Intron_get_calc_macro(gboolean is_5_prime,
                                    gboolean is_forward,
                                    gboolean on_query,
                                    gboolean on_target,
                                    gboolean is_pre){
    register gchar *macro = NULL;
    register gchar *site = g_strdup_printf("ss%d_%s", is_5_prime?5:3,
                                      is_forward?"forward":"reverse");
    register gchar *splice_calc;
    register gchar *a, *b;
    if(on_query && on_target){ /* joint */
        a = Intron_get_calc_macro_splice(TRUE, site);
        b = Intron_get_calc_macro_splice(FALSE, site);
        splice_calc = g_strdup_printf("(%s + %s)", a, b);
        g_free(a);
        g_free(b);
    } else { /* single */
        splice_calc = Intron_get_calc_macro_splice(on_query, site);
        }
    /**/
    if(is_pre){
        macro = g_strdup_printf("(%s + %s)",
                     "id->ias->intron_open_penalty", splice_calc);
    } else {
        if(on_query && on_target){ /* joint */
            a = Intron_get_calc_macro_length(TRUE, splice_calc);
            macro = Intron_get_calc_macro_length(FALSE, a);
            g_free(a);
        } else { /* single */
            macro = Intron_get_calc_macro_length(on_query, splice_calc);
            }
        }
    g_free(splice_calc);
    g_free(site);
    g_assert(macro);
    return macro;
    }

#define Intron_InitFunc(site, chain)                                         \
    static void Intron_init_##site##_##chain(                                \
                             Region *region, gpointer user_data){            \
        register Ungapped_Data *ud = user_data;                              \
        register Intron_Data *id = Ungapped_Data_get_Intron_Data(ud);        \
        Intron_ChainData_init_splice_prediction(id->chain##_data, ud->chain, \
                                           id->ias->sps,                     \
                                           SpliceType_##site,                \
                                           FALSE);                           \
    return;                                                                  \
    }

#define Intron_JointInitFunc(site)                                 \
    static void Intron_joint_init_##site(                          \
                             Region *region, gpointer user_data){  \
    Intron_init_##site##_query(region, user_data);                 \
    Intron_init_##site##_target(region, user_data);                \
    return;                                                        \
    }

static gchar *Intron_get_strand_init_macro(gboolean is_5_prime,
                                           gboolean is_forward,
                                           gboolean on_query){
    register gchar *chain = on_query?"query":"target";
    register gchar *site = g_strdup_printf("ss%d_%s",
            is_5_prime?5:3, is_forward?"forward":"reverse");
    register gchar *macro = g_strdup_printf(
        "Intron_ChainData_init_splice_prediction(id->%s_data, ud->%s, id->ias->sps,"
        "                                        SpliceType_%s,"
        "        ((id->curr_region->%s_start == 0)"
        "        && (id->curr_region->%s_length == ud->%s->len)));",
        chain, chain, site, chain, chain, chain);
    g_free(site);
    return macro;
    }

static gchar *Intron_get_init_macro(gboolean is_5_prime,
                                    gboolean is_forward,
                                    gboolean on_query,
                                    gboolean on_target){
    register gchar *qy_macro, *tg_macro, *macro;
    if(on_query && on_target){
        qy_macro = Intron_get_strand_init_macro(is_5_prime, is_forward,
                                                TRUE);
        tg_macro = Intron_get_strand_init_macro(is_5_prime, is_forward,
                                                FALSE);
        macro = g_strdup_printf("%s;%s", qy_macro, tg_macro);
        g_free(qy_macro);
        g_free(tg_macro);
    } else {
        macro = Intron_get_strand_init_macro(is_5_prime, is_forward,
                                             on_query);
        }
    return macro;
    }

#if 0
#define Intron_InitMacro(site, chain, chainmacro)               \
static gchar *intron_## site ##_## chain ##_init_macro =        \
    "if(!id->"#chain"_data->"#site"){\n"                        \
    "    id->"#chain"_data->"#site" = SplicePrediction_create(" \
    "          id->ias->sps->"# site ",\n"                      \
    "          ud->"#chain"->seq,\n"                            \
    "          ud->"#chain"->len,\n"                            \
    "          %"#chainmacro"S, %"#chainmacro"L);\n"            \
    "    }\n";
/* FIXME: do not give len when doing DP over unknown region
 *        or give hint, eg: computing_full_dp
 */
#endif /* 0 */

#if 0
#define Intron_ExitFunc(site, chain)                                  \
    static void Intron_exit_##site##_##chain(                         \
                       Region *region, gpointer user_data){           \
        register Ungapped_Data *ud = user_data;                       \
        register Intron_Data *id = Ungapped_Data_get_Intron_Data(ud); \
        if(id->chain##_data->site){                                   \
            SplicePrediction_destroy(id->chain##_data->site);         \
            id->chain##_data->site = NULL;                            \
            }                                                         \
        return;                                                       \
        }

#define Intron_JointExitFunc(site)                                    \
    static void Intron_joint_exit_##site(                             \
                       Region *region, gpointer user_data){           \
        Intron_exit_##site##_query(region, user_data);                \
        Intron_exit_##site##_target(region, user_data);               \
        return;                                                       \
        }
#endif /* 0 */

#if 0
#define Intron_ExitMacro(site, chain)                                 \
static gchar *intron_##site##_##chain##_exit_macro =                  \
        "if(id->"#chain"_data->"#site"){\n"                           \
        "    SplicePrediction_destroy(id->"#chain"_data->"#site");\n" \
        "    id->"#chain"_data->"#site" = NULL;"                      \
        "    }\n";
#endif /* 0 */

#if 0
static gchar *Intron_get_strand_exit_macro(gboolean is_5_prime,
                                           gboolean is_forward,
                                           gboolean on_query){
    register gchar *chain = on_query?"query":"target";
    register gchar *site = g_strdup_printf("ss%d_%s",
            is_5_prime?5:3, is_forward?"forward":"reverse");
    return g_strdup_printf(
        "if(id->%s_data->%s){\n"
        "    SplicePrediction_destroy(id->%s_data->%s);\n"
        "    id->%s_data->%s = NULL;"
        "    }\n",
        chain, site, chain, site, chain, site);
    }

static gchar *Intron_get_exit_macro(gboolean is_5_prime,
                                    gboolean is_forward,
                                    gboolean on_query,
                                    gboolean on_target){
    register gchar *qy_macro, *tg_macro, *macro;
    if(on_query && on_target){
        qy_macro = Intron_get_strand_exit_macro(is_5_prime, is_forward,
                                                TRUE);
        tg_macro = Intron_get_strand_exit_macro(is_5_prime, is_forward,
                                                FALSE);
        macro = g_strdup_printf("%s;%s", qy_macro, tg_macro);
        g_free(qy_macro);
        g_free(tg_macro);
    } else {
        macro = Intron_get_strand_exit_macro(is_5_prime, is_forward,
                                             on_query);
        }
    return macro;
    }
#endif /* 0 */

/**/

#define Intron_CalcElements(site, chain, chainmacro, is_pre) \
    Intron_CalcFunc(site, chain, is_pre)                     \
    Intron_InitFunc(site, chain)

    /* Intron_ExitFunc(site, chain) */

#define Intron_JointCalcElements(site, is_pre) \
    Intron_JointCalcFunc(site, is_pre)         \
    Intron_JointInitFunc(site)


    /* Intron_JointExitFunc(site) */

/* pre */
Intron_CalcElements(ss5_forward, query, Q, TRUE)
Intron_CalcElements(ss5_forward, target, T, TRUE)
Intron_JointCalcElements(ss5_forward, TRUE)

/* post */
Intron_CalcElements(ss3_forward, query, Q, FALSE)
Intron_CalcElements(ss3_forward, target, T, FALSE)
Intron_JointCalcElements(ss3_forward, FALSE)

/* pre */
Intron_CalcElements(ss3_reverse, query, Q, TRUE)
Intron_CalcElements(ss3_reverse, target, T, TRUE)
Intron_JointCalcElements(ss3_reverse, FALSE)

/* post */
Intron_CalcElements(ss5_reverse, query, Q, FALSE)
Intron_CalcElements(ss5_reverse, target, T, FALSE)
Intron_JointCalcElements(ss5_reverse, FALSE)

#define Intron_assign_calc_elements(site, chain)       \
    calc_func = Intron_calc_##site##_##chain;          \
    init_func = Intron_init_##site##_##chain;

    /* exit_func = Intron_exit_##site##_##chain; */

/**/

#define Intron_assign_joint_calc_elements(site)    \
    calc_func = Intron_joint_calc_##site;          \
    init_func = Intron_joint_init_##site;

    /* exit_func = Intron_joint_exit_##site; */

/**/

/*
Functions and macros for intron length shadow:

Intron_{query,target}_{start,end}_{func,macro}
*/

#define Intron_StartFunc(chain)                           \
static C4_Score Intron_##chain##_start_func(              \
    gint query_pos, gint target_pos, gpointer user_data){ \
    return (C4_Score)chain##_pos;                         \
    }

#define Intron_StartMacro(chain, chainmacro) \
static gchar *Intron_##chain##_start_macro   \
    = "(C4_Score)(%"#chainmacro"P)";

#define Intron_StartElements(chain, chainmacro) \
    Intron_StartFunc(chain)                     \
    Intron_StartMacro(chain, chainmacro)

#define Intron_EndFunc(chain)                                     \
static void Intron_##chain##_end_func(C4_Score score,             \
    gint query_pos, gint target_pos, gpointer user_data){         \
    register Ungapped_Data *ud = user_data;                       \
    register Intron_Data *id = Ungapped_Data_get_Intron_Data(ud); \
    g_assert(id);                                                 \
    id->chain##_data->curr_intron_start = (gint)score;            \
    return;                                                       \
    }

#define Intron_EndMacro(chain)                            \
static gchar *Intron_##chain##_end_macro                  \
    = "id->"#chain"_data->curr_intron_start = (gint)%SS";

/**/

#define Intron_EndElements(chain) \
    Intron_EndFunc(chain)         \
    Intron_EndMacro(chain)

#define Intron_ShadowElements(chain, chainmacro) \
    Intron_StartElements(chain, chainmacro)      \
    Intron_EndElements(chain)

Intron_ShadowElements(query, Q)
Intron_ShadowElements(target, T)

/**/

static C4_Calc *Intron_add_calc(C4_Model *model,
                                gchar *prefix, gchar *suffix,
                                gboolean use_pre,
                                gboolean is_forward,
                                gboolean on_query, gboolean on_target){
    register gchar *name = g_strconcat(prefix, " ", suffix, NULL);
    register Intron_ArgumentSet *ias
           = Intron_ArgumentSet_create(NULL);
    register SplicePredictor *sp;
    register C4_Calc *calc;
    register C4_CalcFunc calc_func;
    register C4_PrepFunc init_func;
    register gchar *calc_macro, *init_macro;
    register C4_Score bound = 0;
    register gboolean is_5_prime = is_forward?(use_pre?TRUE:FALSE)
                                             :(use_pre?FALSE:TRUE);
    if(is_forward){
        if(use_pre){
            sp = ias->sps->ss5_forward;
            bound = ias->intron_open_penalty;
            if(on_query){
                if(on_target){
                    Intron_assign_joint_calc_elements(ss5_forward);
                } else {
                    Intron_assign_calc_elements(ss5_forward, query);
                    }
            } else {
                Intron_assign_calc_elements(ss5_forward, target);
                }
        } else {
            sp = ias->sps->ss3_forward;
            if(on_query){
                if(on_target){
                    Intron_assign_joint_calc_elements(ss3_forward);
                } else {
                    Intron_assign_calc_elements(ss3_forward, query);
                    }
            } else {
                Intron_assign_calc_elements(ss3_forward, target);
                }
            }
    } else {
        if(use_pre){
            sp = ias->sps->ss3_reverse;
            bound = ias->intron_open_penalty;
            if(on_query){
                if(on_target){
                    Intron_assign_joint_calc_elements(ss3_reverse);
                } else {
                    Intron_assign_calc_elements(ss3_reverse, query);
                    }
            } else {
                Intron_assign_calc_elements(ss3_reverse, target);
                }
        } else {
            sp = ias->sps->ss5_reverse;
            if(on_query){
                if(on_target){
                    Intron_assign_joint_calc_elements(ss5_reverse);
                } else {
                    Intron_assign_calc_elements(ss5_reverse, query);
                    }
            } else {
                Intron_assign_calc_elements(ss5_reverse, target);
                }
            }
        }
    bound += SplicePredictor_get_max_score(sp);
    if(on_query && on_target) /* Double for joint intron */
        bound += SplicePredictor_get_max_score(sp);
    calc_macro = Intron_get_calc_macro(is_5_prime, is_forward,
                                       on_query, on_target, use_pre);
    init_macro = Intron_get_init_macro(is_5_prime, is_forward,
                                       on_query, on_target);
    /*
    exit_macro = Intron_get_exit_macro(is_5_prime, is_forward,
                                       on_query, on_target);
    */
    calc = C4_Model_add_calc(model, name, bound,
                             calc_func, calc_macro,
                             init_func, init_macro,
                             NULL, NULL,
                             C4_Protect_UNDERFLOW);
    g_assert(calc_macro);
    g_assert(init_macro);
    g_free(calc_macro);
    g_free(init_macro);
    g_free(name);
    return calc;
    }

C4_Model *Intron_create(gchar *suffix, gboolean on_query,
                        gboolean on_target, gboolean is_forward){
    register gchar *name, *pre_name, *post_name;
    register C4_Model *model;
    register C4_State *intron_state;
    register gint qy_splice_advance, tg_splice_advance;
    register C4_Calc *pre_calc, *post_calc;
    register C4_Label pre_label, post_label;
    register Intron_ArgumentSet *ias
           = Intron_ArgumentSet_create(NULL);
    g_assert(on_query || on_target);
    name = g_strconcat("intron ", suffix, NULL);
    model = C4_Model_create(name);
    g_free(name);
    if(is_forward){
        pre_name = "5'ss forward";
        post_name = "3'ss forward";
        pre_label = C4_Label_5SS;
        post_label = C4_Label_3SS;
    } else {
        pre_name = "3'ss reverse";
        post_name = "5'ss reverse";
        pre_label = C4_Label_3SS;
        post_label = C4_Label_5SS;
        }
    if(on_query){
        if(on_target){ /* Joint intron */
            qy_splice_advance = 2;
            tg_splice_advance = 2;
        } else {
            qy_splice_advance = 2;
            tg_splice_advance = 0;
            }
    } else { /* on target */
        qy_splice_advance = 0;
        tg_splice_advance = 2;
        }
    /* Add calcs */
    pre_calc = Intron_add_calc(model, pre_name, suffix, TRUE,
                    is_forward, on_query, on_target);
    post_calc = Intron_add_calc(model, post_name, suffix, FALSE,
                    is_forward, on_query, on_target);
    /* Add intron state */
    name = g_strconcat("intron ", suffix, NULL);
    intron_state = C4_Model_add_state(model, name);
    g_free(name);
    /* Add transition START->intron */
    name = g_strconcat("(START) to ", intron_state->name, NULL);
    C4_Model_add_transition(model, name, NULL, intron_state,
       qy_splice_advance, tg_splice_advance, pre_calc, pre_label, NULL);
    g_free(name);
    /* Add transition intron->intron */
    if(on_query){
        name = g_strconcat("query intron loop ", suffix, NULL);
        C4_Model_add_transition(model, name, intron_state, intron_state,
            1, 0, NULL, C4_Label_INTRON, NULL);
        g_free(name);
        }
    if(on_target){
        name = g_strconcat("target intron loop ", suffix, NULL);
        C4_Model_add_transition(model, name, intron_state, intron_state,
            0, 1, NULL, C4_Label_INTRON, NULL);
        g_free(name);
        }
    /* Add transition intron->(END)*/
    name = g_strconcat(intron_state->name, " to (END)", NULL);
    C4_Model_add_transition(model, name,
        intron_state, NULL, qy_splice_advance, tg_splice_advance,
        post_calc, post_label, NULL);
    g_free(name);
    /* Add span state */
    name = g_strconcat("intron span", suffix, NULL);
    if(on_query){
        if(on_target){
            C4_Model_add_span(model, name, intron_state,
                              ias->min_intron, ias->max_intron,
                              ias->min_intron, ias->max_intron);
        } else {
            C4_Model_add_span(model, name, intron_state,
                              ias->min_intron, ias->max_intron, 0, 0);
            }
    } else {
        C4_Model_add_span(model, name, intron_state,
                          0, 0, ias->min_intron, ias->max_intron);
        }
    g_free(name);
    /* Add shadows */
    if(on_query){
        name = g_strconcat("query intron ", suffix, NULL);
        C4_Model_add_shadow(model, name, NULL, NULL,
                Intron_query_start_func, Intron_query_start_macro,
                Intron_query_end_func, Intron_query_end_macro);
        g_free(name);
        }
    if(on_target){
        name = g_strconcat("target intron ", suffix, NULL);
        C4_Model_add_shadow(model, name, NULL, NULL,
                Intron_target_start_func, Intron_target_start_macro,
                Intron_target_end_func, Intron_target_end_macro);
        g_free(name);
        }
    /* Add DP init/exit funcs/macros */
    C4_Model_configure_extra(model,
        intron_init_func, intron_init_macro,
        intron_exit_func, intron_exit_macro);
    C4_Model_append_codegen(model, NULL,
            "register Intron_Data *id = ud->intron_data;\n", NULL);
    C4_Model_close(model);
    return model;
    }

