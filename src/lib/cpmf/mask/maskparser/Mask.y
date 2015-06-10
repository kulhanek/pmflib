/* =============================================================================
AMBER Mask Parser Analyzer
      Copyright (c) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
============================================================================= */

%{
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include "MaskParser.hpp"
#define YYERROR_VERBOSE
%}

%union{
    SGrammar                gValue;
    SInteger                iValue;
    SReal                   rValue;
    SString                 sValue;
    struct SListItem*       itemValue;
    struct SList*           listValue;
    struct SSelection*      selValue;
    struct SExpression*     exprValue;
    };

/* recognized tokens -------------------------------------------------------- */
%token <sValue> STRING
%token <iValue> INUMBER
%token <rValue> RNUMBER
%token <gValue> STAR COMMA RANGE
%token <gValue> RSELECTOR ASELECTOR TSELECTOR
%token <gValue> RLT RGT ALT AGT
%token <gValue> NOT AND OR
%token <gValue> RBRA LBRA
%token <gValue> ERROR
%token <gValue> ORIGIN CBOX LIST COM PLANE

%type <itemValue> item
%type <listValue> list
%type <listValue> sel_list
%type <selValue>  selection
%type <exprValue> expr
%type <exprValue> amber_mask

/* priority of operators ---------------------------------------------------- */
%left OR
%left AND
%left NOT
%left RLT RGT ALT AGT
%left ORIGIN CBOX LIST COM PLANE

/* defined grammar ---------------------------------------------------------- */
%%

amber_mask:
    expr                { PMFTopExpression = $1; $$ = $1; }
    ;

expr:
/* SELECTION ---------------------------------------------------------------- */

    selection {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Selection = $1;
        $$ = p_expr;
        }

/* LOGICAL OPERATORS -------------------------------------------------------- */

    | expr AND expr {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AND; 
        p_expr->LeftExpression = $1;
        p_expr->RightExpression = $3;
        $$ = p_expr;
        }
    | expr OR expr {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_OR;
        p_expr->LeftExpression = $1;
        p_expr->RightExpression = $3;
        $$ = p_expr;
        }
    | NOT expr {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_NOT;
        p_expr->LeftExpression = NULL;
        p_expr->RightExpression = $2;
        $$ = p_expr;
        }

/* BRACKETS ----------------------------------------------------------------- */

    | LBRA expr RBRA {
        $$ = $2;
        }

/* DISTANCE OPERATORS ------------------------------------------------------- */
    | ORIGIN ALT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_ALT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_ORIGIN;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | CBOX ALT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_ALT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_CBOX;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | expr ALT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_ALT;
        p_expr->LeftExpression = $1;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | LIST LBRA expr RBRA ALT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_ALT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | COM LBRA expr RBRA ALT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_ALT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_COM;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | PLANE LBRA expr RBRA ALT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_ALT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_PLANE;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }

    | ORIGIN AGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AGT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_ORIGIN;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | CBOX AGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AGT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_CBOX;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | expr AGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AGT;
        p_expr->LeftExpression = $1;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | LIST LBRA expr RBRA AGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AGT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | COM LBRA expr RBRA AGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AGT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_COM;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | PLANE LBRA expr RBRA AGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AGT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_PLANE;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }

    | ORIGIN RLT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RLT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_ORIGIN;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | CBOX RLT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RLT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_CBOX;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | expr RLT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RLT;
        p_expr->LeftExpression = $1;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | LIST LBRA expr RBRA RLT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RLT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | COM LBRA expr RBRA RLT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RLT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_COM;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | PLANE LBRA expr RBRA RLT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RLT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_PLANE;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }

    | ORIGIN RGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RGT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_ORIGIN;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | CBOX RGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RGT;
        p_expr->LeftExpression = NULL;
        p_expr->Modificator = D_CBOX;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | expr RGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RGT;
        p_expr->LeftExpression = $1;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $3.Number;
        $$ = p_expr;
        }
    | LIST LBRA expr RBRA RGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RGT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_LIST;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | COM LBRA expr RBRA RGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RGT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_COM;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }
    | PLANE LBRA expr RBRA RGT RNUMBER {
        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_RGT;
        p_expr->LeftExpression = $3;
        p_expr->Modificator = D_PLANE;
        p_expr->Distance = $6.Number;
        $$ = p_expr;
        }

/* SELECTORS ---------------------------------------------------------------- */

    | RSELECTOR sel_list ASELECTOR sel_list {
        struct SSelection* p_lsel;
        p_lsel = PMFAllocateSelection(T_RSELECTOR,$2);
        if( p_lsel == NULL ){
            pmf_yyerror("unable to allocate memory for the residue selection");
            YYERROR;
            }
        struct SExpression* p_lexpr;
        p_lexpr = PMFAllocateExpression();
        if( p_lexpr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_lexpr->Selection = p_lsel;

        struct SSelection* p_rsel;
        p_rsel = PMFAllocateSelection(T_ASELECTOR,$4);
        if( p_rsel == NULL ){
            pmf_yyerror("unable to allocate memory for the atom selection");
            YYERROR;
            }
        struct SExpression* p_rexpr;
        p_rexpr = PMFAllocateExpression();
        if( p_rexpr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_rexpr->Selection = p_rsel;

        struct SExpression* p_expr;
        p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AND;
        p_expr->LeftExpression = p_lexpr;
        p_expr->RightExpression = p_rexpr;
        $$ = p_expr;
        }

    | RSELECTOR sel_list TSELECTOR sel_list {
        struct SSelection* p_lsel = PMFAllocateSelection(T_RSELECTOR,$2);
        if( p_lsel == NULL ){
            pmf_yyerror("unable to allocate memory for the residue selection");
            YYERROR;
            }
        struct SExpression* p_lexpr = PMFAllocateExpression();
        if( p_lexpr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_lexpr->Selection = p_lsel;

        struct SSelection* p_rsel = PMFAllocateSelection(T_TSELECTOR,$4);
        if( p_rsel == NULL ){
            pmf_yyerror("unable to allocate memory for the atom selection");
            YYERROR;
            }
        struct SExpression* p_rexpr = PMFAllocateExpression();
        if( p_rexpr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_rexpr->Selection = p_rsel;

        struct SExpression* p_expr = PMFAllocateExpression();
        if( p_expr == NULL ){
            pmf_yyerror("unable to allocate memory for the expression");
            YYERROR;
            }
        p_expr->Operator = O_AND;
        p_expr->LeftExpression = p_lexpr;
        p_expr->RightExpression = p_rexpr;
        $$ = p_expr;
        }
    ;

selection:
    RSELECTOR sel_list {
        struct SSelection* p_selection = PMFAllocateSelection(T_RSELECTOR,$2);
        if( p_selection == NULL ){
            pmf_yyerror("unable to allocate memory for the residue selection");
            YYERROR;
            }
        $$ = p_selection;
        }

    | ASELECTOR sel_list {
        struct SSelection* p_selection = PMFAllocateSelection(T_ASELECTOR,$2);
        if( p_selection == NULL ){
            pmf_yyerror("unable to allocate memory for the atom selection");
            YYERROR;
            }
        $$ = p_selection;
        }

    | TSELECTOR sel_list{
        struct SSelection* p_selection = PMFAllocateSelection(T_TSELECTOR,$2);
        if( p_selection == NULL ){
            pmf_yyerror("unable to allocate memory for the type selection");
            YYERROR;
            }
        $$ = p_selection;
        }
    ;

sel_list:
    list {
        $$ = $1;
        }
    | STAR {
        struct SListItem* p_item = PMFAllocateListItem();
        if( p_item == NULL ){
            pmf_yyerror("unable to allocate memory for the item");
            YYERROR;
            }
        p_item->Index = -1;

        struct SList* p_list = PMFAllocateList();
        if( p_list == NULL ){
            pmf_yyerror("unable to allocate memory for the list");
            YYERROR;
            }
        /* add everything item to the list */
        p_list->FirstItem = p_item;
        p_list->LastItem = p_item;

        $$ = p_list;
        }
    ;

list:
    item {
        /* create list node */
        struct SList* p_list = PMFAllocateList();
        if( p_list == NULL ){
            pmf_yyerror("unable to allocate memory for the list");
            YYERROR;
            }
        /* add item to the list */
        p_list->FirstItem = $1;
        p_list->LastItem = $1;
        $$ = p_list;
        }
    | list COMMA item {
        /* add item to the list */
        $1->LastItem->NextItem = $3;
        $1->LastItem = $3;
        $$ = $1;
        }
    ;

item:
    INUMBER {
        struct SListItem* p_item = PMFAllocateListItem();
        if( p_item == NULL ){
            pmf_yyerror("unable to allocate memory for the item");
            YYERROR;
            }
        p_item->Index = $1.Number;
        p_item->Length = 1;
        $$ = p_item;
        }
    | INUMBER RANGE INUMBER {
        struct SListItem* p_item = PMFAllocateListItem();
        if( p_item == NULL ){
            pmf_yyerror("unable to allocate memory for the item");
            YYERROR;
            }
        p_item->Index = $1.Number;
        p_item->Length = $3.Number - $1.Number + 1;
        $$ = p_item;
        }
    | STRING {
        struct SListItem* p_item = PMFAllocateListItem();
        if( p_item == NULL ){
            pmf_yyerror("unable to allocate memory for the item");
            YYERROR;
            }
        strncpy(p_item->Name,$1.String,4);
        $$ = p_item;
        }
    ;

%%
/* ========================================================================== */

int pmf_parse_mask(const char* p_mask)
{
 pmf_yy_scan_string(p_mask);
 int rvalue = pmf_yyparse();
 pmf_yylex_destroy();
 return(rvalue);
}

