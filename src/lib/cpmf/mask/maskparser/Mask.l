/* =============================================================================
AMBER Mask Lexical Analyzer
      Copyright (c) 2008 Petr Kulhanek, kulhanek@chemi.muni.cz
============================================================================= */

%{
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "MaskParser.hpp"
#include "Mask.tab.h"
%}

%option noyywrap
newline       [\n]
whitesp       [ \t]+
string        [a-zA-z]+[^,&|!() \t@]*
inumber       [0-9]+
rnumber       [0-9]+[.]?[0-9]*

/* rules ==================================================================== */
%%
{inumber} {
    pmf_yylval.iValue.Position = PMFLexPosition;
    if( sscanf(yytext, "%d", &pmf_yylval.iValue.Number) != 1 ){
        pmf_pperror("unable to convert string to integer",PMFLexPosition);
        }
    PMFLexPosition += yyleng;
    return INUMBER;
    }

{rnumber} {
    pmf_yylval.rValue.Position = PMFLexPosition;
    if( sscanf(yytext, "%lf", &pmf_yylval.rValue.Number) != 1 ){
        pmf_pperror("unable to convert string to real number",PMFLexPosition);
        }
    return RNUMBER;
    }

"origin"    { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 6; return ORIGIN; }
"cbox"      { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 4; return CBOX; }
"list"      { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 4; return LIST; }
"com"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 3; return COM; }
"plane"     { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 5; return PLANE; }

{string}  {
    pmf_yylval.sValue.Position = PMFLexPosition;
    if( yyleng > 4 ){
        pmf_pperror("name is too long (max 4 characters)",PMFLexPosition);
        }
    memset(pmf_yylval.sValue.String,' ',4);
    strncpy(pmf_yylval.sValue.String,yytext,yyleng);
    PMFLexPosition += yyleng;
    return STRING;
    }

"*"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return STAR; }
","       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return COMMA; }
"-"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return RANGE; }

":"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return RSELECTOR; }
"@"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return ASELECTOR; }
"@%"      { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 2; return TSELECTOR; }

"<:"      { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 2; return RLT; }
">:"      { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 2; return RGT; }
"<@"      { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 2; return ALT; }
">@"      { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition += 2; return AGT; }

"&"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return AND; }
"|"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return OR; }
"!"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return NOT; }

"("       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return LBRA; }
")"       { pmf_yylval.gValue.Position = PMFLexPosition; PMFLexPosition++; return RBRA; }

{whitesp} { PMFLexPosition += yyleng;}

{newline} {
            pmf_pperror("new line character is not allowed in mask specification",PMFLexPosition);
            return ERROR;
            }
