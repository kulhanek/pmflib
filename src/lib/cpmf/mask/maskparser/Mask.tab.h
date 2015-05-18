/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
/* Put the tokens into the symbol table, so that GDB and other debuggers
   know about them.  */
enum yytokentype {
    STRING = 258,
    INUMBER = 259,
    RNUMBER = 260,
    STAR = 261,
    COMMA = 262,
    RANGE = 263,
    RSELECTOR = 264,
    ASELECTOR = 265,
    TSELECTOR = 266,
    RLT = 267,
    RGT = 268,
    ALT = 269,
    AGT = 270,
    NOT = 271,
    AND = 272,
    OR = 273,
    RBRA = 274,
    LBRA = 275,
    ERROR = 276,
    ORIGIN = 277,
    CBOX = 278,
    LIST = 279,
    COM = 280,
    PLANE = 281
};
#endif
/* Tokens.  */
#define STRING 258
#define INUMBER 259
#define RNUMBER 260
#define STAR 261
#define COMMA 262
#define RANGE 263
#define RSELECTOR 264
#define ASELECTOR 265
#define TSELECTOR 266
#define RLT 267
#define RGT 268
#define ALT 269
#define AGT 270
#define NOT 271
#define AND 272
#define OR 273
#define RBRA 274
#define LBRA 275
#define ERROR 276
#define ORIGIN 277
#define CBOX 278
#define LIST 279
#define COM 280
#define PLANE 281




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 16 "Mask.y"
{
    SGrammar                gValue;
    SInteger                iValue;
    SReal                   rValue;
    SString                 sValue;
    struct SListItem*       itemValue;
    struct SList*           listValue;
    struct SSelection*      selValue;
    struct SExpression*     exprValue;
}
/* Line 1489 of yacc.c.  */
#line 112 "Mask.tab.h"
YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

