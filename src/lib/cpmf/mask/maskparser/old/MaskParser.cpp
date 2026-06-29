/*
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301  USA
// ===============================================================================
*/

#include "MaskParser.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <ErrorSystem.hpp>

using namespace std;

//------------------------------------------------------------------------------

int                 PMFLexPosition     = 0;
struct SExpression* PMFTopExpression   = NULL;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// track info about allocation
vector<SListItem*>      PMFListItemAllocations;
vector<SList*>          PMFListAllocations;
vector<SSelection*>     PMFSelectionAllocations;
vector<SExpression*>    PMFExpressionAllocations;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool pmf_print_expression(FILE* p_fout,struct SExpression* p_expr);
bool pmf_print_selection(FILE* p_fout,struct SSelection* p_sel);
bool pmf_print_list(FILE* p_fout,struct SList* p_list);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// int main(void)
// {
//  init_mask();
//  parse_mask(":1-2@4");
//  print_expression_tree(get_expression_tree());
//  free_mask_tree();
// }

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int pmf_init_mask(void)
{
    pmf_free_mask_tree();
    PMFLexPosition = 1;
    PMFTopExpression = NULL;
    return(0);
}

//------------------------------------------------------------------------------

// int parse_mask(const char* p_mask);

//------------------------------------------------------------------------------

struct SExpression* pmf_get_expression_tree(void) {
    return(PMFTopExpression);
}

//------------------------------------------------------------------------------

int pmf_free_mask_tree(void)
{
    for(unsigned int i=0; i < PMFListItemAllocations.size(); i++) {
        delete PMFListItemAllocations[i];
    }
    PMFListItemAllocations.clear();

    for(unsigned int i=0; i < PMFListAllocations.size(); i++) {
        delete PMFListAllocations[i];
    }
    PMFListAllocations.clear();

    for(unsigned int i=0; i < PMFSelectionAllocations.size(); i++) {
        delete PMFSelectionAllocations[i];
    }
    PMFSelectionAllocations.clear();

    for(unsigned int i=0; i < PMFExpressionAllocations.size(); i++) {
        delete PMFExpressionAllocations[i];
    }
    PMFExpressionAllocations.clear();

    PMFLexPosition = 0;
    PMFTopExpression = NULL;

    return(0);
}

//------------------------------------------------------------------------------

int pmf_print_expression_tree(struct SExpression* p_expr)
{
    if(pmf_print_expression(stdout,p_expr) != true) return(-1);
    return(0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool pmf_print_expression(FILE* p_fout,struct SExpression* p_expr)
{
    if(p_expr == NULL) {
        fprintf(p_fout,"<- NULL expr");
        return(false);
    }

    fprintf(p_fout,"EXP(");
    if(p_expr->Selection != NULL) {
        if(pmf_print_selection(p_fout,p_expr->Selection) == false) return(false);
    } else {
        switch(p_expr->Operator) {
            case O_AND:
                if(pmf_print_expression(p_fout,p_expr->LeftExpression) == false) return(false);
                fprintf(p_fout," AND ");
                if(pmf_print_expression(p_fout,p_expr->RightExpression) == false) return(false);
                break;
            case O_OR:
                if(pmf_print_expression(p_fout,p_expr->LeftExpression) == false) return(false);
                fprintf(p_fout," OR ");
                if(pmf_print_expression(p_fout,p_expr->RightExpression) == false) return(false);
                break;
            case O_NOT:
                fprintf(p_fout,"NOT ");
                if(pmf_print_expression(p_fout,p_expr->RightExpression) == false) return(false);
                break;
            default:
                fprintf(p_fout,"<- unknown operator");
                return(false);
        };
    }
    fprintf(p_fout,")");
    return(true);
}

//------------------------------------------------------------------------------

bool pmf_print_selection(FILE* p_fout,struct SSelection* p_sel)
{
    if(p_sel == NULL) {
        fprintf(p_fout,"<- NULL selection");
        return(false);
    }

    fprintf(p_fout,"SEL");

    char selector;
    switch(p_sel->Type) {
        case T_RSELECTOR:
            selector = 'R';
            break;
        case T_ASELECTOR:
            selector = 'A';
            break;
        case T_TSELECTOR:
            selector = 'T';
            break;
        default:
            fprintf(p_fout,"<- unknown selector");
            return(false);
    };

    fprintf(p_fout,"[%c](",selector);
    if(pmf_print_list(p_fout,p_sel->Items) == false) return(false);
    fprintf(p_fout,")");

    return(true);
}

//------------------------------------------------------------------------------

bool pmf_print_list(FILE* p_fout,struct SList* p_list)
{
    if(p_list == NULL) {
        fprintf(p_fout,"<- NULL list");
        return(false);
    }

    struct SListItem* p_item = p_list->FirstItem;

    while(p_item != NULL) {
        if(p_item->Index < 0) {
            fprintf(p_fout,"*");
        }
        if(p_item->Index > 0) {
            if(p_item->Length > 1) {
                fprintf(p_fout,"%d-%d",p_item->Index,p_item->Index+p_item->Length);
            } else {
                if(p_item->Length == 1) {
                    fprintf(p_fout,"%d",p_item->Index);
                } else {
                    fprintf(p_fout,"<- illegal range");
                    return(false);
                }
            }
        }

        if(p_item->Index == 0) {
            fprintf(p_fout,"%-4s",p_item->Name);
        }

        p_item = p_item->NextItem;
        if(p_item != NULL) fprintf(p_fout,",");
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int pmf_yyerror(const char* p_error)
{
    ES_ERROR(p_error);
    return(0);
}

//------------------------------------------------------------------------------

int pmf_pperror(const char* p_error,int position)
{
    CSmallString error;
    error << position << ": " << p_error;
    ES_ERROR(error);
    return(0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

struct SListItem* PMFAllocateListItem(void) {
    struct SListItem* p_item = new struct SListItem;

    p_item->Index = 0;
    p_item->Length = 0;
    memset(p_item->Name,' ',4);
    p_item->NextItem = NULL;

    PMFListItemAllocations.push_back(p_item);
    return(p_item);
}

//------------------------------------------------------------------------------

struct SList* PMFAllocateList(void) {
    struct SList* p_list = new struct SList;

    p_list->FirstItem = NULL;
    p_list->LastItem = NULL;

    PMFListAllocations.push_back(p_list);
    return(p_list);
}

//------------------------------------------------------------------------------

struct SSelection* PMFAllocateSelection(enum SType type,struct SList* p_list) {
    struct SSelection* p_sel = new struct SSelection;

    p_sel->Type = type;
    p_sel->Items = p_list;

    PMFSelectionAllocations.push_back(p_sel);
    return(p_sel);
}

//------------------------------------------------------------------------------

struct SExpression* PMFAllocateExpression(void) {
    struct SExpression* p_expr = new struct SExpression;

    p_expr->Operator = O_NONE;
    p_expr->Distance = 0.0;
    p_expr->Selection = NULL;
    p_expr->LeftExpression = NULL;
    p_expr->RightExpression = NULL;

    PMFExpressionAllocations.push_back(p_expr);
    return(p_expr);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
