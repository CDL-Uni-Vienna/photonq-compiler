%{
#include <stdio.h>
#include <math.h>

#define YY_DECL int yylex()

#include "compiler.tab.h"

%}

/* Token Definitions */

integer [1-9]+[0-9]*|0
float   [0-9]+\.[0-9]+
pi      [-+]?pi|[-+]?[0-9].[0-9]+[*]pi

/* Regular Expressions defining the Tokens for PhotonQ Compiler */
/* pattern + action */
%%

{integer}	{yylval.ival = yytext; return INT;}
{float} 	{yylval.fval = yytext; return FLOAT;}
{pi}    { yylval.pi = yytext; return PI; }
"rz"    { return RZ; }
"rx"    { return RX; }
"h"     { return HAD; }
"cz"    { return CZ; }
"cx"    { return CX; }
"["     { return LEFTBRACK; }
"]"     { return RIGHTBRACK; }
"("     { return LEFTPARENTH; }
")"     { return RIGHTPARENTH; }
"q"     { return QUBIT; }
","     { return COMMA; }
";"     { return SEMICOLON; }
\n      { return EOL; }
[ \t]   { /* ignore whitespace */ }
.       { printf("Unknown character \n"); }

%%
