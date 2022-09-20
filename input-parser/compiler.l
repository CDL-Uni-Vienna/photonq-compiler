%{
#include <stdio.h>
#include <math.h> // needed for call to atof() below

#define YY_DECL int yylex()

#include "compiler.tab.h"

%}

/* Token Definitions */

integer [1-9]+[0-9]*|0
float   [0-9]+\.[0-9]+

/* Regular Expressions defining the Tokens for PhotonQ Compiler */
/* pattern + action */
%%

{integer}	{yylval.ival = atoi(yytext); return INT;}
{float} 	{yylval.fval = atof(yytext); return FLOAT;}
"pi"    { return PI; }
"rz"    { return RZ; }
"rx"    { return RX; }
"h"     { return HAD; }
"cz"    { return CZ; }
"+"     { return ADD; }
"-"     { return SUB; }
"*"     { return MUL; }
"/"     { return DIV; }
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
