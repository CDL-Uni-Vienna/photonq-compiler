%option noyywrap

%{
#include <stdio.h>
#include <string.h>
#include <math.h>

#define YY_DECL int yylex()

#include "compiler.tab.h"

%}

/* Token Definitions */

integer [1-9]+[0-9]*|0
float   [0-9]+\.[0-9]+

/* Regular Expressions defining the Tokens for PhotonQ Compiler */
/* pattern + action */
%%

{float} 	{yylval = atof(yytext); return FLOAT;}
{integer}	{yylval = atoi(yytext); return INT;}
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
// int yywrap(void){}