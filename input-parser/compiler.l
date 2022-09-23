%{
#include <stdio.h>

#define YY_DECL int yylex()

/* compiler.tab.h includes definitions of tokens and yylval */
#include "compiler.tab.h"
%}

/* Token Definitions */

int     [-+]?[1-9]+[0-9]*|0
float   [-+]?[0-9]+\.[0-9]+
pi      [-+]?((({int}|{float})[*\/])*pi([*]({int}|{float})|[\/]([1-9]+[0-9]*|[0-9]+\.[0-9]*[1-9]+[0-9]*))*)


/* Regular Expressions (pattern description) defining the Tokens for PhotonQ Compiler */
/* The semantic value of tokens is stored in the corresponding union member of yylval */

%%

{int}	{ yylval.ival = yytext; return INT;}
{float} { yylval.fval = yytext; return FLOAT;}
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
