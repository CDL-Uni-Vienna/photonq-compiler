%{
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int yylex();
extern int yyparse();
extern FILE* yyin;
extern FILE* yyout;

void yyerror(const char* s);
%}

/* Declaration of data types that semantic values may have */

%union {
    const char* ival;
    const char* fval;
    const char* pi;
}

/* Token declaration with no precedence or associativity */

%token<ival> INT
%token<fval> FLOAT 
%token<pi> PI
%token RZ RX HAD CZ CX
%token ADD SUB MUL DIV LEFTBRACK RIGHTBRACK LEFTPARENTH RIGHTPARENTH
%token QUBIT COMMA SEMICOLON EOL

/* Declare token type that is left associative */

%left ADD SUB
%left MUL DIV

/* Type declaration of semantic value for non-terminal symbols with 
value types specified in %union section above */

%type<ival> qubit_
%type<fval> arg


%%
/* This section defines the PhotonQ compiler grammar rules

<non-terminal symbol>: <various terminal and non-terminal 
symbols put together> with optional {C code actions}

Distinct rules are be joined with the vertical-bar | */


input: /* empty */
	   | input line
;

line: EOL
    | gate EOL { }
;

gate: rx_gate { }
	  | rz_gate { }
	  | h_gate { }
	  | cz_gate { }
      | cx_gate { }
;

rx_gate: RX arg qubit_ SEMICOLON { 
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%s", $3);
                fprintf(yyout, "\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%s", $3);
                fprintf(yyout, "\n");
                fprintf(yyout, "rz(");
                fprintf(yyout, "%s", $2);
                fprintf(yyout, "\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%s", $3);
                fprintf(yyout, "\n");
                }
;

rz_gate: RZ arg qubit_ SEMICOLON { 
                fprintf(yyout, "rz(");
                fprintf(yyout, "%s", $2);
                fprintf(yyout, "\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%s", $3);
                fprintf(yyout, "\n");
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%s", $3);
                fprintf(yyout, "\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%s", $3);
                fprintf(yyout, "\n");
                }
;

h_gate: HAD qubit_ SEMICOLON { 
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%s", $2);
                fprintf(yyout, "\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%s", $2);
                fprintf(yyout, "\n");
                }
;

cz_gate: CZ qubit_ COMMA qubit_	SEMICOLON {
                fprintf(yyout, "cz q[");
                fprintf(yyout, "%s", $2);
                fprintf(yyout, "\n");
                }
;

cx_gate: CX qubit_ COMMA qubit_	SEMICOLON {
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%s", $4);
                fprintf(yyout, "\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%s", $4);
                fprintf(yyout, "\n");
                fprintf(yyout, "cz q[");
                fprintf(yyout, "%s", $2);
                fprintf(yyout, "\n");
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%s", $4);
                fprintf(yyout, "\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%s", $4);
                fprintf(yyout, "\n");
                }
;

qubit_: QUBIT LEFTBRACK INT RIGHTBRACK { $$ = $3; }
;

arg: LEFTPARENTH INT RIGHTPARENTH { $$ = $2; }
| LEFTPARENTH FLOAT RIGHTPARENTH { $$ = $2; }
| LEFTPARENTH PI RIGHTPARENTH { $$ = $2; }
;


%%

int main() {
	yyin = fopen("./input.qasm", "r");
    yyout = fopen("./output.qasm", "w+");
	do {
		yyparse(); /* Yacc/Bison parser */
	} while(!feof(yyin));
    fclose(yyout);
	return 0;
}

void yyerror(const char* s) {
	fprintf(stderr, "Parse error: %s\n", s);
	exit(1);
}