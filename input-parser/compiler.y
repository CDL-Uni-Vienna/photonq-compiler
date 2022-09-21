%{

#include <stdio.h>  // For I/O
#include <stdlib.h> // For malloc here and symbol table
#include <math.h>

extern int yylex();
extern int yyparse();
extern FILE* yyin;
extern FILE* yyout;

void yyerror(const char* s);
%}

%union {
    int ival;
    float fval;
}

%token<ival> INT
%token<fval> FLOAT 
%token PI
%token RZ RX HAD CZ
%token ADD SUB MUL DIV LEFTBRACK RIGHTBRACK LEFTPARENTH RIGHTPARENTH
%token QUBIT COMMA SEMICOLON EOL

%left ADD SUB
%left MUL DIV

%type<ival> exp qubit_
%type<fval> expf arg

%%

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
;

rx_gate: RX arg qubit_ SEMICOLON { 
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                fprintf(yyout, "rz(");
                fprintf(yyout, "%f", $2);
                fprintf(yyout, ") q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                }
;

rz_gate: RZ arg qubit_ SEMICOLON { 
                fprintf(yyout, "rz(");
                fprintf(yyout, "%f", $2);
                fprintf(yyout, ") q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%i", $3);
                fprintf(yyout, "];\n");
                }
;

h_gate: HAD qubit_ SEMICOLON { 
                fprintf(yyout, "rz(0) q[");
                fprintf(yyout, "%i", $2);
                fprintf(yyout, "];\n");
                fprintf(yyout, "h q[");
                fprintf(yyout, "%i", $2);
                fprintf(yyout, "];\n");
                }
;

cz_gate: CZ qubit_ COMMA qubit_	SEMICOLON {
                fprintf(yyout, "cz q[");
                fprintf(yyout, "%i", $2);
                fprintf(yyout, "], q[");
                fprintf(yyout, "%i", $4);
                fprintf(yyout, "];\n");
                }
;

qubit_: QUBIT LEFTBRACK INT RIGHTBRACK { $$ = $3; }
;

arg: LEFTPARENTH expf RIGHTPARENTH { $$ = $2; }
;

expf: FLOAT { $$ = $1; }
| PI { $$ = M_PI; }
| expf ADD expf { $$ = $1 + $3; }
| expf SUB expf { $$ = $1 - $3; }
| expf MUL expf { $$ = $1 * $3; }
| expf DIV expf { $$ = $1 / $3; }
| SUB expf { $$ = 0 - $2 ; }
| exp ADD expf { $$ = $1 + $3; }
| exp SUB expf { $$ = $1 - $3; }
| exp MUL expf { $$ = $1 * $3; }
| exp DIV expf { $$ = $1 / $3; }
| expf ADD exp { $$ = $1 + $3; }
| expf SUB exp { $$ = $1 - $3; }
| expf MUL exp { $$ = $1 * $3; }
| expf DIV exp { $$ = $1 / $3; }
;

exp: INT { $$ = $1; }
| exp ADD exp { $$ = $1 + $3; }
| exp SUB exp { $$ = $1 - $3; }
| exp MUL exp { $$ = $1 * $3; }
| exp DIV exp { $$ = $1 / $3; }
| SUB exp { $$ = 0 - $2 ; }
;


%%

int main() {
	yyin = fopen("./test.qasm", "r");
    yyout = fopen("./output.qasm", "a");
	do {
		yyparse();
	} while(!feof(yyin));
    fclose(yyout);
	return 0;
}

void yyerror(const char* s) {
	fprintf(stderr, "Parse error: %s\n", s);
	exit(1);
}