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
    char* str;
}

%token<ival> INT
%token<str> PI
%token RZ RX HAD CZ
%token ADD SUB MUL DIV LEFTBRACK RIGHTBRACK LEFTPARENTH RIGHTPARENTH
%token QUBIT COMMA SEMICOLON EOL

%left ADD SUB
%left MUL DIV

%type<ival> qubit_
%type<str> pi_ arg

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
                fprintf(yyout, "%s", $2);
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
                fprintf(yyout, "%s", $2);
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

arg: LEFTPARENTH pi_ RIGHTPARENTH { $$ = (char *)$2; }
;

pi_: PI { }
;





%%

int main() {
	yyin = fopen("./input.qasm", "r");
    yyout = fopen("./output.qasm", "w+");
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