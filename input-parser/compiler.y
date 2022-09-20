%{

#include <stdio.h>  // For I/O
#include <stdlib.h> // For malloc here and symbol table
#include <string.h> // For strcmp in symbol table

#define YYDEBUG 1   // For Debugging
int errors;         // Error Count

extern int yylex();
extern int yyparse();
extern FILE* yyin;

void yyerror(const char* s);
%}


%token INT FLOAT PI
%token RZ RX HAD CZ
%token ADD SUB MUL DIV LEFTBRACK RIGHTBRACK LEFTPARENTH RIGHTPARENTH
%token QUBIT COMMA SEMICOLON EOL

%%

input: /* empty */
	   | input line
;

line: EOL
    | gate EOL { printf("\tResult: %f\n", $1);}
;

gate: rx_gate                 		{ $$ = $1; }
	  | rz_gate                 	{ $$ = $1; }
	  | h_gate                 		{ $$ = $1; }
	  | cz_gate                     { $$ = $1; }
;

rx_gate: RX arg qubit_ SEMICOLON { printf("rz GATE"); }
;

rz_gate: RZ arg qubit_ SEMICOLON { printf("rz GATE"); }
;

h_gate: HAD qubit_ SEMICOLON { printf("HADAMARD GATE"); }
;

cz_gate: CZ qubit_ COMMA qubit_	SEMICOLON { printf("CZ GATE"); }
;

qubit_: QUBIT LEFTBRACK INT RIGHTBRACK { printf("QUBIT"); }
;

arg: LEFTPARENTH exp_ RIGHTPARENTH { printf("ARG"); }
;

exp_: exp { $$ = $1; }
| exp_ MUL exp_ { $$ = $1 * $3; }
| exp_ DIV exp_ { $$ = $1 / $3; }
| SUB exp_ { $$ = 0 - $1 ; }
;

exp: INT { $$ = $1; }
| FLOAT { $$ = $1; }
| PI { $$ = $1; }
| exp ADD exp { $$ = $1 + $3; }
| exp SUB exp { $$ = $1 - $3; }
| SUB exp { $$ = 0 - $1 ; }
;

%%

int main() {
	yyin = stdin;

	do {
		yyparse();
	} while(!feof(yyin));

	return 0;
}

void yyerror(const char* s) {
	fprintf(stderr, "Parse error: %s\n", s);
	exit(1);
}