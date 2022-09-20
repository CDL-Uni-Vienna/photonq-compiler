# makefile
compiler:	compiler.l compiler.y
		bison -d compiler.y
		flex compiler.l
		cc -o $@ compiler.tab.c lex.yy.c -lfl
