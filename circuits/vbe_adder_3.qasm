OPENQASM 2.0;
include "qelib1.inc";

gate ccz a,b,c
{
  h c;
  h c;
  cx b,c; tdg c;
  cx a,c; t c;
  cx b,c; tdg c;
  cx a,c; t b; t c; h c;
  cx a,b; t a; tdg b;
  cx a,b;
  h c;
}

qreg q[10];
h q[3];
h q[6];
h q[9];
ccz q[1], q[2], q[3];
cx q[1], q[2];
ccz q[0], q[2], q[3];
h q[3];
ccz q[4], q[5], q[6];
cx q[4], q[5];
ccz q[3], q[5], q[6];
h q[6];
ccz q[7], q[8], q[9];
cx q[7], q[8];
ccz q[6], q[8], q[9];
cx q[6], q[8];
h q[9];
h q[6];
ccz q[3], q[5], q[6];
cx q[4], q[5];
ccz q[4], q[5], q[6];
cx q[3], q[5];
cx q[4], q[5];
h q[6];
h q[3];
ccz q[0], q[2], q[3];
cx q[1], q[2];
ccz q[1], q[2], q[3];
cx q[0], q[2];
cx q[1], q[2];
h q[3];