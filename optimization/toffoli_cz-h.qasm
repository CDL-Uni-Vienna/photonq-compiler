OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
cz q[1],q[2];
rz(0) q[2];
h q[2];
rz(-pi/4) q[2];
h q[2];
cz q[0],q[2];
rz(0) q[2];
h q[2];
rz(pi/4) q[2];
h q[2];
cz q[1],q[2];
rz(pi/4) q[1];
h q[1];
rz(0) q[1];
h q[1];
rz(0) q[1];
h q[1];
rz(0) q[2];
h q[2];
rz(-pi/4) q[2];
h q[2];
cz q[0],q[2];
cz q[0],q[1];
rz(0) q[0];
h q[0];
rz(pi/4) q[0];
h q[0];
rz(0) q[1];
h q[1];
rz(-pi/4) q[1];
h q[1];
cz q[0],q[1];
rz(0) q[1];
h q[1];
rz(0) q[2];
h q[2];
rz(pi/4) q[2];
h q[2];
