# photonq-compiler

## The PhotonQ Compiler converts quantum circuits expressed in OpenQASM 2.0, into an optimized graph-like representation that can be directly translated to hardware instructions on a (photonic) one-way quantum computer. The compiler is made up of three modules:
<br>

## 1. An input parser that reads [OpenQASM 2.0](https://arxiv.org/abs/1707.03429) source files and transpiles the input into a circuit consisting only of declared reference gates $\{HR_z, CZ\}$
## 2. A graph extraction and optimization module based on [ZX-calculus](https://arxiv.org/abs/1904.04735)
## 3. A hardware-dependent instruction generator that generates the photonic one-way computer's explicit hardware instructions 
<br>

# Parser

## The QASM parser is built with Unix text processing tools Flex and Bison and converts a pre-transpiled QASM circuit into a quantum circuit consisting only of $\{HR_z, CZ\}$ gates.
<br>

# Graph extraction and optimization

## This step constructs a graph representation of the circuit and applies methods from ZX-Calculus to perform optimization of the graph.
<br>

# Hardware instructions

## The last part involves execution on a discrete variable photonic platform. We currently include a simulation of the considered experimental photonic setup using [*Perceval*](https://perceval.quandela.net/docs/index.html).
