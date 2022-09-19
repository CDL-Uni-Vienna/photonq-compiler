# photonq-compiler

## The PhotonQ Compiler converts OpenQASM 2.0, the de-facto standard quantum circuit language, into a graph-like representation that can be directly translated to hardware instructions on a photonic one-way quantum computer. The compiler is made up of three modules:
<br>

## 1. an input parser that reads OpenQASM 2.0 source files and transpiles the input into a circuit consisting only of declared reference gates
## 2. a graph extractor and optimizer
## 3. an instruction generator that generates the photonic one-way computer's explicit hardware instructions 
