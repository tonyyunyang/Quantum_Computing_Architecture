Version 1.0

qubits 2

.prepare
    prep_z q[0:1]
    map q[0], qx
    map q[1], qy

.init
    X qy

display

.uf
    CNOT qx, qy
    X qy

display