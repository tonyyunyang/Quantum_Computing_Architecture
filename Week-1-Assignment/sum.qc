Version 1.0

qubits 3

.prepare
    prep_z q[0:2]

.init
    X q[1]

.sum
    CNOT q[0], q[2]
    CNOT q[1], q[2]

.measurement
    measure q[2]
    display