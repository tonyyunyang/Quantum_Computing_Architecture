Version 1.0

qubits 3

.prepare
    prep_z q[0:2]

.init
    Rx q[0], pi/2

.epr
    H q[1]
    CNOT q[1], q[2]

.encode
    CNOT q[0], q[1]
    H q[0]
    measure q[0]
    measure q[1]

.decode
    c-x b[1], q[2]
    c-z b[0], q[2]

.measurement
    measure q[2]

display