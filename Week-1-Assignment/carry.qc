Version 1.0

qubits 4

.prepare
    prep_z q[0:3]

.init
    {X q[1] | X q[2]}

.carry
    Toffoli q[1], q[2], q[3]
    CNOT q[1], q[2]
    Toffoli q[0], q[2], q[3]

.measurement
    measure q[3]

display