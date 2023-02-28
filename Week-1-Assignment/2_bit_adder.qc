Version 1.0

qubits 7

.prepare
    prep_z q[0:6]

.init
    {X q[2] | X q[4]}

display

.carry_1
    Toffoli q[1], q[2], q[3]
    CNOT q[1], q[2]
    Toffoli q[0], q[2], q[3]

.carry_2
    Toffoli q[4], q[5], q[6]
    CNOT q[4], q[5]
    Toffoli q[3], q[5], q[6]

.sum_1
    CNOT q[4], q[5]
    
    CNOT q[3], q[5]
    CNOT q[4], q[5]

.sum_2
    Toffoli q[0], q[2], q[3]
    CNOT q[1], q[2]
    Toffoli q[1], q[2], q[3]

    CNOT q[0], q[2]
    CNOT q[1], q[2]

display