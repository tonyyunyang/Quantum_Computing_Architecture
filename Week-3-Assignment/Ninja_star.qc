Version 1.0

qubits 17

.init
    prep_z q[0:16]
    H q[0:8]
    map q[0], D1
    map q[1], D2
    map q[2], D3
    map q[3], D4
    map q[4], D5
    map q[5], D6
    map q[6], D7
    map q[7], D8
    map q[8], D9
    map q[9], X1
    map q[10], Z2
    map q[11], X2
    map q[12], Z1
    map q[13], Z3
    map q[14], X3
    map q[15], Z4
    map q[16], X4
    display

.surface_code_cycle
    H X1
    H X2
    H X3
    H X4

    CNOT D1, Z1
    CNOT X2, D2
    CNOT D3, Z2
    CNOT D5, Z3
    CNOT X3, D6
    CNOT X4, D8

    CNOT X2, D1
    CNOT D2, Z2
    CNOT D4, Z3
    CNOT X3, D5
    CNOT D6, Z4
    CNOT X4, D7

    CNOT X1, D3
    CNOT D4, Z1
    CNOT X2, D5
    CNOT D6, Z2
    CNOT D8, Z3
    CNOT X3, D9

    CNOT X1, D2
    CNOT X2, D4
    CNOT D5, Z2
    CNOT D7, Z3
    CNOT X3, D8
    CNOT D9, Z4

    H X1
    H X2
    H X3
    H X4
    display

.syndromes_code
    Measure X1
    Measure X2
    Measure X3
    Measure X4
    Measure Z1
    Measure Z2
    Measure Z3
    Measure Z4
    display
    display_binary