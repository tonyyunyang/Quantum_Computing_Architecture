Version 0.5

qubits 13

.prepare
    prep_z q[0:12]
    display

.encode
    H q[0]
    H q[1]
    H q[2]
    CNOT q[6], q[5]
    CNOT q[6], q[4]
    CNOT q[0], q[3]
    CNOT q[0], q[5]
    CNOT q[0], q[6]
    CNOT q[1], q[3]
    CNOT q[1], q[4]
    CNOT q[1], q[6]
    CNOT q[2], q[3]
    CNOT q[2], q[4]
    CNOT q[2], q[5]
    display

.inject_error
    Z q[6]

.detect_error
    CNOT q[0], q[12]
    CNOT q[1], q[12]
    CNOT q[2], q[12]
    CNOT q[3], q[12]

    CNOT q[0], q[11]
    CNOT q[1], q[11]
    CNOT q[4], q[11]
    CNOT q[5], q[11]

    CNOT q[0], q[10]
    CNOT q[2], q[10]
    CNOT q[4], q[10]
    CNOT q[6], q[10]
    display

    H q[7]
    CNOT q[7], q[0]
    CNOT q[7], q[1]
    CNOT q[7], q[2]
    CNOT q[7], q[3]
    H q[7]

    H q[8]
    CNOT q[8], q[0]
    CNOT q[8], q[1]
    CNOT q[8], q[4]
    CNOT q[8], q[5]
    H q[8]
    
    H q[9]
    CNOT q[9], q[0]
    CNOT q[9], q[2]
    CNOT q[9], q[4]
    CNOT q[9], q[6]
    H q[9]
    display

    measure q[7]
    measure q[8]
    measure q[9]
    measure q[10]
    measure q[11]
    measure q[12]