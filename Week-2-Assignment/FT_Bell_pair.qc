Version 1.0

qubits 14 

.prepare
    prep_z q[0:13]
    map q[0], q0
    map q[1], q1
    map q[2], q2
    map q[3], q3
    map q[4], q4
    map q[5], q5
    map q[6], q6
    map q[7], q7
    map q[8], q8
    map q[9], q9
    map q[10], q10
    map q[11], q11
    map q[12], q12
    map q[13], q13

.encode_qA
    H q0
    H q1
    H q2
    CNOT q6, q5
    CNOT q6, q4
    CNOT q0, q3
    CNOT q0, q5
    CNOT q0, q6
    CNOT q1, q3
    CNOT q1, q4
    CNOT q1, q6
    CNOT q2, q3
    CNOT q2, q4
    CNOT q2, q5
    X q0
    X q1
    X q2
    X q3
    X q4
    X q5
    X q6

.encode_qB
    H q7
    H q8
    H q9
    CNOT q13, q12
    CNOT q13, q11
    CNOT q7, q10
    CNOT q7, q12
    CNOT q7, q13
    CNOT q8, q10
    CNOT q8, q11
    CNOT q8, q13
    CNOT q9, q10
    CNOT q9, q11
    CNOT q9, q12

.logical_H
    H q0
    H q1
    H q2
    H q3
    H q4
    H q5
    H q6

.logical_CNOT
    CNOT q0, q7
    CNOT q1 , q8
    CNOT q2 , q9
    CNOT q3 , q10
    CNOT q4 , q11
    CNOT q5 , q12
    CNOT q6 , q13

.decoding
    X q6
    X q5
    X q4
    X q3
    X q2
    X q1
    X q0
    CNOT q2, q5
    CNOT q2, q4
    CNOT q2, q3
    CNOT q1, q6
    CNOT q1, q4
    CNOT q1, q3
    CNOT q0, q6
    CNOT q0, q5
    CNOT q0, q3
    CNOT q6, q4
    CNOT q6, q5
    H q2
    H q1
    H q0

    CNOT q9, q12
    CNOT q9, q11
    CNOT q9, q10
    CNOT q8, q13
    CNOT q8, q11
    CNOT q8, q10
    CNOT q7, q13
    CNOT q7, q12
    CNOT q7, q10
    CNOT q13, q11
    CNOT q13, q12
    H q9
    H q8
    H q7
    display

.measurement
    Measure q6
    Measure q13