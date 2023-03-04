Version 1.0

qubits 9

.prepare
    prep_z q[0:8]
    display

.init
    rx q[0], 1.3
    ry q[0], 0.3
    display

.phaseflip_encode
    CNOT q[0], q[3]
    CNOT q[0], q[6]
    {H q[0] | H q[3] | H q[6]}
    display

.bitflip_encode
    {CNOT q[0], q[1] | CNOT q[3], q[4] | CNOT q[6], q[7]}
    {CNOT q[0], q[2] | CNOT q[3], q[5] | CNOT q[6], q[8]}
    display

.inject_error
    X q[0]
    display

.bitflip_decode_and_correct
    {CNOT q[0], q[1] | CNOT q[3], q[4] | CNOT q[6], q[7]}
    {CNOT q[0], q[2] | CNOT q[3], q[5] | CNOT q[6], q[8]}
    {Toffoli q[2], q[1], q[0] | Toffoli q[5], q[4], q[3] | Toffoli q[8], q[7], q[6]}
    display

.phaseflip_decode_and_correct
    {H q[0] | H q[3] | H q[6]}
    CNOT q[0], q[3]
    CNOT q[0], q[6]
    Toffoli q[6], q[3], q[0]
    display