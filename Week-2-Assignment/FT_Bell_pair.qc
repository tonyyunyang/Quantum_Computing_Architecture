Version 1.0

qubits 14 

.prepare
    prep_z q[0:13]
    display
    X q[13]

.encode_qB
    {H q[0] | H q[1] | H q[2]}
    {CNOT q[6], q[5] | CNOT q[6], q[4]}
    {CNOT q[0], q[3] | CNOT q[0], q[5] | CNOT q[0], q[6]}
    {CNOT q[1], q[3] | CNOT q[1], q[4] | CNOT q[1], q[6]}
    {CNOT q[2], q[3] | CNOT q[2], q[4] | CNOT q[2], q[5]}
    display

.encode_qA
    {H q[7] | H q[8] | H q[9]}
    {CNOT q[13], q[12] | CNOT q[13], q[11]}
    {CNOT q[7], q[10] | CNOT q[7], q[12] | CNOT q[7], q[13]}
    {CNOT q[8], q[10] | CNOT q[8], q[11] | CNOT q[8], q[13]}
    {CNOT q[9], q[10] | CNOT q[9], q[11] | CNOT q[9], q[12]}
    display