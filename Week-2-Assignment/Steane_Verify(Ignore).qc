Version 1.0

qubits 13

.prepare
    prep_z q[0:12]
    display

.encode
    {H q[4] | H q[5] | H q[6]}
    {CNOT q[0], q[1] | CNOT q[0], q[2]}
    {CNOT q[6], q[3] | CNOT q[6], q[1] | CNOT q[6], q[0]}
    {CNOT q[5], q[3] | CNOT q[5], q[2] | CNOT q[5], q[0]}
    {CNOT q[4], q[3] | CNOT q[4], q[2] | CNOT q[4], q[1]}
    display

.inject_error

.detect_error
    {CNOT q[6], q[7] | CNOT q[4], q[7] | CNOT q[2], q[7] | CNOT q[0], q[7]}
    {CNOT q[6], q[8] | CNOT q[5], q[8] | CNOT q[2], q[8] | CNOT q[1], q[8]}
    {CNOT q[6], q[9] | CNOT q[5], q[9] | CNOT q[4], q[9] | CNOT q[3], q[9]}

    {H q[10] | H q[11] | H q[12]}
    {CNOT q[10], q[6] | CNOT q[10], q[4] | CNOT q[10], q[2] | CNOT q[10], q[0]}
    {CNOT q[11], q[6] | CNOT q[11], q[5] | CNOT q[11], q[2] | CNOT q[11], q[1]}
    {CNOT q[12], q[6] | CNOT q[12], q[5] | CNOT q[12], q[4] | CNOT q[12], q[3]}
    {H q[10] | H q[11] | H q[12]}
    
    #{Measure q[7] | Measure q[8] | Measure q[9] | Measure q[10] | Measure q[11] | Measure q[12]}

    display