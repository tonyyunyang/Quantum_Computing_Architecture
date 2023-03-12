Version 0.5

qubits 9

.init
    prep_z q[0:8]
    map q[0], psi
    map q[1], o1
    map q[2], o2
    map q[3], o3
    map q[4], o4
    map q[5], o5
    map q[6], o6
    map q[7], o7
    map q[8], o8
    display

.entangle
    CNOT psi, o3
    CNOT psi, o6
    {H psi | H o3 | H o6}
    {CNOT psi, o1 | CNOT o3, o4 | CNOT o6, o7}
    {CNOT psi, o2 | CNOT o3, o5 | CNOT o6, o8}
    display