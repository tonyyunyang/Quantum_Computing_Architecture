Version 0.5

qubits 9

.init
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
    {prep_z psi | prep_z o3 | prep_z o6}
    {CNOT psi, o3 | CNOT psi, o6}
    {H psi | H o3 | H o6 | prep_z o1 | prep_z o4 | prep_z o7 | prep_z o2 | prep_z o5 | prep_z o8}
    {CNOT psi, o1 | CNOT o3, o4 | CNOT o6, o7 | CNOT psi, o2 | CNOT o3, o5 | CNOT o6, o8}
    display