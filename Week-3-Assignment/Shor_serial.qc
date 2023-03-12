Version 0.5

qubits 17

.init
    prep_z q[0:16]
    map q[0], psi
    map q[1], o1
    map q[2], o2
    map q[3], o3
    map q[4], o4
    map q[5], o5
    map q[6], o6
    map q[7], o7
    map q[8], o8

    map q[9], a1
    map q[10], a2
    map q[11], a3
    map q[12], a4
    map q[13], a5
    map q[14], a6
    map q[15], a7
    map q[16], a8
    display

.encode
    CNOT psi, o3
    CNOT psi, o6
    {H psi | H o3 | H o6}
    {CNOT psi, o1 | CNOT o3, o4 | CNOT o6, o7}
    {CNOT psi, o2 | CNOT o3, o5 | CNOT o6, o8}
    display

.detection
    CNOT psi, a1
    CNOT o1, a1
    CNOT o1, a2
    CNOT o2, a2
    CNOT o3, a3
    CNOT o4, a3
    CNOT o4, a4
    CNOT o5, a4
    CNOT o6, a5
    CNOT o7, a5
    CNOT o7, a6
    CNOT o8, a6

    {H psi | H o1 | H o2 | H o3 | H o4 | H o5 | H o6 | H o7 | H o8}
    CNOT psi, a7
    CNOT o1, a7
    CNOT o2, a7
    CNOT o3, a7
    CNOT o4, a7
    CNOT o5, a7
    CNOT o3, a8
    CNOT o4, a8
    CNOT o5, a8
    CNOT o6, a8
    CNOT o7, a8
    CNOT o8, a8
    {H psi | H o1 | H o2 | H o3 | H o4 | H o5 | H o6 | H o7 | H o8}
    display

    {Measure a1 | Measure a2 | Measure a3 | Measure a4 | Measure a5 | Measure a6 | Measure a7 | Measure a8}
    display