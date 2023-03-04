Version 1.0

qubits 1

.prepare
    prep_z q[0]

.init
    rx q[0], 1.3
    ry q[0], 0.3
    display

.test1
    H q[0]
    display

.error
    X q[0]
    z q[0]
    display

.test2
    H q[0]
    X q[0]
    display