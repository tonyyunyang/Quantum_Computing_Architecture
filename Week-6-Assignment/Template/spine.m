function [U, probabilities] = spine(fR, f0, dt, signal, Nplot)
% SPINE   Simulates the quantum operation of qubit(s) in response to an
%         applied signal, and optionally provides a plot (note: slower!)
%
% Arguments:
%
% fR      The Rabi frequency (in Hz) for each of the qubits.
%         If multiple values are provided, multiple qubits are simulated.
%
% f0      The Larmor frequency (in Hz) for each of the qubits.
%         If multiple values are provided, multiple qubits are simulated.
%
% dt      The fixed time step used for the signal. The signal should be
%         sufficiently oversampled for accurate simulation results.
%
% signal  The amplitude of the signal applied to the qubits (in V).
%         (amplitude of 1V -> oscillation at the specified Rabi frequency)
%
% Nplot   when  0: no plot is provided, this is fastest
%         when >0: every 'Nplot' points is plotted in a Bloch sphere as the
%                  state the qubit would have if it were initialized to the
%                  ground state. The plot is with reference to a frame
%                  rotating with f0.
%
% Return values:
%
% U       The 2x2 unitary operation performed by each of the qubits
%         U(:, :, n) gives the 2x2 unitary for qubit n
%
% probabilities
%         The measurement probability over time during the operation, along
%         the x, y and z axis.
%         E.g. probabilities(3, :, n) gives the probability along z
%         (the computational basis) over the simulated time
%

    % Determine the number of signal points
    N = length(signal);
    if (dt > 1/max(f0)/20)
        error("Timestep is too large for accurate simulation results!");
    end

    % Determine the number of qubits to simulate (single signal applied to all)
    Nq = length(fR);
    if (Nq ~= length(f0))
        error("The length of 'fR' should match the length of 'f0'!");
    end

    % Setup
    U = complex(zeros(2, 2, Nq));
    probabilities = zeros(3, 1, Nq);
    if (Nplot > 0)
        probabilities = zeros(3, floor(N/Nplot), Nq);
    end
    
    % Run in parallel if possible (Parallel Computing Toolbox required!)
    parfor nq=1:Nq
    
        % Setup
        wR = (2 * pi * fR(nq));
        w0 = (2 * pi * f0(nq));
        
        % Start from the identity operation
        dU = complex(zeros(2, 2));
        Uq = complex(eye(2));
        pq = [0; 0; 0];
        if (Nplot > 0)
            pq = zeros(3, floor(N/Nplot));
        end
        for n = 1:N

            % dt * H
            dtH = [-dt * w0 / 2.0, dt * wR * signal(n);
                    dt * wR * signal(n), dt * w0 / 2.0];

            % dU = expm(-1i * dt * H)
            a = -dtH(1,1);
            b = -dtH(1,2);
            tmp = sqrt(a*a + b*b);
            dU(1,1) = cos(tmp) + 1i * a / tmp * sin(tmp);
            dU(1,2) = 0 + 1i * b / tmp * sin(tmp);
            dU(2,1) = 0 + 1i * b / tmp * sin(tmp);
            dU(2,2) = cos(tmp) - 1i * a / tmp * sin(tmp);

            % U = dU * U
            Uq = dU * Uq;
            
            % Point to plot?
            if (Nplot > 0)
                if (mod(n, Nplot) == 0)
                    
                    % Save the coordinates of the point to plot
                    theta = dt*w0*n;
                    Px = abs(Uq(1,1)+(cos(theta)+1j*sin(theta))*Uq(2,1))^2/2;
                    Py = abs(Uq(1,1)+(sin(theta)-1j*cos(theta))*Uq(2,1))^2/2;
                    Pz = abs(Uq(1,1))^2;
                    pq(:, round(n / Nplot)) = [Px;Py;Pz];
                    
                end
            end

        end
    
        % Move to the rotating frame (unwrap the Larmor precessions)
        % dU = expm(1i * dt * H)
        a = dt * w0 * (N-1) / 2.0;
        dU(1,1) = cos(a) + 1i * sin(a);
        dU(1,2) = 0;
        dU(2,1) = 0;
        dU(2,2) = cos(a) - 1i * sin(a);
        Uq = dU * Uq;
        
        % Save the result of the individual simulation
        U(:, :, nq) = Uq;
        probabilities(:, :, nq) = pq;
        
    end

    % Plot the result in a Bloch sphere
    if (Nplot > 0)
        
        % Draw a simple Bloch sphere
        figure();
        hold on;
        axis off;
        view(-25, 25);
        camva(6);
        theta = linspace(0, pi, 50);
        for phi = 0:0.1*pi:2*pi
            plot3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), ':', 'color', 0.3*[1 1 1]);
        end
        phi = linspace(0, 2*pi, 100);
        for theta = 0:0.1*pi:pi
            plot3(sin(theta)*cos(phi), sin(theta)*sin(phi), ones(1, 100)*cos(theta), ':', 'color', 0.3*[1 1 1]);
        end
        plot3([-1.1, 1.1], [0, 0], [0, 0], '-', 'color', 0*[1 1 1]);
        plot3([0, 0], [-1.1, 1.1], [0, 0], '-', 'color', 0*[1 1 1]);
        plot3([0, 0], [0, 0], [-1.1, 1.1], '-', 'color', 0*[1 1 1]);
        xlim([-1.1, 1.1]);
        ylim([-1.1, 1.1]);
        zlim([-1.1, 1.1]);
        text(0, 0, 1, '|0>', 'fontsize', 18);
        text(0, 0, -1, '|1>', 'fontsize', 18);
    
        % Plot the trajectories
        for nq=1:Nq
            x = 2*probabilities(1, :, nq) - 1;
            y = 2*probabilities(2, :, nq) - 1;
            z = 2*probabilities(3, :, nq) - 1;
            plot3(x, y, z); hold on;
        end
        
    end

end
