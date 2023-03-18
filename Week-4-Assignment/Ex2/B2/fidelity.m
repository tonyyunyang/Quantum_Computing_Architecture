function F = fidelity(U, Uideal)
% FIDELITY  Calculates the fidelity of the operation U against the ideal
%           operation given in Uideal
%
% Arguments:
%
% U         The 2x2 unitary describing the real, non-ideal operation
%
% Uideal    The 2x2 unitary describing the ideal operation
%
% Return values:
%
% F         The fidelity of the operation
%

    % Return the process fidelity
    F = abs(trace(Uideal' * U))^2/4;

end
