function [c, ceq] = constraints(x)

% Givens

    [T,p,G,Z,efficiency,cost] = newParameter(x);

    T0 = 293.15;
    p0 = 1;

    % Inequalities
    c = [];
    c(1) = T0 - T(2);
    c(2) = T(2) - T(3);
    c(3) = T(4) - T(3);

    c(4) = p0 - p(2);
    c(5) = p(3) - p(2);
    c(6) = p(4) - p(3);

    c(7) = 0 - G(1);

    % Equalities
    ceq = [];
    ceq(1) = p0 - p(4);

end
