function [C,Ceq] = rocket_constraint(X, g)
%ROCKET_CONSTRAINT Nonlinear inequality constraint to minimize mass.
C = [];
Ceq = g(X);
end

