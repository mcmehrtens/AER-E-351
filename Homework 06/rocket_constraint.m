function [C,Ceq] = rocket_constraint(deltaV_T, c, X, g)
%ROCKET_CONSTRAINT Nonlinear inequality constraint to minimize mass.
C = g(deltaV_T, c, X);
Ceq = [];
end

