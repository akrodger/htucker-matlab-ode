function [f_out] = imp_midpoint_contraction(N,f0,f_guess,dt)
%    Implicit Midpoint contraction mapping by Bram Rodgers
%    Original Draft 01 Nov, 2021
%
%
%    Description of This Function:
%        This function computes a contraction mapping step for
%        the implicit midpoint method of order 2.
%
%	Argument List:
%        N       : A (possibly nonlinear) time independent vector field.
%        f0      : The current state of your ODE
%        f_guess : The guess for the next state of your ODE
%        dt      : The temporal step size.
%
%    Return List:
%        f_out   : A refined estimate of the next state guess.
    f_out = f0 + dt*N(((0.5*f0)+(0.5*f_guess)));
end
