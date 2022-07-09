function [f_out] = imp_euler_contraction(N,f0,f_guess,dt)
%
%    Implicit Euler contraction mapping by Bram Rodgers
%    Original Draft 29 Oct, 2021
%
%
%    Description of This Function:
%        This function returns the contraction mapping output
%        for Implicit Euler method.
%
%	Argument List:
%        N       : A (possibly nonlinear) time independent vector field.
%        f0      : The current state of your ODE
%        f_guess : The guess for the next state of your ODE
%        dt      : The temporal step size.
%
%    Return List:
%        f_out   : A refined estimate of the next state guess.
%
    f_out = f0 + dt*N(f_guess);     %Formula for Euler with contraction update.
end
