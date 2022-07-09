function [f1] = imp_midpoint_update_ht(N,f0,dt,err_tol)
%    Tensor Train Implicit Midpoint Update by Bram Rodgers
%    Original Draft 29 Oct, 2021
%
%
%    Description of This Function:
%        This function runs the adaptive rank fixed point iteration
%        to get the next value from implicit midpoint method up
%        to a specified tolerance.
%
%	
%	Argument List:
%        N       : A (possibly nonlinear) time independent vector field.
%                        It must be able to input and output in TT format.
%        f0      : The current state of your ODE in Tensor Train format.
%        dt      : The time step for your ODE solution
%        err_tol : The tolerance by which to stop contraction mapping
%                    iterations.
%    Return List:
%        f1      : Value of your ODE at the next time step. (Order 2 accurate.)
    q = 0.5;
    Q = @(f_guess) imp_midpoint_contraction(N,f0,f_guess,dt);
    Q_appx = @(f_guess,eps_err) ...
                imp_midpoint_contraction_ht(N,f0,f_guess,dt,eps_err);
    %f1 = exp_midpoint_ht(N,f0,dt,1e+5,1e+6,1e+6);
    f1 = arfpie(Q,Q_appx,q,@tt_mse,f0,err_tol);
end

