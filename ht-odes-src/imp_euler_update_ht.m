function [f1] = imp_euler_update_ht(N,f0,dt,err_tol)
%    HTucker Implicit Euler Update by Bram Rodgers
%    Original Draft 29 Oct, 2021
%
%
%    Description of This Function:
%        This function runs the adaptive rank fixed point iteration
%        to get the next value from implicit Euler method up
%        to a specified tolerance.
%
%	Argument List:
%        N       : A (possibly nonlinear) time independent vector field.
%                        It must be able to input and output in HT format.
%        f0      : The current state of your ODE in HTucker format.
%        dt      : The time step for your ODE solution
%        err_tol : The tolerance by which to stop contraction mapping
%                    iterations.
%    Return List:
%        f1      : Value of your ODE at the next time step. (Order 1 accurate.)
    q = 0.9;
    Q = @(f_guess) imp_euler_contraction(N,f0,f_guess,dt);
    Q_appx = @(f_guess,eps_err) ...
                imp_euler_contraction_ht(N,f0,f_guess,dt,eps_err);
    f_est = exp_midpoint_ht(N,f0,...
                                    dt,1e+6,1e+7,1e+7);
    f1 = arfpie(Q,Q_appx,q,@ht_mse,f0,err_tol);
end

