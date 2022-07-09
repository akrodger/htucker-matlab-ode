function [f_out] =  imp_midpoint_contraction_tt(N,f0,f_guess,dt,err_tol)
%
%    HT format Implicit Midpoint contraction by Bram Rodgers
%    Original Draft 01 Nov, 2021
%
%
%    Description of This Function:
%        This function computes an approximate contraction mapping step for
%        the implicit midpoint method of order 2. It is guaranteed to be at
%        most eps_er in error from the output of imp_midpoint_contraction.m 
%
%	Argument List:
%        N       : A (possibly nonlinear) time independent vector field.
%                        It must be able to input and output in HT format.
%        f0      : The current state of your ODE in HTucker format.
%        f_guess : The guess for the next state of your ODE
%                        in HTucker format.
%        dt      : The temporal step size.
%        eps_err : An error to achieve so that the output is at most eps_err
%                        away from implicit_euler.m
%
%    Return List:
%        f_out   : A refined estimate of the next state guess in HTucker
%                       format.
%
    a = 1.0/3.0;
    b = a;
    c = 1.0-a-b;
    e_1 = a*err_tol;
    e_2 = b*err_tol/dt;
    e_3 = c*err_tol;
    f_out = f0+f_guess;
    %relNrm = norm(f_out);
    ht_opts.max_rank = max(size(f_guess))^(floor(ndims(f_guess)/2));
    ht_opts.abs_eps = e_1;
    ht_opts.rel_eps = e_1;
    f_out = 0.5*truncate_nonorthog(f_out,ht_opts);
    f_out = orthog(f_out);
    N_half = N(f_out);
    %relNrm = norm(N_half);
    ht_opts.abs_eps = e_2;
    ht_opts.rel_eps = e_2;
    N_half = dt*truncate_nonorthog(N_half,ht_opts);
    f_out = f0 + N_half;
    %relNrm = norm(f_out);
    ht_opts.abs_eps = e_3;
    ht_opts.rel_eps = e_3;
    f_out = truncate_nonorthog(f_out,ht_opts);
    f_out = orthog(f_out);
end
