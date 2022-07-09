function [f_out] = imp_euler_contraction_ht(N,f0,f_guess,dt,eps_err)
%    HTucker Implicit Euler contraction mapping by Bram Rodgers
%    Original Draft 29 Oct, 2021
%
%
%    Description of This Function:
%        This function returns an approximate contraction mapping output
%        for Implicit Euler method. It is guaranteed to be at most eps_er
%        in errror from the output of imp_euler_contraction.m
%
%	Argument List:
%        N       : A (possibly nonlinear) time independent vector field.
%                        It must be able to input and output in TT format.
%        f0      : The current state of your ODE in HTucker format.
%        f_guess : The guess for the next state of your ODE
%                        in HTucker format.
%        dt      : The temporal step size.
%        eps_err : An error to achieve so that the output is at most eps_errr
%                        away from implicit_euler.m
%
%    Return List:
%        f_out   : A refined estimate of the next state guess in Tensor
%                        Train format.
    a = 0.5;                            %A number between in (0,1)
    b = 1-a;                            %The convex complement of that number
    f_out = N(f_guess);                 %Output of vector field
    ht_opts.max_rank = max(size(f_guess))^(floor(ndims(f_guess)/2));
    ht_opts.abs_eps = a*eps_err/(dt);
    ht_opts.rel_eps = a*eps_err/(dt);
    f_out = truncate_std(f_out,ht_opts);  %Approximate vector field
    f_out = f0 + dt*f_out;              %Update guess
    ht_opts.abs_eps = b*eps_err;
    ht_opts.rel_eps = b*eps_err;
    f_out = truncate_std(f_out,ht_opts);     %Approximate guess to low rank.
end
