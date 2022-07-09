function [f1] = exp_heun_o3_ht(vField,f0,dt,A1,A2,A3,A4,A5)
%
%    Explicit Step-Truncation Heun Order 3 by Bram Rodgers
%    Original Draft 16 Feb, 2022
%
%
%    Description of This Function:
%        Apply the 3rd order Heun method to the problem
%            dfdt = vField(f)
%        to get 1 step in time forward.
%        
%        i.e. f1 = f(dt) + O(dt*dt*dt*dt)
%
%	Argument List:
%        vField:     vector field to take htensor object. no time dependence
%        f0    :     state at currect time
%        dt    :     time step size
%        A1    :     local error coefficient for stage 1
%        A2    :     local error coefficient for stage 2
%        A3    :     local error coefficient for stage 3
%        A4    :     local error coefficient for vector field weighted average
%        A5    :     local error coefficient for final truncation
%    Return List:
%        f1    :     order 3 approximation of f(dt)
%
    htSize = size(f0);
    ht_opts.max_rank = floor(sqrt(prod(htSize)));%exact if square tensor. 
    ht_opts.abs_eps = (dt^3)*A1;
    ht_opts.rel_eps = (dt^3)*A1;
    k1 = vField(f0);
    k1 = truncate_nonorthog(k1,ht_opts);
    k2 = vField(f0+((dt/3.0)*k1));
    ht_opts.abs_eps = (dt^2)*A2;
    ht_opts.rel_eps = (dt^2)*A2;
    k2 = truncate_nonorthog(k2,ht_opts);
    k3 = vField(f0+((2.0*dt/3.0)*k2));
    ht_opts.abs_eps = (dt^3)*A3;
    ht_opts.rel_eps = (dt^3)*A3;
    k3 = truncate_nonorthog(k3,ht_opts);
    k1 = (0.25*k1)+(0.75*k3);
    ht_opts.abs_eps = (dt^3)*A4;
    ht_opts.rel_eps = (dt^3)*A4;
    k1 = truncate_nonorthog(k1,ht_opts);
    f1 = f0 + dt*k1;
    ht_opts.abs_eps = (dt^4)*A5;
    ht_opts.rel_eps = (dt^4)*A5;
    f1 = truncate_nonorthog(f1,ht_opts);
    f1 = orthog(f1);
end
