function [x_out] = arfpie(F,F_appx,q,metric,x_0,err_tol)
%    Adaptive Rank Fixed Point Iteration Estimate by Bram Rodgers
%    Original Draft 27 Oct, 2021
%
%
%    Description of This Function:
%        This function follows a contraction mapping iteration
%        by approximating the contraction map and then updating
%        in a modified direction. The iteration is
%        
%        y[j+1] = F(x[j])
%        eps[j] = (1-q) * metric(y[j+1],x[j])/b
%        x[j+1] = F_appx(x[j],eps[j])
%        
%        where 0 <= q < 1 is a Lipschitz constant of F
%        and b>2 is a real number.
%        
%        In this immplementation, we set b = 2.1
%
%	Argument List:
%        F:       A function handle for a contraction mapping.
%        F_appx:  A function handle for an approximation of
%                   a contraction satisfying
%                     metric(F(z),F_appx(z,err)) <= err
%        q:       A Lipschitz constant for F smaller than 1.
%        metric:  A function handle for computing distance between
%                 two objects. expected use:
%                     d = metric(x,y)
%                 d is a positive real, and x,y are 2 objects of the same type.
%        x_0:     initial guess of the fixed point
%        err_tol: The stopping tolerance of the contraction.
%
%    Return List:
%        x_out:   An object satisfying:
%                     metric(x_out,F(x_out)) <= err_tol
    b       = 5;                        %An arbitrary parameter over 2
    y_jp1   = F(x_0);                   %Initialize contraction output
    err_now = metric(y_jp1,x_0);        %Get current error
    eps_j   = err_now*(1.0-q)/b;        %Find approximation error
    x_out   = F_appx(x_0,eps_j);        %Update fixed point guess
    c       = 2.0/b;                    %Truncation coefficient
    g_c_q   = (c-(c*q))+q;              %Convergence rate
    err_bdd = err_now*(g_c_q/(1-g_c_q))*((c*((1-q)/(1+q)))+1.0);
    err_chk = err_now;                  %variable for checking error.
    while(err_chk > err_tol)       %Loop until the error is low    
        y_jp1   = F(x_out);                 %Update contraction output
        err_old = err_now;
        err_now = metric(y_jp1,x_out);      %Get current error
        eps_j   = err_now*(1.0-q)/b;        %Find approximation error
        x_out   = F_appx(x_out,eps_j);      %Update fixed point guess
        err_bdd = err_bdd*g_c_q;            %Update proven error bound.
        err_chk = min([err_bdd err_now]);   %Stop early if error is low.
    end
end
