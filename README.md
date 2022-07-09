# htucker-matlab-ode
A toolbox for solving initial value problems (called IVPs or ODEs) using
HTucker Tensors.

# What is this?
This is the experimental development code for a sequence of papers on a class
of ODE solving algorithms named Step Truncation Methods. These algorithms aim
to mitigate the curse of dimensionality in solutions to initial value problems
by replacing an exponential storage cost with a polynomial one whenever the
solution has quickly decaying singular values. The algorithms have a distinct
advantage of being convergent with small time step, which machine learning
methods for similar methods typically lack.

Step truncation methods are the work of Abram Rodgers, Alec Dektor, and
Daniele Venturi. Below are links to the preprints of the papers proving many
the important related theorems along with a number of applications to partial
differential equations.

Explicit Time Stepping:
https://arxiv.org/abs/2008.00155

Implicit Time Stepping:
https://arxiv.org/abs/2207.01962

# Required Software
Any version of Matlab compatible with the HTucker toolbox. 

Matlab HTucker Toolbox:
https://www.epfl.ch/labs/anchp/index-html/software/htucker/

# To Install

    1) Place the ht-ode-src file wherever you want on your computer.
    2) Open Matlab and navigate to the "ht-ode-src" folder.
    3) Run the command "addpath(pwd); savepath"

# How do I use it?
The different routines in this library of functions perform different parts of
a single time update. The relevant functions:

Explicit Methods, these don't call a rootfinder:

    An order 2 local approximation of the flow map:
        function [f1] = exp_midpoint_ht(vField,f0,dt,A1,A2,A3)
    An order 3 local approximation of the flow map:
        function [f1] = exp_heun_o3_ht(vField,f0,dt,A1,A2,A3,A4,A5)

Implicit Methods, these call a nonlinear rootfinding method developed in a
the Implicit Methods paper above:

    An order 1 local approximation of the flow map:
        function [f1] = imp_euler_update_ht(N,f0,dt,err_tol)
    An order 2 local approximation of the flow map:
        function [f1] = imp_midpoint_update_ht(N,f0,dt,err_tol)

This code contains a number of internal parameters, some of which may be played
with to change approximate IVP solution.

In order to solve a problem over many time steps, place these routines into
a while loop and halt when a desired final time is hit.

You can also use the matlab "help" command to read the small amount of
documentation for each routine. For further information, read the PDF files
linked above.

