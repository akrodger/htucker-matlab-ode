# htucker-matlab-ode
A toolbox for solving initial value problems using HTucker Tensors.

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

