% The idea here is to assume we start out with some kind of maxwellian
% distibution. We specify a time. During this time, we send the
% distribution through anumber of Fermi reflections in a shrinking flux
% tube, each one of them accelerating the particle. After the given time,
% we compare the distribution to some measued distribution, minimising
% their difference, and in such a way finding the starting distribution.
% The Fermi process will be executed within the cost function.
% Alternatively, we can also keep the time as a parameter to be optimized.
% Potential parameters would be, initial flux tube length, shortening speed
% (based on observed plasma flow), initial distribution parameters.