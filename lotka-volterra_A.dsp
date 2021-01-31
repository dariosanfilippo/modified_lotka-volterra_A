// =============================================================================
//      Lotka-Volterra complex oscillator (A)
// ============================================================================= 
//
//  This is a complex double oscillator based on the Lotka-Volterra equations.
//  The model is structurally-stable through hyperbolic tangent function
//  saturators and allows for parameters in unstable ranges to explore 
//  different dynamics. Furthermore, this model includes DC-blockers in the 
//  feedback paths to counterbalance a tendency towards fixed-point attractors 
//  – thus enhancing complex behaviours – and obtain signals suitable for audio.
//  Besides the common parameters in an L-V model, this system includes a
//  saturating threshold determining the positive and negative bounds in the
//  equations, while the output peaks are within the [-1.0; 1.0] range.
//
//  The system can be triggered by an impulse or by a constant of arbitrary
//  values for deterministic and reproducable behaviours. Alternatively,
//  the oscillator can be fed with external inputs to be used as a nonlinear
//  distortion unit.
//
// =============================================================================

import("stdfaust.lib");

declare name "Lotka-Volterra complex oscillator (A)";
declare author "Dario Sanfilippo";
declare copyright "Copyright (C) 2021 Dario Sanfilippo 
    <sanfilippo.dario@gmail.com>";
declare version "1.0";
declare license "GPL v3.0 license";

// Lotka-Volterra differential equations:
//      dx/dt = ax - bxy
//      dy/dt = cxy - gy
//
// Discrete model:
//      x[n] = x[n - 1] + dt(ax[n - 1] - bx[n - 1]y[n - 1])
//      y[n] = y[n - 1] + dt(gx[n - 1]y[n - 1] - cy[n - 1])
lotkavolterra(L, a, b, c, g, dt, x_0, y_0) =    prey_level(out * (x / (L * 2.0))) , 
                                                pred_level(out * (y / (L * 2.0)))
    letrec {
        'x = fi.highpass(1, 10, tanh(L, (x_0 + x + dt * (a * x - b * x * y))));
        'y = fi.highpass(1, 10, tanh(L, (y_0 + y + dt * (g * x * y - c * y))));
    };

// tanh() saturator with adjustable saturating threshold
tanh(l, x) = l * ma.tanh(x / l);

// smoothing function for click-free parameter variations using 
// a one-pole low-pass with a 20-Hz cut-off frequency.
smooth(x) = fi.pole(pole, x * (1.0 - pole))
    with {
        pole = exp(-2.0 * ma.PI * 20.0 / ma.SR);
    };

// GUI parameters
prey_level(x) = attach(x , abs(x) : ba.linear2db : 
    levels_group(hbargraph("[5]Prey[style:dB]", -60, 0)));
pred_level(x) = attach(x , abs(x) : ba.linear2db : 
    levels_group(hbargraph("[6]Predator[style:dB]", -60, 0)));
prey_group(x) = vgroup("[1]Prey", x);
pred_group(x) = vgroup("[2]Predator", x);
global_group(x) = vgroup("[3]Global", x);
levels_group(x) = hgroup("[4]Levels (dB)", x);
a = prey_group(hslider("Growth rate", 4, -10, 10, .000001) : smooth);
b = prey_group(hslider("Interaction parameter", 1, -10, 10, .000001) : smooth);
c = pred_group(hslider("Extinction rate", 2, -10, 10, .000001) : smooth);                
g = pred_group(hslider("Interaction parameter", 1, -10, 10, .000001) : smooth);
dt = global_group(
    hslider("[4]dt (time delta)[scale:exp]", 0.1, 0.000001, 10, .000001) : 
        smooth);
input(x) = global_group(nentry("[3]Input value", 1, 0, 10, .000001) <: 
    _ * impulse + _ * checkbox("[1]Constant inputs") + 
        x * checkbox("[0]External inputs"));
impulse = checkbox("[2]Impulse inputs") <: _ - _' : abs;
limit = global_group(
    hslider("[5]Saturation limit[scale:exp]", 4, 1, 1024, .000001) : smooth);
out = global_group(hslider("[6]Output scaling[scale:exp]", 0, 0, 1, .000001) : 
    smooth);

process(x1, x2) = lotkavolterra(limit, a, b, c, g, dt, input(x1), input(x2));

