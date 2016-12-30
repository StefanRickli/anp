clear all;
close all;

s = tf('s');
f = (s+1)/(s*(s+2)*(s+3));

args.tf_obj = f;

args.radii.R =              10;
args.radii.auto_main_R =	true;

args.angles.crop_inf_transition =  7;                  % [°], roundoff-parameter of the D-curve, usually not necessary to change
args.angles.min_angle_contribution_at_R = 65;          % [°]
args.angles.detour =        45;                 % [°]

args.weights.pole =  1;
args.weights.zero =  1/3;


a = anp_tf_processor;

a.set_all_params(args)
