function out = trapping_angle(B_ref,theta_ref,B)

theta = asind(sqrt(B.data/B_ref)*sind(theta_ref));

out = irf.ts_scalar(B.time,[theta 180-theta]);