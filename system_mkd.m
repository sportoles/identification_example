function [ dxdt ] = system_mkd( t, x, M, K, B, f_0, control_law )
%SYSTEM_MKD Summary of this function goes here
%   Detailed explanation goes here
    
%     M = 10.0;
%     B = 0.1;
%     K = 2.0;
%     f_0 = 0;

    f_claw = control_law(t,x);
    %f_claw = 0.0;

    dxdt_1 = x(2);
    dxdt_2 = -K/M*x(1) - B/M*x(2) + f_0 + f_claw;
    dxdt = [dxdt_1; dxdt_2];

end
