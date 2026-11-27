function field = gps_field(r, theta, par, sol)
    % n=0 contributions (uniform, axisymmetric)
    sr_n0  = -sol.q0 * ones(size(r));
    sth_n0 = -sol.q0 * ones(size(r));

    % n=1 contributions (Michell A*r^3*cos(theta) stress function) 
    sr_n1  = 2 * sol.A * r .* cos(theta);
    sth_n1 = 6 * sol.A * r .* cos(theta);
    trt_n1 = 2 * sol.A * r .* sin(theta);

    % Total in-plane stresses 
    field.sigma_r  = sr_n0  + sr_n1;
    field.sigma_th = sth_n0 + sth_n1;
    field.tau_rth  = trt_n1;                                                 % n=0 contributes zero shear

    % Axial strain (kinematic input from plane sections) 
    field.eps_z = sol.eps0 + par.kappa * r .* cos(theta);

    % Axial stress (GPS: sigma_z = E_c*eps_z + nu_c*(sigma_r + sigma_th)) 
    field.sigma_z = par.E_c * field.eps_z + par.nu_c * (field.sigma_r + field.sigma_th);

    % Component breakdown for plots
    field.sigma_r_n0  = sr_n0;
    field.sigma_r_n1  = sr_n1;
    field.sigma_th_n0 = sth_n0;
    field.sigma_th_n1 = sth_n1;

    % FRP-interface quantities (only meaningful at r = R)
    field.fl     = -field.sigma_r;     
    field.eh_frp = sol.eh_axisym - sol.eh_bend * cos(theta);
end
