function sol = gps_solver(par)
    % Dimensionless groups
    sol.lambda = par.E_f * par.t_f / (par.E_c * par.R);
    sol.X      = (1 + par.nu_c) * (3 - 4*par.nu_c);

    % n=0 mode: self-consistent closed form in terms of eps0
    sol.eps0 = par.eps0;
    sol.q0   = -par.nu_c * par.E_c * sol.lambda * par.eps0 / (1 + (1 + par.nu_c)*(1 - 2*par.nu_c)*sol.lambda);

    % Equivalent applied axial stress, back-derived for reference/comparison
    sol.sigma_c = 2*par.nu_c*sol.q0 - par.E_c*par.eps0;

    % n=1 mode (unaffected by the eps0-vs-sigma_c choice: depends only on kappa)
    sol.A        = par.nu_c * par.E_c * sol.lambda * par.kappa / (2*(1 + sol.X*sol.lambda));
    sol.q1_tilde = -2 * sol.A * par.R;

    % Confinement extremes
    sol.fl_ECF  = sol.q0 - sol.q1_tilde;   % theta = pi:  cos(pi) = -1
    sol.fl_tens = sol.q0 + sol.q1_tilde;   % theta = 0:   cos(0)  = +1

    % FRP hoop strains
    sol.eh_axisym = sol.q0 / (par.E_c * sol.lambda);
    sol.eh_bend   = par.nu_c * par.kappa * par.R / (1 + sol.X*sol.lambda);   % amplitude
    sol.eh_max    = sol.eh_axisym + sol.eh_bend;  % at theta = pi

    % For comparison: previous (CPS-mixed, incorrect) formula
    sol.q0_PDF      = (par.nu_c * sol.lambda / (1 + (1 - par.nu_c^2)*sol.lambda)) * sol.sigma_c;
    sol.eh_axi_PDF  = par.nu_c * sol.sigma_c / (par.E_c * (1 + (1 - par.nu_c^2)*sol.lambda));
end
