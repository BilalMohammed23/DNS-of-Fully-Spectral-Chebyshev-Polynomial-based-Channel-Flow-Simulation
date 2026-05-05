function dt_max = dt_channel(max_u, max_v, dx, dy_min, gamma)
    eps0 = 1e-12;
    CFL = 0.2;
    dt_x = CFL * dx / max(max_u, eps0);
    dt_y = CFL * dy_min / max(max_v, eps0);
    dt_conv = min(dt_x, dt_y);

    dt_diff = 0.5 / ( gamma * (1/dx^2 + 1/dy_min^2) );

    dt_max = min(dt_conv, dt_diff);
end