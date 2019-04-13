function plotState(fig, tspan, sigma, omega)
    global T;
    x0 = 10;
    y0 = 10;
    width=1800;
    height=900;
    set(fig,'units','points','position',[x0,y0,width,height])
    set(fig, 'Visible', 'on')
    set(fig, 'defaultaxesfontsize', 16)

    subplot(2,1,1)
    xlabel('Time t [s]')
    ylabel('MRP \sigma')
    hold on
    grid minor
    plot(tspan/T, sigma(1, :), 'Linewidth', 1.5)
    plot(tspan/T, sigma(2, :), 'Linewidth', 1.5)
    plot(tspan/T, sigma(3, :), 'Linewidth', 1.5)
    legend('\sigma_1', '\sigma_2', '\sigma_3', 'location', 'eastoutside');

    subplot(2,1,2)
    xlabel('Time t [s]')
    ylabel('Angular Velocity \omega [deg/s]')
    hold on
    grid minor
    plot(tspan/T, rad2deg(omega(1, :)), 'Linewidth', 1.5)
    plot(tspan/T, rad2deg(omega(2, :)), 'Linewidth', 1.5)
    plot(tspan/T, rad2deg(omega(3, :)), 'Linewidth', 1.5)
    % threshold
    plot(tspan/T, 3*ones(1,length(tspan)),'k')
    plot(tspan/T, -3*ones(1,length(tspan)),'k')
    leg = legend('\omega_1', '\omega_2', '\omega_3', '3 deg/s threshold', 'location', 'eastoutside');
end