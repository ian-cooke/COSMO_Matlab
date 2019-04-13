function plotControl(fig, tspan, u, m)
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
    ylabel('Control Torque u [N-m]')
    hold on
    grid minor
    plot(tspan/T, u(1, :), 'Linewidth', 2)
    plot(tspan/T, u(2, :), 'Linewidth', 2)
    plot(tspan/T, u(3, :), 'Linewidth', 2)
    leg = legend('u_1', 'u_2', 'u_3', 'location', 'eastoutside');

    subplot(2,1,2)
    xlabel('Time t [s]')
    ylabel('Magnetic Dipole m [T]')
    hold on
    grid minor
    plot(tspan/T, m(1, :), 'Linewidth', 2)
    plot(tspan/T, m(2, :), 'Linewidth', 2)
    plot(tspan/T, m(3, :), 'Linewidth', 2)
    leg = legend('m_1', 'm_2', 'm_3', 'location', 'eastoutside');
end