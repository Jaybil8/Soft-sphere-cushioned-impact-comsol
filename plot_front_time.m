vel_list = [0.050881433
        0.100894492
        0.245026071
        0.398172834
        0.491437564
        0.743051348
        0.98612192
        1.239968288
        1.432426556
        1.6
        2.034619632
        2.430985452
        3.375204002
        4.528002089
        5.424159102
        ];

path = pwd;
color_map = plasma;
color_len = floor(length(color_map)/1.2);
n_vel = length(vel_list);
color_idx2 = color_len;

mesh_size = 10e-6;
radius = 7.1e-3;
rho_i = 1140;
nu_i = 0.47;
mu_m = 1.81e-5;
rho_m = 1200; %wrong but it does not enter in the simulation
h0 = 12e-6;
impactor_E =250e3;
steps_before_impact = 60;

count = 0;
fig = figure(1);
fontsize = 14;
for i=1:n_vel
    % soft impactor
    velocity = vel_list(i);
    stk = rho_i*velocity*radius/(12*mu_m);
    G = impactor_E/(2*(1+nu_i));
    cs = sqrt(G/rho_i);
    phi = velocity*stk^(1/3)/cs;
    
    color_idx2 = color_idx2 - floor(color_len/n_vel);
    if mod(count,2)==0
        [tV_r, r_R] = compute_nearest_point_contact(velocity,mesh_size, radius, rho_i, nu_i, mu_m, rho_m, h0, impactor_E, steps_before_impact);
        
        plot(tV_r, r_R, 'LineWidth', 2, 'DisplayName',strcat('$\phi=', num2str(phi, '%.1f'),'$' ), 'color', color_map(color_idx2, :))
        xlabel('$t V/R$', 'Interpreter','latex', 'FontSize',fontsize)
        ylabel('$r^*/R$', 'Interpreter','latex', 'FontSize',fontsize)
        hold on
        x2 = linspace(min(tV_r), max(tV_r));
        a2 = (r_R(end)-r_R(1))/(sqrt(tV_r(end))-sqrt(tV_r(1)));
        b2 = r_R(1)-a2*tV_r(1)^.5;
        %p = [sqrt(times2scaled)'  ones(size(r1'))]\r1';
        %plot(x2, a2*x2.^.5+b2, '--','color', color_map(color_idx2, :), 'LineWidth', 2, 'HandleVisibility','off')
        %plot(x2, yp(1)*x2.^.5+yp(2), '--','color', color_map(color_idx2, :), 'LineWidth', 2, 'HandleVisibility','off')
        hold on
        legend('Interpreter','latex', 'Location','best')
        %axes('position',[.1 .775 .25 .25])
        %box on % put box around new pair of axes
        %index_x2 = (x2 < 9) & (x2 > 13); % range of t near perturbation
        axis tight
    end
    count = count +1;
end

plot(x2, (2*x2).^.5, 'k:', 'LineWidth', 2, 'DisplayName', '$\sqrt{2tV/R}$')
plot(x2, x2.^.5, 'k--', 'LineWidth', 1, 'DisplayName', '$\sqrt{tV/R}$')
set(gca, 'FontSize',fontsize)

cd Automated_figures\
exportgraphics(gcf,strcat('mesh_size_', num2str(mesh_size), '_E_', num2str(impactor_E), '_T_', num2str(steps_before_impact), 'radius_time.png'),'Resolution',300, 'BackgroundColor','none') 

cd(path)


function [tV_R, r_R] = compute_nearest_point_contact(velocity,mesh_size, radius, rho_i, nu_i, mu_m, rho_m, h0, impactor_E, steps_before_impact)

    name = strcat('velocity', strrep(num2str(velocity),'.','_'), '_mesh_size', num2str(mesh_size),'_nu_', num2str(nu_i), '_h0_', num2str(h0), '_E_', num2str(impactor_E), '_T_', num2str(steps_before_impact));
    cd Automated_data\
    cd(name)
    
    times = readmatrix('times.txt');
    r_min = readmatrix('r_min.txt'); 
    first = find(r_min>1e-5,1);
    tip_min = readmatrix(strcat('tip', num2str(first),'.txt'));

    h1 = min(tip_min)/10^3; %conversion from mm to m
    times2 = times(first-1:end)-times(first-1);
    r1 = r_min(first-1:end)/radius;
    times2scaled = times2 *velocity/radius;       
    
    cd ..\..
    tV_R =  times2scaled;
    r_R = r1;

end
