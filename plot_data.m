function [outputArg1,outputArg2] = plot_data()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Experimental_data\'

h0 = 1.2e-5; %initial gap
vel_list_soft = csvread("soft_impact.csv");
vel_list_hard = csvread("rigid_impact.csv");
vel_list_soft_sim = [vel_list_soft(1:9,1); 1.6; vel_list_soft(10:end,1)];
n_vel = size(vel_list_soft_sim);
r0 = zeros(n_vel(1), 2);
h_final = zeros(n_vel(1), 2);
r0ell = zeros(n_vel(1), 2);
%color_map = readmatrix("cividis.txt");
color_map = plasma;
color_len = floor(length(color_map)/1.2);
mu_air = 1.81e-5;
rho = 1140;
radius = 7.1e-3;
poisson = .47;
phi_hard = zeros(n_vel(1), 1);
phi_soft = zeros(n_vel(1), 1);
ell_hard = zeros(n_vel(1), 1);
ell_soft = zeros(n_vel(1), 1);

for ii = 1:2
    if ii ==1
        stiffness = 'Hard';
        E = 1.1e+6;
    else
        stiffness = 'Soft';
        E = 250e+3;
    end
    G_solid = E/(2*(1+poisson));
    cs = sqrt(G_solid/rho);
    cp = sqrt(E*(1-poisson)/((1+poisson)*(1-2*poisson)*rho));
    %ell_factor = radius^(5/6)*(mu_air/rho)^(1/6);
    phi_factor = 1/(cs*(12*mu_air/(rho*radius))^(1/3)); %choose between cs or cp
    ell_factor = (radius^4*12*mu_air/G_solid)^(1/5);

    height_factor = radius^(2/5)*(mu_air/rho)^(3/5);
    time_factor = radius^(3/5)*(mu_air/rho)^(2/5);
    %{
    phi_factor = 1/(3*cs*(3*mu_air/(rho*radius))^(1/3)); %choose between cs or cp
    ell_factor = radius^(2/3)*(3*mu_air/rho)^(1/3);
    height_factor = 3*radius^(1/3)*(mu_air/rho)^(2/3);
    time_factor = radius^(3/5)*(mu_air/rho)^(2/5);
    %}
    count = 0;
    color_idx2 = color_len;
    for i= 1:n_vel(1)
        %replace "nh" with "higher_drop" or "lower_drop"
        cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_data_nh\'
        vel = vel_list_soft_sim(i);
        T_impact = h0/vel;
        data_folder = strcat(stiffness, '_',  strrep(num2str(vel),'.','_'));
        cd(data_folder)
        times = readmatrix('times.txt');
        r_min = readmatrix('r_min.txt'); 
        first = find(r_min>1e-5,1);
        tip_min = readmatrix(strcat('tip', num2str(first),'.txt'));
        r0(i, ii) = r_min(end);
        color_idx = color_len;
        n = length(times);
        n_plots = 7;
        plot_steps = floor(n/n_plots);
        final_profile = readmatrix(strcat('tip', num2str(length(r_min)),'.txt'));
        h_final(i, ii) = min(final_profile);

        ell = ell_factor*vel^(1/5);
        phi = phi_factor*vel^(4/3);
        %ell = ell_factor*vel^(-1/3);
        
        %height_scale = height_factor*vel^(-3/5);
        %pressure_scale = rho^(6/5)*vel^(11/5)*(radius/mu_air)^(1/5);
        height_scale = height_factor*vel^(-2/3);
        pressure_scale = rho^(6/5)*vel^(11/5)*(radius/mu_air)^(1/5);
        %time_scale = time_factor*vel^(-8/5);
        time_scale = time_factor*vel^(-7/5);
        if ii ==1
                phi_hard(i) = phi;
                ell_hard(i) = ell;
                r0ell(i, ii) = r0(i, ii)/ell_hard(i);
            else
                phi_soft(i) = phi;
                ell_soft(i) = ell;
                r0ell(i, ii) = r0(i, ii)/ell_soft(i);
            end
        
            

        %{
        fig = figure();
        for j = mod(n-1,plot_steps)+1 : plot_steps : n
    
            coor_name = strcat('coor', num2str(j), '.txt');
            tip_name = strcat('tip', num2str(j), '.txt');
            pres_name = strcat('pressure', num2str(j), '.txt');
            coor = readmatrix(coor_name);
            tip = readmatrix(tip_name);
            pres = readmatrix(pres_name);
            m = length(coor);
            ax1 = subplot(2,1,1);
            %tt = num2str(times(j)/T_impact, '%.1f');
            tt = num2str((times(j)-T_impact)/time_scale, '%.2f');
            %plot(coor, pres(1:m), 'LineWidth',1.5, 'DisplayName', strcat('t=', tt,' T'), 'color', color_map(color_idx, :))
            plot(coor/ell, pres(1:m)/pressure_scale, 'LineWidth',1.5, 'DisplayName', strcat('t=', tt,' T'), 'color', color_map(color_idx, :))
            %ylabel('$p_m \, [Pa]$','Interpreter','latex') 
            ylabel('$p_m/\hat{p}$','Interpreter','latex')
            xticks([])
            hold on
            ax2 = subplot(2,1,2);
            %plot(coor, tip, 'LineWidth',1.5, 'color', color_map(color_idx, :))
            plot(coor/ell, tip/10^3/height_scale, 'LineWidth',1.5, 'color', color_map(color_idx, :))
            %ylim([10^-4, 1])
            %ylabel('$h \, [mm]$', 'Interpreter','latex')
            %xlabel('$r \, [m]$','Interpreter','latex')
            ylabel('$h / \hat{h}$', 'Interpreter','latex')
            xlabel('$r / \ell$','Interpreter','latex')
            hold on
            color_idx = color_idx - floor(color_len/n_plots);
        end
        set(gca, 'YScale', 'log')
        ax1.Position(2)=(ax2.Position(2)+ax2.Position(4))*1.1;
        linkaxes([ax1 ax2],'x')
        legend(ax1)
        fontsize(fig, 18, 'points')
        set(fig, 'Position', [488   342   560   840])
        cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_regimes'
        exportgraphics(gcf,strcat('normalized_impact_', stiffness, '_',num2str(vel) ,'plasma.png'),'Resolution',300, 'BackgroundColor','none') 
        %saveas(gcf, )
        close
        
        %}
           
           
           color_idx2 = color_idx2 - floor(color_len/n_vel(1));
           if mod(count,2)==0
               %{
            fig2 = figure(101);
            %cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_data_nh\'
            cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_data_nh\'
            h1 = min(tip_min)/10^3; %conversion from mm to m
            times2 = times(first-1:end)-times(first-1);
            r1 = r_min(first-1:end);
            times2scaled = times2 *vel/h1;
            
            x2 = linspace(min(times2scaled), max(times2scaled));
            a2 = (r1(end)-r1(1))/(sqrt(times2scaled(end))-sqrt(times2scaled(1)));
            b2 = r1(1)-a2*times2scaled(1)^.5;
            plot(times2scaled, r1, 'LineWidth', 1, 'DisplayName',strcat('$\phi=', num2str(phi, '%.1f'),'$' ), 'color', color_map(color_idx2, :))
            xlabel('$t \cdot V/h_1$', 'Interpreter','latex')
            ylabel('$r^* \, [mm]$', 'Interpreter','latex')
            hold on
            %p = [sqrt(times2scaled)'  ones(size(r1'))]\r1';
            yp = polyfit(sqrt(times2scaled), r1, 1);
            %plot(x2, a2*x2.^.5+b2, '--','color', color_map(color_idx2, :), 'LineWidth', 2, 'HandleVisibility','off')
            plot(x2, yp(1)*x2.^.5+yp(2), '--','color', color_map(color_idx2, :), 'LineWidth', 2, 'HandleVisibility','off')
            
            %axes('position',[.1 .775 .25 .25])
            %box on % put box around new pair of axes
            %index_x2 = (x2 < 9) & (x2 > 13); % range of t near perturbation
            %axis tight
            %} 
           end
           count = count +1;
        
end
%{
xlim([0.1 25])
legend('Location', 'best','Interpreter','latex')
fontsize(fig2, 12, 'points')
set(fig2, 'Position', [488   342   600   500])
cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh\'
exportgraphics(gcf, strcat('radius_evolution_', stiffness, '_',num2str(vel) ,'.png'), 'Resolution',300, 'BackgroundColor','none')
close
%}
writematrix([phi_hard'; r0ell(:,1)'], 'Hard_plate_r0_phi')
writematrix([phi_soft'; r0ell(:,2)'], 'Soft_plate_r0_phi')
end


%cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh\'
cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh\'
fig = figure();
writematrix([vel_list_soft_sim(:, 1), r0(:,2)*10^3], 'Simulation_soft.csv')
writematrix([vel_list_soft_sim(:, 1), r0(:,1)*10^3], 'Simulation_hard.csv')
writematrix([vel_list_soft_sim(:, 1), h_final(:,2)*10^3], 'Simulation_soft_final_height.csv')
writematrix([vel_list_soft_sim(:, 1), h_final(:,1)*10^3], 'Simulation_hard_final_height.csv')
plot(vel_list_soft(:, 1), vel_list_soft(:,2), 'o','MarkerFaceColor',color_map(1, :),'LineWidth',1.5, 'Color',color_map(1, :))
hold on
plot(vel_list_hard(:, 1), vel_list_hard(:,2), 'x','MarkerFaceColor',color_map(1, :), 'LineWidth',1.5,  'Color',color_map(1, :))
%plot(vel_list_soft_sim(:, 1), r0(:,2)*10^3, 'diamond','MarkerFaceColor',color_map(floor(color_len/1.2), :),'LineWidth',1.5,  'Color',color_map(1, :))
%plot(vel_list_soft_sim(:, 1), r0(:,1)*10^3, 's', 'MarkerFaceColor', color_map(floor(color_len/1.2), :)','LineWidth',1.5,  'Color',color_map(1, :))
plot(vel_list_soft_sim(:, 1), r0(:,2)*10^3, 'o','LineWidth',1.5,  'Color',color_map(end/2, :))
plot(vel_list_soft_sim(:, 1), r0(:,1)*10^3, 'x','LineWidth',1.5,  'Color',color_map(end/2, :))
xlabel('$V \, [m/s]$','Interpreter','latex')
ylabel('$r_{contact} \, [mm]$','Interpreter','latex')
ylim([0, 1.4])
legend('Soft (Experiment)','Hard (Experiment)', 'Soft (Simulation)', 'Hard (Simulation)')
fontsize(fig, 16, 'points')
%cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh'
cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh'
exportgraphics(gcf, strcat('trapped_radius_comparison_color', stiffness,'.png'), 'Resolution',300, 'BackgroundColor','none')

%}
% SCALED PLOT PHI ELL
%{
cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh\'
fig = figure();
%plot(vel_list_soft(:, 1), vel_list_soft(:,2), 'diamond','MarkerFaceColor',color_map(1, :),'LineWidth',1.5, 'Color',color_map(1, :))
%hold on
%plot(vel_list_hard(:, 1), vel_list_hard(:,2), 's','MarkerFaceColor',color_map(1, :), 'LineWidth',1.5,  'Color',color_map(1, :))
loglog(phi_soft, r0ell(:,2), 'diamond','MarkerFaceColor',color_map(floor(color_len/1.2), :),'LineWidth',1.5,  'Color',color_map(1, :))

hold on
loglog(phi_hard, r0ell(:,1), 's', 'MarkerFaceColor', color_map(floor(color_len/1.2), :)','LineWidth',1.5,  'Color',color_map(1, :))
writematrix([phi_hard'; r0ell(:,1)'], 'Hard_plate_r0_phi')
writematrix([phi_soft'; r0ell(:,2)'], 'Soft_plate_r0_phi')
xlabel('$\phi$','Interpreter','latex')
ylabel('$r_{contact}/\ell $','Interpreter','latex')
legend('Soft (Simulation)', 'Hard (Simulation)', 'Location','best')
fontsize(fig, 16, 'points')
cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh'
exportgraphics(gcf, strcat('hicks_trapped_radius_scaled_phi_ell', stiffness,'.png'), 'Resolution',300, 'BackgroundColor','none')
%}
end