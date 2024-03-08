function fig = plot_last_step(velocity,mesh_size, radius, rho_i, nu_i, mu_m, rho_m, h0, impactor_E, steps_before_impact)
% This function plots data from comsol simulationscolor_map = plasma;

path = pwd;
name = strcat('velocity', strrep(num2str(velocity),'.','_'), '_mesh_size', num2str(mesh_size),'_nu_', num2str(nu_i), '_h0_', num2str(h0), '_E_', num2str(impactor_E), '_T_', num2str(steps_before_impact));

%color_map = readmatrix("cividis.txt");

cd Automated_data
cd(name)
times = readmatrix('times.txt');
r_min = readmatrix('r_min.txt'); 
first = find(r_min>1e-5,1);
tip_min = readmatrix(strcat('tip', num2str(first),'.txt'));
r0 = r_min(end);
last_step = length(times);
final_profile = readmatrix(strcat('tip', num2str(length(r_min)),'.txt'));
h_final = min(final_profile);
T_impact = h0/velocity; % time it would take to make contact without cushioning
fig = figure();  
last_idx = num2str(last_step);
coor_name = strcat('coor', last_idx, '.txt');
tip_name = strcat('tip', last_idx, '.txt');
pres_name = strcat('pressure', last_idx, '.txt');
coor = readmatrix(coor_name);
coor = 10^3*coor; % units in mm
tip = readmatrix(tip_name);
tip = 10^3*tip; %units in um
pres = readmatrix(pres_name);
m = length(coor);
ax1 = subplot(2,1,1);
tt = num2str((times(last_step)-T_impact), '%.2e');
plot(coor, pres(1:m), 'LineWidth',2.0, 'DisplayName', strcat('t=', tt,' [s]'))
ylabel('$p_m \, [Pa]$','Interpreter','latex')
xticks([])
hold on
ax2 = subplot(2,1,2);
plot(coor, tip, 'LineWidth',2.0)
ylabel('$h \, [\mu m]$', 'Interpreter','latex')
xlabel('$r \, [mm]$','Interpreter','latex')
ylim([-1e-1, 10])
ax1.Position(2)=(ax2.Position(2)+ax2.Position(4))*1.1;
linkaxes([ax1 ax2],'x')
legend(ax1)
fontsize(fig, 18, 'points')
set(fig, 'Position', [488   342   560   840])
cd ..\..
mkdir Automated_figures
cd Automated_figures
exportgraphics(gcf,strcat('End_',name,'.png'),'Resolution',300, 'BackgroundColor','none') 
%saveas(gcf, )
close
cd(path) %ruturn to path from where script was launched
        
        %{
           
           
           color_idx2 = color_idx2 - floor(color_len/n_vel(1));
           if mod(count,2)==0
               
            fig2 = figure(101);
            %cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_data_nh\'
            cd ..\Automated_data_nh\
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
%{
xlim([0.1 25])
legend('Location', 'best','Interpreter','latex')
fontsize(fig2, 12, 'points')
set(fig2, 'Position', [488   342   600   500])
cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_nh\'
exportgraphics(gcf, strcat('radius_evolution_', stiffness, '_',num2str(vel) ,'.png'), 'Resolution',300, 'BackgroundColor','none')
close

writematrix([phi_hard'; r0ell(:,1)'], 'Hard_plate_r0_phi')
writematrix([phi_soft'; r0ell(:,2)'], 'Soft_plate_r0_phi')
end
%}
%{
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