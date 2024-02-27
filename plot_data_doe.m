function [outputArg1,outputArg2] = plot_data()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_exp = 15; %30;
phi = round(logspace(-1,1,n_exp),2); %round(logspace(-2,1.5,n_exp),2);
min_interval = 0; %-.5; %variations of plus or minus 50%
max_interval = 0; %.5;
n_param = 6;
rng(3) % setting seed for random number
variaton = min_interval + (max_interval-min_interval).*rand(n_param, n_exp); %uniform distribution of 
vel = zeros(1,n_exp);
ell = zeros(1,n_exp);
poisson = .47;

if poisson == .47
    poisson_value = 'poisson_047_';
else
    poisson_value = '';
end

%phi2 = zeros(1,n_exp); to check computation of phi is correct
for i=1:n_exp
    mu_air = 1.81e-5*(1+variaton(1,i));
    rho_solid = 1140*(1+variaton(2,i));
    rho_air = 1.2*(1+variaton(3,i));
    radius = 7.1e-3*(1+variaton(4,i));
    %poisson = .33*(1+variaton(5,i));
    E_solid = 500e+3*(1+variaton(6,i));
    shear_modulus = E_solid/(2*(1+poisson));
    cs = (shear_modulus/rho_solid)^.5;
    p_wave_modulus = E_solid*(1-poisson)/((1+poisson)*(1-2*poisson));
    vel(i) = (phi(i)*cs)^(5/6)*(mu_air/(radius*rho_solid))^(1/6);
    ell(i) = radius^(4/5)*(mu_air/(rho_solid*vel(i)))^(1/5);
    h0(i) = radius^(2/5)*(mu_air/(rho_solid*vel(i)))^(3/5);
    %phi2(i) = vel(i)^(4/3)*(radius/mu_air)^(1/3)*rho_solid^(5/6)/p_wave_modulus^.5;
end

h0 = 1.2e-5; %initial gap
r0 = zeros(n_exp, 1);
%color_map = readmatrix("cividis.txt");
color_map = plasma;
color_len = floor(length(color_map)/1.2);

count = 0;
color_idx2 = color_len;
    for i= 1:n_exp
        cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_data_doe\'
       
        T_impact = h0/vel(i);
        data_folder = strcat(poisson_value,'phi',  strrep(num2str(phi(i)),'.','_'));
        cd(data_folder)
        times = readmatrix('times.txt');
        r_min = readmatrix('r_min.txt'); 
        first = find(r_min>1e-5,1);
        tip_min = readmatrix(strcat('tip', num2str(first),'.txt'));
        r0(i) = r_min(end);
        color_idx = color_len;
        n = length(times);
        n_plots = 7;
        plot_steps = floor(n/n_plots);
        
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
            tt = num2str(times(j)/T_impact, '%.1f');
            plot(coor, pres(1:m), 'LineWidth',1.5, 'DisplayName', strcat('t=', tt,' T'), 'color', color_map(color_idx, :))
            ylabel('$p_m \, [Pa]$','Interpreter','latex') 
            xticks([])
            hold on
            ax2 = subplot(2,1,2);
            plot(coor, tip, 'LineWidth',1.5, 'color', color_map(color_idx, :))
            ylim([10^-4, 1])
            ylabel('$h \, [mm]$', 'Interpreter','latex')
            xlabel('$r \, [m]$','Interpreter','latex')
            hold on
            color_idx = color_idx - floor(color_len/n_plots);
        end
        set(gca, 'YScale', 'log')
        ax1.Position(2)=(ax2.Position(2)+ax2.Position(4))*1.1;
        linkaxes([ax1 ax2],'x')
        legend(ax1)
        fontsize(fig, 18, 'points')
        set(fig, 'Position', [488   342   560   840])
        cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_doe'
        exportgraphics(gcf,strcat('impact', poisson_value,'phi',num2str(phi(i)) ,'plasma.png'),'Resolution',300, 'BackgroundColor','none') 
        %saveas(gcf, )
        close
        %}
        r0(i) = r0(i)/ell(i);
        color_idx2 = color_idx2 - floor(color_len/n_exp);
           if mod(count,2)==1 
            fig2 = figure(101);
            cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_data_doe\'
            h1 = min(tip_min)/10^3; %conversion from mm to m
            times2 = times(first-1:end)-times(first-1);
            r1 = r_min(first-1:end);
            times2scaled = times2 *vel(i)/h1;
            
            x2 = linspace(min(times2scaled), max(times2scaled));
            a2 = (r1(end)-r1(1))/(sqrt(times2scaled(end))-sqrt(times2scaled(1)));
            b2 = r1(1)-a2*times2scaled(1)^.5;
            plot(times2scaled, r1, 'LineWidth', 3, 'DisplayName',strcat('$\phi_2=', num2str(phi(i), '%.2f'),'$' ), 'color', color_map(color_idx2, :))
            xlabel('$t \cdot V/h_1$', 'Interpreter','latex')
            ylabel('$r^* \, [mm]$', 'Interpreter','latex')
            hold on
            %p = [sqrt(times2scaled)'  ones(size(r1'))]\r1';
            yp = polyfit(sqrt(times2scaled), r1, 1);
            %plot(x2, a2*x2.^.5+b2, '--','color', color_map(color_idx2, :), 'LineWidth', 2, 'HandleVisibility','off')
            plot(x2, yp(1)*x2.^.5+yp(2), '--','color', color_map(color_idx2, :), 'LineWidth', 2, 'HandleVisibility','off')
            %} 
           end
           count = count +1;
        end
        


xlim([0.1 25])
legend('Location', 'best','Interpreter','latex')
fontsize(fig2, 12, 'points')
set(fig2, 'Position', [488   342   600   500])
cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_doe\'
exportgraphics(gcf, strcat('radius_evolution_', poisson_value ,'phi_',num2str(phi(i)) ,'.png'), 'Resolution',300, 'BackgroundColor','none')
close

% SCALED PLOT PHI ELL
cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_doe\'
fig = figure();
loglog(phi, r0, 'diamond','MarkerFaceColor',color_map(floor(color_len/1.2), :),'LineWidth',1.5,  'Color',color_map(1, :))
writematrix([phi;round(r0,2)'], strcat(poisson_value,'r0_phi.txt'))
xlabel('$\phi_p$','Interpreter','latex')
ylabel('$r_{contact}/\ell $','Interpreter','latex')
legend('DoE', 'Location','best')
fontsize(fig, 16, 'points')
cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_figures_doe'
exportgraphics(gcf, strcat('trapped_radius_scaled_phi_ell_doe',poisson_value,'.png'), 'Resolution',300, 'BackgroundColor','none')
%}
end
