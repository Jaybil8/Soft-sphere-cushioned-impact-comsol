cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\'
color_map = plasma;
color_len = floor(length(color_map)/1.2);
cividis = readmatrix('cividis.txt');
%start_name = {'Soft', 'Hard','DoE', 'poisson_01', 'poisson_02', 'poisson_03', 'poisson_04'};
start_name = {'updated2_phi_ell_Soft', 'updated2_phi_ell_Hard','updated2_phi_poisson_047', 'updated2_phi_poisson_01', 'updated2_phi_poisson_02', 'updated2_phi_poisson_03', 'updated2_phi_poisson_04'};
symbol = {'o', 'x', '<', '-o','-o','-o','-o'};
end_name = '_r0_phi.txt';
%legend_names = {'Soft', 'Hard', 'Perturbation', '\nu=0.1', '\nu=0.2', '\nu=0.3', '\nu=0.4'};
legend_names = {'Soft', 'Hard', '\nu=0.47 (perturbed)', '\nu=0.1', '\nu=0.2', '\nu=0.3', '\nu=0.4'};

fig = figure();
for i=1:7
    data = readmatrix(strcat(start_name{i},end_name));  
    if i<=3
        if i ==1
            loglog(data(1,:), data(2,:), symbol{i}, 'MarkerFaceColor',color_map(floor(color_len/1.2), :), 'LineWidth',1.5, 'Color',color_map(1, :),'DisplayName', legend_names{i})
        elseif i==2
            loglog(data(1,:), data(2,:), symbol{i}, 'LineWidth',1.5, 'Color',color_map(1, :),'DisplayName', legend_names{i})
        else
            loglog(data(1,:), data(2,:), symbol{i}, 'LineWidth',1.5, 'Color',color_map(1, :),'DisplayName', legend_names{i})
        end
            %plot(data(1,:), data(2,:), symbol{i}, 'MarkerFaceColor',color_map(floor(color_len/1.2), :), 'LineWidth',1.5, 'Color',color_map(1, :),'DisplayName', legend_names{i})
    else 
        loglog(data(1,:), data(2,:), symbol{i},'color',cividis(end-i*30,:),'MarkerEdgeColor','k' ,'LineWidth',1.5,'DisplayName', legend_names{i})
        %loglog(data(1,:), data(2,:), symbol{i},'LineWidth',3,'DisplayName', legend_names{i})
    end
    hold on
    
end

%loglog(phi_soft, r0ell(:,2), 'diamond','MarkerFaceColor',color_map(floor(color_len/1.2), :),'LineWidth',1.5,  'Color',color_map(1, :))

xlabel('$\phi$','Interpreter','latex')
ylabel('$r_{contact}/ \ell $','Interpreter','latex')
%xlim([0, 3])
legend('Location','best')
fontsize(fig, 16, 'points')
ax = gca;
leftcolor = cividis(1,:);
rightcolor = cividis(end,:);
p = addgradient(ax,leftcolor, rightcolor);  % Red-green gradient to bottom plot
set(p,'FaceAlpha',.3)  % Make transparent

cd 'C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Updated_figures\'
exportgraphics(gcf, strcat('trapped_radius_updated_phi_ell_poisson.png'), 'Resolution',300, 'BackgroundColor','none')
%}