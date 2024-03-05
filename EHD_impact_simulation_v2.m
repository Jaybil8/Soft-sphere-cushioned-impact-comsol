function out = model
%
% General comsol simulation set-up for the ball impact
%
% Model exported on Mar 28 2023, 22:43 by COMSOL 6.0.0.318.
path = pwd;

%global tip_data pres_data

import com.comsol.model.*
import com.comsol.model.util.*

vel_list_hard =[0.061799254
        0.096527445
        0.259345269
        0.396189685
        0.498723812
        0.746407497
        0.990023249
        1.237537657
        1.438670225
        1.6
        2.01449168
        2.462853002
        3.308041441
        4.502934572
        5.399556706];

vel_list_soft = [0.050881433
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
%vel_list_soft = [1.6];

mesh_size = 20e-6;
vertical_ref = 0.05;
radial_ref = 0.85;
radius = 7.1e-3;
fluid_mu = 1.81e-5;
h0 = 1.2e-5;
relative_tolerance = 1e-5;
neohookean = 1; % 0 to deactivate neohookean
n_simulations_hard = length(vel_list_soft);
n_simulations_soft = length(vel_list_soft);

for jj = [1,2]
    if jj == 1
        stiffness = 'Hard';
        n_simulations = n_simulations_hard;
        vel_list = vel_list_soft;
    else
        stiffness = 'Soft';
        n_simulations = n_simulations_soft;
        vel_list = vel_list_soft;
    end

for kk = 1:n_simulations
    V = vel_list(kk);
    velocity = num2str(V);

%{
prompt = "Specify the value of velocity in meters per second";
velocity = input(prompt,"s");
if str2num(velocity)>10 || str2num(velocity)<1e-5
    error('Enter a value of velocity between 1e-5 and 10 m/s')
end

stiffness = input('Type "Soft" or "Hard" to simulate soft or hard materials respectively','s');
%}
switch stiffness
    case 'Soft'
        impactor_density = '1140 [kg/m^3]';
        impactor_E = '250 [kPa]';
        c_s = 8.3;
    case 'Hard'
        impactor_density = '1140 [kg/m^3]';
        impactor_E = '1100 [kPa]';
        c_s = 17.5;
    otherwise
        error('Invalid value of the input. Run the program again')
end

model_name = strcat('EHD_ball_impact_', stiffness, velocity, 'ms.mph');

% some paramters of the simulation which may need tuning depending on the
% regime studied


%Decrease tolerance for more accurate, but more computationally expensive
%results

model = ModelUtil.create(model_name);

model.label(model_name);

model.title('Ball_impact_mediated_by_a_fluid');

model.description('We try to capture the physics and the different regimes of a solid ball impacting on a rigid surface with the presence of a fluid all around');

mkdir('Automated_simulations_lower_drop');
mkdir('Automated_data_lower_drop');
model.modelPath(strcat(path,'\Automated_simulatons_lower_drop'));
name = strcat(stiffness, "_", strrep(velocity,'.','_'));

model.param.set('fluid_density', '1.2 [kg/m^3]', 'density of the fluid');
model.param.set('ball_density', impactor_density, 'density of the ball');
model.param.set('velocity', strcat(velocity, '[m/s]'), 'initial velocity of impact');
model.param.set('ball_radius', strcat(num2str(radius), '[m]'), 'radius of the ball');
model.param.set('ball_E', impactor_E, 'young modulus of the ball');
model.param.set('ball_nu', '0.47', 'poisson ratio of the ball (almost incompressible)');
model.param.set('fluid_mu', strcat(num2str(fluid_mu),'[Pa*s]'), 'fluid viscosity');
model.param.set('mesh_size', strcat(num2str(mesh_size), '[m]'), 'minimum size of edge');
model.param.set('delta', '(fluid_mu/(ball_density*velocity*ball_radius))^(1/3)', 'small parameter');
model.param.set('height', 'ball_radius*delta^2', 'height at which the pressure makes a leading order contribution to the deformation');
model.param.set('h0', '.625E-5 [m]', 'initial height at the start of the smulation (check 5 times is good enough)');
model.param.set('length', 'ball_radius*delta', 'initial length of deformation in the radial direction');
model.param.set('time_to_impact', 'h0/velocity', 'time it would take the ball to impact if there was no fluid cushioning');
model.param.set('total_time', '20*time_to_impact', 'time of the simulation to capture dynamics till contact');
% change to have good cfl if wave speed smaller than bc speed
model.param.set('CFL', '0.1', 'CFL number to ensure time accuracy (stability is already provided by the implicit solver)');
% if(Phi>1, CFL*mesh_size/(velocity/delta), CFL*mesh_size/shear_wave_speed)
model.param.set('time_step', 'time_to_impact/7.5', 'time_to_impact/10');
model.param.set('shear_modulus', 'ball_E/(2*(1+ball_nu))');
model.param.set('shear_wave_speed', 'sqrt(shear_modulus/ball_density)*0.95');
model.param.set('ratio_impact_wave', 'velocity/shear_wave_speed');
model.param.set('Phi', 'ratio_impact_wave/delta');
model.param.set('pressure_scale', 'ball_density*velocity^2/delta'); %pressure scale from the scaling
model.param.set('disp_scale', '0.02[mm]');
model.param.set('h_center', 'ball_radius*delta^2');
model.param.set('time_scale', 'delta^2');
model.param.set('c_p', 'sqrt(ball_E/(3*(1-ball_nu)*ball_density))');
model.param.set('phi2', 'velocity/(delta*c_p)');
model.param.set('vertical_ref', vertical_ref, 'minimum size of edge');
model.param.set('radial_ref', radial_ref, 'minimum size of edge');

comp1 = model.component.create('comp1', true);

geom1 = comp1.geom.create('geom1', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('evl2', 'Table');

geom1.axisymmetric(true);
geom1.label('Ball');
geom1.lengthUnit('mm');
geom1.create('c1', 'Circle');
geom1.feature('c1').set('pos', {'0' 'h0+ball_radius'});
geom1.feature('c1').set('rot', 270);
geom1.feature('c1').set('r', 'ball_radius');
geom1.feature('c1').set('angle', 180);
geom1.create('pare1', 'PartitionEdges');
geom1.feature('pare1').setIndex('param', num2str(radial_ref), 0);
geom1.feature('pare1').selection('edge').set('c1(1)', 1);
geom1.run('fin');
geom1.create('pare2', 'PartitionEdges');
geom1.feature('pare2').setIndex('param', num2str(vertical_ref), 0);
geom1.feature('pare2').selection('edge').set('fin(1)', 1);
geom1.run;

model.view.create('view2', 3);
model.view.create('view3', 3);

%define minimum over the leading edge operator
min_op1 = comp1.cpl.create('minop1', 'Minimum');
min_op1.selection.geom('geom1', 1);
comp1.cpl('minop1').selection.set([4]);

%Defining the physics, solid and fluid mechanics
solid = comp1.physics.create('solid', 'SolidMechanics', 'geom1');
solid.create('bndl1', 'BoundaryLoad', 1);
solid.feature('bndl1').selection.set([4 6]);
film = comp1.physics.create('tffs', 'ThinFilmFlowEdge', 'geom1');
film.selection.set([4 6]);

%meshing
mesh1 = comp1.mesh.create('mesh1');
mesh1.create('fq1', 'FreeQuad');
quad_mesh = mesh1.feature('fq1');
quad_mesh.selection.geom('geom1', 2);
quad_mesh.selection.set([1]);
quad_mesh.create('dis1', 'Distribution');
quad_mesh.create('dis2', 'Distribution');
quad_mesh.feature('dis1').selection.set([1]);
quad_mesh.feature('dis2').selection.set([4]);
mesh1.feature('size').set('hauto', 4);
mesh1.feature('size').set('table', 'cfd');
quad_mesh.feature('dis1').label('vertical_refinement');
quad_mesh.feature('dis1').set('numelem', 'floor(ball_radius*vertical_ref/mesh_size)');
quad_mesh.feature('dis2').label('radial_refinement');
quad_mesh.feature('dis2').set('numelem', 'floor(3.14*ball_radius*(1-radial_ref)/(2*mesh_size))');
% Generate mesh
mesh1.run;

model.result.table('tbl1').comments('Line Minimum 1');
model.result.table('evl2').label('Evaluation 2D');
model.result.table('evl2').comments('Interactive 2D values');

%{
comp1.view('view1').axis.set('xmin', -0.2765345871448517);
comp1.view('view1').axis.set('xmax', 2.405282974243164);
comp1.view('view1').axis.set('ymin', -0.6919764280319214);
comp1.view('view1').axis.set('ymax', 0.5952959060668945);
%}

solid.prop('ShapeProperty').set('order_displacement', 2);
solid.prop('EquationForm').set('form', 'Transient');
%linear elasticity solid
linear_elasticity = solid.feature('lemm1');
linear_elasticity.set('E_mat', 'userdef');
linear_elasticity.set('E', 'ball_E');
linear_elasticity.set('nu_mat', 'userdef');
linear_elasticity.set('nu', 'ball_nu');
linear_elasticity.set('rho_mat', 'userdef');
linear_elasticity.set('rho', 'ball_density');
%initial condition solid
solid.feature('init1').set('ut', {'0'; '0'; '-velocity'});
solid.feature('init1').label('Initial_velocity');

% Creation of NeoHookean model
solid.create('hmm1', 'HyperelasticModel', 2);
solid.feature('hmm1').label('Neo_hokean');
solid.feature('hmm1').selection.set([1]);
solid.feature('hmm1').set('IsotropicOption', 'Enu');
solid.feature('hmm1').set('E_mat', 'userdef');
solid.feature('hmm1').set('E', 'ball_E');
solid.feature('hmm1').set('nu_mat', 'userdef');
solid.feature('hmm1').set('nu', 'ball_nu');
solid.feature('hmm1').set('rho_mat', 'userdef');
solid.feature('hmm1').set('rho', 'ball_density');

if neohookean == 0
    solid.feature('hmm1').active(false);
end
% neo-hookean material model is deactivated

% Fluid pressure BC on the solid
solid.feature('bndl1').set('FperArea_src', 'root.comp1.tffs.fwallr');
solid.feature('bndl1').label('Fluid_pressure_acting_on_ball');

film.prop('EquationForm').set('form', 'Transient');
film.prop('ReferencePressure').set('pref', '0[atm]');
thin_film = film.feature('ffp1');
thin_film.set('hw1', 'ball_radius + h0 -sqrt(ball_radius^2- r^2)');
thin_film.set('TangentialWallVelocity', 'FromDeformation');
thin_film.set('uw_src', 'root.comp1.u');
thin_film.set('mure_mat', 'userdef');
thin_film.set('mure', 'fluid_mu');
thin_film.set('rho_mat', 'userdef');
thin_film.set('rho', 'fluid_density');

study1 = model.study.create('std1');
study1.create('time', 'Transient');

sol1 = model.sol.create('sol1');
sol1.study('std1');
sol1.attach('std1');
sol1.create('st1', 'StudyStep');
sol1.create('v1', 'Variables');
sol1.create('t1', 'Time');
sol1.feature('t1').create('fc1', 'FullyCoupled');
sol1.feature('t1').create('st1', 'StopCondition');

sol1.feature('t1').feature.remove('fcDef');

study1.feature('time').set('tlist', 'range(0,total_time/400,total_time)');
study1.feature('time').set('usertol', true);
study1.feature('time').set('rtol', num2str(relative_tolerance));
% always leave geometric non linearity active
study1.feature('time').set('geometricNonlinearity', true);
study1.feature('time').set('plot', true);
study1.feature('time').set('plotfreq', 'tsteps');
study1.feature('time').set('probesel', 'none');

% creating datasets for 2D plots and evaluation of the minimum radius
% during deformation
model.result.dataset.create('rev1', 'Revolve2D');
model.result.numerical.create('min1', 'MinLine');
model.result.numerical('min1').selection.set([4]);
model.result.numerical('min1').set('probetag', 'none');

%Enable progress bar
ModelUtil.showProgress(true);

sol1.attach('std1');
sol1.feature('st1').label('Compile Equations: Time Dependent');
sol1.feature('st1').set('keeplog', true);
sol1.feature('v1').label('Dependent Variables 1.1');
sol1.feature('v1').set('clist', {'range(0,total_time/400,total_time)' '5.0E-9[s]'});
sol1.feature('v1').set('keeplog', true);
sol1.feature('v1').feature('comp1_pfilm').set('scalemethod', 'manual');
sol1.feature('v1').feature('comp1_pfilm').set('scaleval', 'pressure_scale');
sol1.feature('v1').feature('comp1_pfilm').set('resscalemethod', 'manual');
sol1.feature('v1').feature('comp1_pfilm').set('resscaleval', 'pressure_scale');
sol1.feature('v1').feature('comp1_u').set('scalemethod', 'manual');
sol1.feature('v1').feature('comp1_u').set('scaleval', 'disp_scale');
sol1.feature('v1').feature('comp1_u').set('resscalemethod', 'manual');
sol1.feature('v1').feature('comp1_u').set('resscaleval', 'disp_scale');
sol1.feature('t1').label('Time-Dependent Solver 1.1');
sol1.feature('t1').set('tlist', 'range(0,total_time/400,total_time)');
sol1.feature('t1').set('rtol', num2str(relative_tolerance));
sol1.feature('t1').set('atolglobalfactor', '.05');
sol1.feature('t1').set('atolfactor', {'comp1_pfilm' '0.01' 'comp1_u' '0.01'});
sol1.feature('t1').set('tstepsbdf', 'manual');
sol1.feature('t1').set('endtimeinterpolation', true);
sol1.feature('t1').set('timestepbdf', 'time_step');
sol1.feature('t1').set('eventtol', 1);
sol1.feature('t1').set('stabcntrl', true);
sol1.feature('t1').set('rescaleafterinitbw', true);
sol1.feature('t1').set('keeplog', true);
sol1.feature('t1').feature('dDef').label('Direct 1');
sol1.feature('t1').feature('dDef').set('linsolver', 'pardiso');
sol1.feature('t1').feature('aDef').label('Advanced 1');
sol1.feature('t1').feature('aDef').set('storeresidual', 'solvingandoutput');
sol1.feature('t1').feature('aDef').set('convinfo', 'detailed');
sol1.feature('t1').feature('aDef').set('cachepattern', true);
sol1.feature('t1').feature('fc1').label('Fully Coupled 1.1');
sol1.feature('t1').feature('fc1').set('dtech', 'auto');
sol1.feature('t1').feature('fc1').set('maxiter', 25);
sol1.feature('t1').feature('fc1').set('termonres', true);

%Stop conditions if gap is smaller than 100nm or timestep 1ns 
sol1.feature('t1').feature('st1').label('Stop Condition 1.1');
sol1.feature('t1').feature('st1').set('stopcondterminateon', {'true' 'true'});
sol1.feature('t1').feature('st1').set('stopcondActive', {'on' 'on'});
sol1.feature('t1').feature('st1').set('stopconddesc', {'Stop if times step is too small' 'Stop if gap is smaller than 100nm'});
sol1.feature('t1').feature('st1').set('stopcondarr', {'1/timestep > 1e12' 'comp1.minop1(root.z) < 1e-7 [m]'});
sol1.feature('t1').feature('st1').set('storestopcondsol', 'stepafter');
% Save file
mphsave(model_name, strcat("\Automated_simulations_lower_drop\",name, "_ms.mph"))

% Solve command
sol1.runAll;

%-----------------POST-PROCESSING------------------

model.result.dataset('rev1').set('revangle', 90);
model.result.numerical('min1').set('table', 'tbl1');
model.result.numerical('min1').set('expr', {'z'});
model.result.numerical('min1').set('unit', {'mm'});
model.result.numerical('min1').set('descr', {'z-coordinate'});
model.result.numerical('min1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result.numerical('min1').set('includepos', true);
model.result.numerical('min1').setResult;

% -----------------PLOTS-----------------------------
tip_profile1 = model.result.create('tip_bw', 'PlotGroup1D');
tip_profile1.label('Tip profile');
tip_profile1.set('looplevelinput', {'all'});
%tip_profile1.set('looplevelindices', {'range(100,1,105)'});
tip_profile1.set('titletype', 'manual');
tip_profile1.set('title', ['Profile during impact']);
tip_profile1.set('xlabel', 'r [m]');
tip_profile1.set('ylabel', 'z-coordinate (mm)');
tip_profile1.set('ylog', true);
tip_profile1.set('xlabelactive', false);
tip_profile1.set('ylabelactive', false);
tip_plot1 = tip_profile1.create('lngr1', 'LineGraph');
tip_plot1.set('xdata', 'expr');
tip_plot1.selection.set([4]);
tip_plot1.set('expr', 'z');
tip_plot1.set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
tip_plot1.set('xdataexpr', 'if(r<1.5[mm],r,none)');
tip_plot1.set('xdataunit', '');
tip_plot1.set('xdatadescractive', true);
tip_plot1.set('xdatadescr', 'r [m]');
tip_plot1.set('linestyle', 'cycle');
tip_plot1.set('linecolor', 'cyclereset');
tip_plot1.set('linewidth', 3);
tip_plot1.set('linemarker', 'cycle');
tip_plot1.set('legend', true);
tip_plot1.set('resolution', 'normal');

data_folder = strcat(path,'\Automated_data_lower_drop\',name);
[status, msg, msgID] = mkdir(data_folder);

tip_data = mphplot(model, 'tip_bw', 'rangenum', 1,'createplot','off');

n_times = length(tip_data{1});
times = zeros(1, n_times);

pressure_profile1 = model.result.create('pressure_bw', 'PlotGroup1D');
pressure_profile1.label('Pressure_profile');
pressure_profile1.set('looplevelinput', {'all'});
%pressure_profile1.set('looplevel', [108]);
pressure_profile1.set('titletype', 'manual');
pressure_profile1.set('titlecolor', 'black');
pressure_profile1.set('title', 'Pressure during impact');
pressure_profile1.set('xlabel', 'r [mm]');
pressure_profile1.set('ylabel', 'Physical pressure (Pa)');
pressure_profile1.set('xlabelactive', false);
pressure_profile1.set('ylabelactive', false);
pressure_plot1 = pressure_profile1.create('lngr1', 'LineGraph');
pressure_plot1.set('xdata', 'expr');
pressure_plot1.selection.set([4]);
pressure_plot1.set('expr', 'tffs.p');
pressure_plot1.set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
pressure_plot1.set('xdataexpr', 'if(r<1.5[mm],r,none)');
pressure_plot1.set('xdataunit', '');
pressure_plot1.set('xdatadescractive', true);
pressure_plot1.set('xdatadescr', 'r [mm]');
pressure_plot1.set('linecolor', 'black');
pressure_plot1.set('linestyle', 'cycle');
pressure_plot1.set('linewidth', 3);
pressure_plot1.set('legend', true);
pressure_plot1.set('resolution', 'normal');
pres_data = mphplot(model, 'pressure_bw', 'rangenum', 1, 'createplot','off');


study1.feature('time').set('plotgroup', 'pressure_bw');
sol1.feature('t1').set('plot', true);
sol1.feature('t1').set('plotgroup', 'pressure_bw');
sol1.feature('t1').set('plotfreq', 'tsteps');
sol1.feature('t1').set('probesel', 'none');
residual_disp = model.result.create('residual_displacement', 'PlotGroup2D');
residual_disp.label('Residual_displacement')
residual_disp.create('surf1', 'Surface');
residual_disp.feature('surf1').set('expr', 'residual(solid.disp)');
residual_disp.label('Residual');
residual_disp.set('looplevel', [1]);
residual_disp.feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
residual_disp.feature('surf1').set('resolution', 'normal');

residual_pres = model.result.create('residual_pressure', 'PlotGroup1D');
residual_pres.label('Residual_pressure')
residual_pres.create('lngr1', 'LineGraph');
residual_pres.feature('lngr1').set('xdata', 'expr');
residual_pres.feature('lngr1').selection.set([4]);
residual_pres.feature('lngr1').set('expr', 'residual(tffs.p)');
residual_pres.label('Residual Pressure');
residual_pres.set('looplevelinput', {'last'});
residual_pres.set('xlabel', 'r-coordinate (mm)');
residual_pres.set('ylabel', 'residual(tffs.p)');
residual_pres.set('xlabelactive', false);
residual_pres.set('ylabelactive', false);
residual_pres.feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
residual_pres.feature('lngr1').set('xdataexpr', 'r');
residual_pres.feature('lngr1').set('xdatadescr', 'r-coordinate');
residual_pres.feature('lngr1').set('resolution', 'normal');

radius_time = model.result.create('Radius_time', 'PlotGroup1D');
radius_time.label('Table_plot radius-time')
radius_time.create('tblp1', 'Table');
radius_time.set('data', 'none');
radius_time.set('titletype', 'manual');
radius_time.set('title', 'Radius of trapped air as a function of time');
radius_time.set('xlabel', 'Time (s)');
radius_time.set('ylabel', 'x (mm)');
radius_time.set('xlabelactive', false);
radius_time.set('ylabelactive', false);
radius_time.feature('tblp1').set('xaxisdata', 1);
radius_time.feature('tblp1').set('plotcolumninput', 'manual');
radius_time.feature('tblp1').set('plotcolumns', [3]);
radius_time.feature('tblp1').set('linewidth', 3);

surf_plot = model.result.create('surf_plot', 'PlotGroup2D');
surf_plot.create('surf1', 'Surface');
surf_plot.create('con1', 'Contour');
surf_plot.create('surf2', 'Surface');
surf_plot.create('surf3', 'Surface');
surf_plot.feature('con1').set('expr', 'solid.disp');
surf_plot.feature('surf2').set('expr', 'solid.sz');
surf_plot.feature('surf3').set('expr', 'solid.el33');
surf_plot.set('looplevel', [1]);
surf_plot.set('symmetryaxis', true);
surf_plot.set('frametype', 'spatial');
surf_plot.feature('surf1').active(false);
surf_plot.feature('surf1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
surf_plot.feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
surf_plot.feature('surf1').set('resolution', 'normal');
surf_plot.feature('con1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
surf_plot.feature('con1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
surf_plot.feature('con1').set('levelmethod', 'levels');
surf_plot.feature('con1').set('levels', '10^{range(0,.1,2)}');
surf_plot.feature('con1').set('colortabletrans', 'nonlinear');
surf_plot.feature('con1').set('resolution', 'normal');
surf_plot.feature('surf2').active(false);
surf_plot.feature('surf2').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
surf_plot.feature('surf2').set('resolution', 'normal');
surf_plot.feature('surf3').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
surf_plot.feature('surf3').set('resolution', 'normal');
model.result.create('tip_profile19', 'PlotGroup3D');

%{
plot_3d = model.result().create('3D_plot', 'Volume');
plot_3d.feature('vol1').set('expr', 'z');
plot_3d.set('looplevel', [1]);
plot_3d.set('frametype', 'spatial');
plot_3d.feature('vol1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
plot_3d.feature('vol1').set('resolution', 'normal');
%}

gif1 = model.result.export.create('anim1', 'Animation');
gif2 = model.result.export.create('anim2', 'Animation');
gif3 = model.result.export.create('anim3', 'Animation');
gif4 = model.result.export.create('anim4', 'Animation');

tip_profile2 = model.result.create('tip_profile_color', 'PlotGroup1D');
tip_profile2.create('lngr1', 'LineGraph');
tip_profile2.feature('lngr1').set('xdata', 'expr');
tip_profile2.feature('lngr1').selection.set([4]);
tip_profile2.feature('lngr1').set('expr', 'z');
tip_profile2.label('Tip profile color');
tip_profile2.set('looplevelinput', {'manual'});
tip_profile2.set('looplevelindices', {strcat('range(',num2str(mod(n_times,7)),',',num2str(fix(n_times/7)),',', num2str(fix(n_times)),')')}');
tip_profile2.set('titletype', 'manual');
tip_profile2.set('title', ['Profile during impact']);
tip_profile2.set('xlabel', 'r [m]');
tip_profile2.set('ylabel', 'z-coordinate (mm)');
tip_profile2.set('ylog', true);
tip_profile2.set('xlabelactive', false);
tip_profile2.set('ylabelactive', false);
tip_profile2.feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
tip_profile2.feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
tip_profile2.feature('lngr1').set('xdataunit', '');
tip_profile2.feature('lngr1').set('xdatadescractive', true);
tip_profile2.feature('lngr1').set('xdatadescr', 'r [m]');
tip_profile2.feature('lngr1').set('linewidth', 3);
tip_profile2.feature('lngr1').set('legend', true);
tip_profile2.feature('lngr1').set('resolution', 'normal');

pressure_profile2 = model.result.create('pressure_profile_color', 'PlotGroup1D');
pressure_profile2.create('lngr1', 'LineGraph');
pressure_profile2.feature('lngr1').set('xdata', 'expr');
pressure_profile2.feature('lngr1').selection.set([4]);
pressure_profile2.feature('lngr1').set('expr', 'tffs.p');
pressure_profile2.label('Pressure_profile color');
pressure_profile2.set('looplevelinput', {'manual'});
pressure_profile2.set('looplevelindices', {strcat('range(',num2str(mod(n_times,7)),',',num2str(fix(n_times/7)),',', num2str(fix(n_times)),')')}');
pressure_profile2.set('titletype', 'manual');
pressure_profile2.set('titlecolor', 'black');
pressure_profile2.set('title', 'Pressure during impact');
pressure_profile2.set('xlabel', 'r [m]');
pressure_profile2.set('ylabel', 'Physical pressure (Pa)');
pressure_profile2.set('xlabelactive', false);
pressure_profile2.set('ylabelactive', false);
pressure_profile2.feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
pressure_profile2.feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
pressure_profile2.feature('lngr1').set('xdataunit', '');
pressure_profile2.feature('lngr1').set('xdatadescractive', true);
pressure_profile2.feature('lngr1').set('xdatadescr', 'r [m]');
pressure_profile2.feature('lngr1').set('linewidth', 3);
pressure_profile2.feature('lngr1').set('legend', true);
pressure_profile2.feature('lngr1').set('resolution', 'normal');

%gifs path
data_path = strcat("\Automated_data_lower_drop\", name, '_ms');

% Pressure profile evolution 
gif1.label('Pressure_gif');
gif1.set('plotgroup', 'pressure_profile_color');
gif1.set('giffilename', strcat(data_path, 'pressure.gif'));
gif1.set('solnumtype', 'inner');
gif1.set('timeinterp', true);
gif1.set('t', 'range(0, 7.7e-6/20,7.7e-6)');
gif1.set('framesel', 'all');
gif1.set('frametime', 0.3);
gif1.set('fontsize', '20');
gif1.set('colortheme', 'globaltheme');
gif1.set('customcolor', [1 1 1]);
gif1.set('background', 'color');
gif1.set('gltfincludelines', 'on');
gif1.set('title1d', 'on');
gif1.set('legend1d', 'on');
gif1.set('logo1d', 'on');
gif1.set('options1d', 'on');
gif1.set('title2d', 'on');
gif1.set('legend2d', 'on');
gif1.set('logo2d', 'on');
gif1.set('options2d', 'off');
gif1.set('title3d', 'on');
gif1.set('legend3d', 'on');
gif1.set('logo3d', 'on');
gif1.set('options3d', 'off');
gif1.set('axisorientation', 'on');
gif1.set('grid', 'on');
gif1.set('axes1d', 'on');
gif1.set('axes2d', 'on');
gif1.set('showgrid', 'on');

%Leading edge profile evolution
gif2.label('Tip_gif');
gif2.set('plotgroup', 'tip_profile_color');
gif1.set('giffilename', strcat(data_path, 'tip.gif'));
gif2.set('solnumtype', 'inner');
gif2.set('timeinterp', true);
gif2.set('t', 'range(0, 3.5e-6/20,3.5e-6)');
gif2.set('framesel', 'all');
gif2.set('frametime', 0.3);
gif2.set('synchronize', false);
gif2.set('fontsize', '20');
gif2.set('colortheme', 'globaltheme');
gif2.set('customcolor', [1 1 1]);
gif2.set('background', 'color');
gif2.set('gltfincludelines', 'on');
gif2.set('title1d', 'on');
gif2.set('legend1d', 'on');
gif2.set('logo1d', 'on');
gif2.set('options1d', 'on');
gif2.set('title2d', 'on');
gif2.set('legend2d', 'on');
gif2.set('logo2d', 'on');
gif2.set('options2d', 'off');
gif2.set('title3d', 'on');
gif2.set('legend3d', 'on');
gif2.set('logo3d', 'on');
gif2.set('options3d', 'off');
gif2.set('axisorientation', 'on');
gif2.set('grid', 'on');
gif2.set('axes1d', 'on');
gif2.set('axes2d', 'on');
gif2.set('showgrid', 'on');

gif3.set('plotgroup', 'tip_profile19');
gif3.set('target', 'player');
gif3.set('solnumtype', 'inner');
gif3.set('solnum', [1]);
gif3.set('maxframes', 7);
gif3.set('showframe', 7);
gif3.set('shownparameter', '1.1827E-5');
gif3.set('fontsize', '20');
gif3.set('colortheme', 'globaltheme');
gif3.set('customcolor', [1 1 1]);
gif3.set('background', 'color');
gif3.set('gltfincludelines', 'on');
gif3.set('title1d', 'on');
gif3.set('legend1d', 'on');
gif3.set('logo1d', 'on');
gif3.set('options1d', 'on');
gif3.set('title2d', 'on');
gif3.set('legend2d', 'on');
gif3.set('logo2d', 'on');
gif3.set('options2d', 'off');
gif3.set('title3d', 'on');
gif3.set('legend3d', 'on');
gif3.set('logo3d', 'on');
gif3.set('options3d', 'off');
gif3.set('axisorientation', 'on');
gif3.set('grid', 'on');
gif3.set('axes1d', 'on');
gif3.set('axes2d', 'on');
gif3.set('showgrid', 'off');

%{
gif4.set('plotgroup', 'tip_profile18');
gif1.set('giffilename', strcat(data_path, 'pressure.gif'));
gif4.set('solnumtype', 'inner');
gif4.set('solnum', [1]);
gif4.set('frametime', 0.3);
gif4.set('synchronize', false);
gif4.set('fontsize', '20');
gif4.set('colortheme', 'globaltheme');
gif4.set('customcolor', [1 1 1]);
gif4.set('background', 'color');
gif4.set('gltfincludelines', 'on');
gif4.set('title1d', 'on');
gif4.set('legend1d', 'on');
gif4.set('logo1d', 'on');
gif4.set('options1d', 'on');
gif4.set('title2d', 'on');
gif4.set('legend2d', 'on');
gif4.set('logo2d', 'on');
gif4.set('options2d', 'off');
gif4.set('title3d', 'on');
gif4.set('legend3d', 'on');
gif4.set('logo3d', 'on');
gif4.set('options3d', 'off');
gif4.set('axisorientation', 'on');
gif4.set('grid', 'on');
gif4.set('axes1d', 'on');
gif4.set('axes2d', 'on');
gif4.set('showgrid', 'off');
%}

count = 0;
r_min = zeros(1,n_times);
figure()
for i=1:n_times
    [coor, idx] = sort(tip_data{1}{i}.p);
    times(i) = str2num(tip_data{1}{1}.legend{i}(1:end-2));
    writematrix(coor',strcat(data_folder,'\coor',num2str(i),'.txt'),'Delimiter','tab')
    writematrix(tip_data{1}{i}.d(idx),strcat(data_folder,'\tip',num2str(i),'.txt'),'Delimiter','tab')

    [coor2, idx2] = sort(pres_data{1}{i}.p);
    pres_time = pres_data{1}{i};
    writematrix(coor2',strcat(data_folder,'\coor_pres',num2str(i),'.txt'),'Delimiter','tab')
    writematrix(pres_time.d(idx2),strcat(data_folder,'\pressure',num2str(i),'.txt'),'Delimiter','tab')
    tip2 = tip_data{1}{i}.d(idx);
    r_min(i) = coor(find(tip2==min(tip2),1));
    if mod(count,n_times/7) == 0
        subplot(2,1,1)
        plot(coor,tip2, 'DisplayName',num2str(times(i)))
        hold on
        subplot(2,1,2)
        plot(coor2,pres_data{1}{i}.d(idx2))
        hold on
    end
    count = count +1;
end
subplot(2,1,1)
set(gca, 'YScale', 'log')
legend

writematrix(times',strcat(data_folder,'\times.txt'),'Delimiter','tab')
writematrix(r_min',strcat(data_folder,'\r_min.txt'),'Delimiter','tab')

figure()
plot(times,r_min)

%Save file
mphsave(model_name, strcat("\Automated_simulations_lower_drop\", name, "_ms.mph"))

end
end
out = model;
