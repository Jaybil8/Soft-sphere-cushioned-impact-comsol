function out = model
%
% General comsol simulation set-up for the ball impact
%
% Model exported on Mar 28 2023, 22:43 by COMSOL 6.0.0.318.
path = pwd;

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');
[status, msg, msgID] = mkdir('Automated_simulatons')
[status, msg, msgID] = mkdir('Automated_data')
model.modelPath(strcat(path,'\Automated_simulatons'));

prompt = "Specify the value of velocity in meters per second";
velocity = input(prompt,"s");
if str2num(velocity)>10 || str2num(velocity)<1e-5
    error('Enter a value of velocity between 1e-5 and 10 m/s')
end

EHD_label = strcat('EHD_ball_impact_', velocity, 'm_s.mph');

stiffness = input('Type "Soft" or "Hard" to simulate soft or hard materials respectively','s');
switch stiffness
    case 'Soft'
        impactor_density = '1140 [kg/m^3]';
        impactor_E = '250 [kPa]';
    case 'Hard'
        impactor_density = '1358 [kg/m^3]';
        impactor_E = '1100 [kPa]';
    otherwise
        error('Invalid value of the input. Run the program again')
end

% some paramters of the simulation which may need tuning depending on the
% regime studied
mesh_elements_radius = 200;
mesh_elements_axis = 80;
relative_tolerance = 1e-4;

model.label(EHD_label);

model.title('Ball_impact_mediated_by_a_fluid');

model.description('We try to capture the physics and the different regimes of a solid ball impacting on a rigid surface with the presence of a fluid all around');

model.param.set('fluid_density', '1.2 [kg/m^3]', 'density of the fluid');
model.param.set('ball_density', impactor_density, 'density of the ball');
model.param.set('velocity', strcat(velocity, '[m/s]'), 'initial velocity of impact');
model.param.set('ball_radius', '7.1[mm]', 'radius of the ball');
model.param.set('ball_E', impactor_E, 'young modulus of the ball');
model.param.set('ball_nu', '0.47', 'poisson ratio of the ball (almost incompressible)');
model.param.set('fluid_mu', '1.18e-5[Pa*s]', 'fluid viscosity');
model.param.set('delta', '(fluid_mu/(ball_density*velocity*ball_radius))^(1/3)', 'small parameter');
model.param.set('height', 'ball_radius*delta^2', 'height at which the pressure makes a leading order contribution to the deformation');
model.param.set('h0', '1.1524E-5 [m]', 'initial height at the start of the smulation (check 5 times is good enough)');
model.param.set('length', 'ball_radius*delta', 'initial length of deformation in the radial direction');
model.param.set('time_to_impact', 'h0/velocity', 'time it would take the ball to impact if there was no fluid cushioning');
model.param.set('total_time', '10*time_to_impact', 'time of the simulation to capture dynamics till contact');
model.param.set('time_step', '5e-6', 'time_to_impact/30');
model.param.set('shear_modulus', 'ball_E/(2*(1+ball_nu))');
model.param.set('shear_wave_speed', 'sqrt(shear_modulus/ball_density)*0.95');
model.param.set('ratio_impact_wave', 'velocity/shear_wave_speed');
model.param.set('Phi', 'ratio_impact_wave/delta');
model.param.set('pressure_scale', '1E1[Pa]', '1E4[Pa]');
model.param.set('disp_scale', '0.02[mm]');
model.param.set('h_center', 'ball_radius*delta^2');
model.param.set('time_scale', 'delta^2');
model.param.set('c_p', 'sqrt(ball_E/(3*(1-ball_nu)*ball_density))');
model.param.set('phi2', 'velocity/(delta*c_p)');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('evl2', 'Table');

model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').label('Ball');
model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('pos', {'0' 'h0+ball_radius'});
model.component('comp1').geom('geom1').feature('c1').set('rot', 270);
model.component('comp1').geom('geom1').feature('c1').set('r', 'ball_radius');
model.component('comp1').geom('geom1').feature('c1').set('angle', 180);
model.component('comp1').geom('geom1').create('pare1', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare1').setIndex('param', '0.85', 0);
model.component('comp1').geom('geom1').feature('pare1').selection('edge').set('c1(1)', 1);
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').create('pare2', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare2').setIndex('param', '.2', 0);
model.component('comp1').geom('geom1').feature('pare2').selection('edge').set('fin(1)', 1);
model.component('comp1').geom('geom1').run;

model.view.create('view2', 3);
model.view.create('view3', 3);

%define minimum over the leading edge operator
model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl('minop1').selection.geom('geom1', 1);
model.component('comp1').cpl('minop1').selection.set([4]);

model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
model.component('comp1').physics('solid').create('bndl1', 'BoundaryLoad', 1);
model.component('comp1').physics('solid').feature('bndl1').selection.set([4 6]);
model.component('comp1').physics.create('tffs', 'ThinFilmFlowEdge', 'geom1');
model.component('comp1').physics('tffs').selection.set([4 6]);

model.component('comp1').mesh('mesh1').create('fq1', 'FreeQuad');
model.component('comp1').mesh('mesh1').feature('fq1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('fq1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('fq1').create('dis1', 'Distribution');
model.component('comp1').mesh('mesh1').feature('fq1').create('dis2', 'Distribution');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis2').selection.set([4]);

model.result.table('tbl1').comments('Line Minimum 1');
model.result.table('evl2').label('Evaluation 2D');
model.result.table('evl2').comments('Interactive 2D values');

%{
model.component('comp1').view('view1').axis.set('xmin', -0.2765345871448517);
model.component('comp1').view('view1').axis.set('xmax', 2.405282974243164);
model.component('comp1').view('view1').axis.set('ymin', -0.6919764280319214);
model.component('comp1').view('view1').axis.set('ymax', 0.5952959060668945);
%}

model.component('comp1').physics('solid').prop('ShapeProperty').set('order_displacement', 2);
model.component('comp1').physics('solid').prop('EquationForm').set('form', 'Transient');
model.component('comp1').physics('solid').feature('lemm1').set('E_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('E', 'ball_E');
model.component('comp1').physics('solid').feature('lemm1').set('nu_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('nu', 'ball_nu');
model.component('comp1').physics('solid').feature('lemm1').set('rho_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('rho', 'ball_density');
model.component('comp1').physics('solid').feature('init1').set('ut', {'0'; '0'; '-velocity'});
model.component('comp1').physics('solid').feature('init1').label('Initial_velocity');

% Creation of NeoHookean model
model.component('comp1').physics('solid').create('hmm1', 'HyperelasticModel', 2);
model.component('comp1').physics('solid').feature('hmm1').label('Neo_hokean');
model.component('comp1').physics('solid').feature('hmm1').selection.set([1]);
model.component('comp1').physics('solid').feature('hmm1').set('IsotropicOption', 'Enu');
model.component('comp1').physics('solid').feature('hmm1').set('E_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('E', 'ball_E');
model.component('comp1').physics('solid').feature('hmm1').set('nu_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('nu', 'ball_nu');
model.component('comp1').physics('solid').feature('hmm1').set('rho_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('rho', 'ball_density');
model.component('comp1').physics('solid').feature('hmm1').active(false);
% neo-hookean material model is deactivated

model.component('comp1').physics('solid').feature('bndl1').set('FperArea_src', 'root.comp1.tffs.fwallr');
model.component('comp1').physics('solid').feature('bndl1').label('Fluid_pressure_acting_on_ball');
model.component('comp1').physics('tffs').prop('EquationForm').set('form', 'Transient');
model.component('comp1').physics('tffs').prop('ReferencePressure').set('pref', '0[atm]');
model.component('comp1').physics('tffs').feature('ffp1').set('hw1', 'ball_radius + h0 -sqrt(ball_radius^2- r^2)');
model.component('comp1').physics('tffs').feature('ffp1').set('TangentialWallVelocity', 'FromDeformation');
model.component('comp1').physics('tffs').feature('ffp1').set('uw_src', 'root.comp1.u');
model.component('comp1').physics('tffs').feature('ffp1').set('mure_mat', 'userdef');
model.component('comp1').physics('tffs').feature('ffp1').set('mure', 'fluid_mu');
model.component('comp1').physics('tffs').feature('ffp1').set('rho_mat', 'userdef');
model.component('comp1').physics('tffs').feature('ffp1').set('rho', 'fluid_density');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 4);
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis1').label('vertical_refinement');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis1').set('numelem', mesh_elements_axis);
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis2').label('radial_refinement');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis2').set('numelem', mesh_elements_radius);
% Generate mesh
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').create('st1', 'StopCondition');

model.sol('sol1').feature('t1').feature.remove('fcDef');

model.study('std1').feature('time').set('tlist', 'range(0,total_time/400,total_time)');
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', num2str(relative_tolerance));
% always leave geometric non linearity active
model.study('std1').feature('time').set('geometricNonlinearity', true);
model.study('std1').feature('time').set('plot', true);
model.study('std1').feature('time').set('plotfreq', 'tsteps');
model.study('std1').feature('time').set('probesel', 'none');

model.result.dataset.create('rev1', 'Revolve2D');
model.result.numerical.create('min1', 'MinLine');
model.result.numerical('min1').selection.set([4]);
model.result.numerical('min1').set('probetag', 'none');
model.result.create('pg10', 'PlotGroup1D');
model.result.create('pg11', 'PlotGroup1D');
model.study('std1').feature('time').set('plotgroup', 'pg11');
model.result.create('pg13', 'PlotGroup2D');
model.result.create('pg14', 'PlotGroup1D');
model.result.create('pg15', 'PlotGroup1D');
model.result.create('pg16', 'PlotGroup1D');
model.result.create('pg17', 'PlotGroup1D');
model.result.create('pg18', 'PlotGroup2D');
model.result.create('pg19', 'PlotGroup3D');
model.result('pg10').create('lngr1', 'LineGraph');
model.result('pg10').feature('lngr1').set('xdata', 'expr');
model.result('pg10').feature('lngr1').selection.set([4]);
model.result('pg10').feature('lngr1').set('expr', 'z');
model.result('pg11').create('lngr1', 'LineGraph');
model.result('pg11').feature('lngr1').set('xdata', 'expr');
model.result('pg11').feature('lngr1').selection.set([4]);
model.result('pg11').feature('lngr1').set('expr', 'tffs.p');
model.result('pg13').create('surf1', 'Surface');
model.result('pg13').feature('surf1').set('expr', 'residual(solid.disp)');
model.result('pg14').create('lngr1', 'LineGraph');
model.result('pg14').feature('lngr1').set('xdata', 'expr');
model.result('pg14').feature('lngr1').selection.set([4]);
model.result('pg14').feature('lngr1').set('expr', 'residual(tffs.p)');
model.result('pg15').create('tblp1', 'Table');
model.result('pg16').create('lngr1', 'LineGraph');
model.result('pg16').feature('lngr1').set('xdata', 'expr');
model.result('pg16').feature('lngr1').selection.set([4]);
model.result('pg16').feature('lngr1').set('expr', 'z');
model.result('pg17').create('lngr1', 'LineGraph');
model.result('pg17').feature('lngr1').set('xdata', 'expr');
model.result('pg17').feature('lngr1').selection.set([4]);
model.result('pg17').feature('lngr1').set('expr', 'tffs.p');
model.result('pg18').create('surf1', 'Surface');
model.result('pg18').create('con1', 'Contour');
model.result('pg18').create('surf2', 'Surface');
model.result('pg18').create('surf3', 'Surface');
model.result('pg18').feature('con1').set('expr', 'u');
model.result('pg18').feature('surf2').set('expr', 'solid.sz');
model.result('pg18').feature('surf3').set('expr', 'solid.el33');
model.result('pg19').create('vol1', 'Volume');
model.result('pg19').feature('vol1').set('expr', 'z');
model.result.export.create('anim1', 'Animation');
model.result.export.create('anim2', 'Animation');
model.result.export.create('anim3', 'Animation');
model.result.export.create('anim4', 'Animation');


model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol1').feature('st1').set('keeplog', true);
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clist', {'range(0,total_time/400,total_time)' '5.0E-9[s]'});
model.sol('sol1').feature('v1').set('keeplog', true);
model.sol('sol1').feature('v1').feature('comp1_pfilm').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_pfilm').set('scaleval', 'pressure_scale');
model.sol('sol1').feature('v1').feature('comp1_u').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u').set('scaleval', 'disp_scale');
model.sol('sol1').feature('t1').label('Time-Dependent Solver 1.1');
model.sol('sol1').feature('t1').set('tlist', 'range(0,total_time/400,total_time)');
model.sol('sol1').feature('t1').set('rtol', num2str(relative_tolerance));
model.sol('sol1').feature('t1').set('atolglobalfactor', '.05');
model.sol('sol1').feature('t1').set('atolfactor', {'comp1_pfilm' '0.01' 'comp1_u' '0.01'});
model.sol('sol1').feature('t1').set('tstepsbdf', 'manual');
model.sol('sol1').feature('t1').set('endtimeinterpolation', false);
model.sol('sol1').feature('t1').set('timestepbdf', 'time_step');
model.sol('sol1').feature('t1').set('eventtol', 1);
model.sol('sol1').feature('t1').set('stabcntrl', true);
model.sol('sol1').feature('t1').set('rescaleafterinitbw', true);
model.sol('sol1').feature('t1').set('plot', true);
model.sol('sol1').feature('t1').set('plotgroup', 'pg11');
model.sol('sol1').feature('t1').set('plotfreq', 'tsteps');
model.sol('sol1').feature('t1').set('probesel', 'none');
model.sol('sol1').feature('t1').set('keeplog', true);
model.sol('sol1').feature('t1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('t1').feature('aDef').set('storeresidual', 'solvingandoutput');
model.sol('sol1').feature('t1').feature('aDef').set('convinfo', 'detailed');
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('t1').feature('fc1').set('dtech', 'auto');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 25);
model.sol('sol1').feature('t1').feature('fc1').set('termonres', true);

%Stop conditions if gap is smaller than 100nm or timestep 1ns 
model.sol('sol1').feature('t1').feature('st1').label('Stop Condition 1.1');
model.sol('sol1').feature('t1').feature('st1').set('stopcondterminateon', {'true' 'true'});
model.sol('sol1').feature('t1').feature('st1').set('stopcondActive', {'on' 'on'});
model.sol('sol1').feature('t1').feature('st1').set('stopconddesc', {'Stop if times step is too small' 'Stop if gap is smaller than 100nm'});
model.sol('sol1').feature('t1').feature('st1').set('stopcondarr', {'timestep < 1eint(-9' 'comp1.minop1(root.z)<1e-7'});

%Enable progress bar
ModelUtil.showProgress(true);

%Solve command
model.sol('sol1').runAll;

model.result.dataset('rev1').set('revangle', 90);
model.result.numerical('min1').set('table', 'tbl1');
model.result.numerical('min1').set('expr', {'z'});
model.result.numerical('min1').set('unit', {'mm'});
model.result.numerical('min1').set('descr', {'z-coordinate'});
model.result.numerical('min1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result.numerical('min1').set('includepos', true);
model.result.numerical('min1').setResult;
model.result('pg10').label('Tip profile');
model.result('pg10').set('looplevelinput', {'manualindices'});
model.result('pg10').set('looplevelindices', {'range(100,1,105)'});
model.result('pg10').set('titletype', 'manual');
model.result('pg10').set('title', ['Profile during impact\n']);
model.result('pg10').set('xlabel', 'r [m]');
model.result('pg10').set('ylabel', 'z-coordinate (mm)');
model.result('pg10').set('ylog', true);
model.result('pg10').set('xlabelactive', false);
model.result('pg10').set('ylabelactive', false);
model.result('pg10').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg10').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('pg10').feature('lngr1').set('xdataunit', '');
model.result('pg10').feature('lngr1').set('xdatadescractive', true);
model.result('pg10').feature('lngr1').set('xdatadescr', 'r [m]');
model.result('pg10').feature('lngr1').set('linestyle', 'cycle');
model.result('pg10').feature('lngr1').set('linecolor', 'cyclereset');
model.result('pg10').feature('lngr1').set('linewidth', 3);
model.result('pg10').feature('lngr1').set('linemarker', 'cycle');
model.result('pg10').feature('lngr1').set('legend', true);
model.result('pg10').feature('lngr1').set('resolution', 'normal');
model.result('pg11').label('Pressure_profile');
model.result('pg11').set('looplevelinput', {'manual'});
model.result('pg11').set('looplevel', [108]);
model.result('pg11').set('titletype', 'manual');
model.result('pg11').set('titlecolor', 'black');
model.result('pg11').set('title', 'Pressure during impact');
model.result('pg11').set('xlabel', 'r [mm]');
model.result('pg11').set('ylabel', 'Physical pressure (Pa)');
model.result('pg11').set('xlabelactive', false);
model.result('pg11').set('ylabelactive', false);
model.result('pg11').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg11').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('pg11').feature('lngr1').set('xdataunit', '');
model.result('pg11').feature('lngr1').set('xdatadescractive', true);
model.result('pg11').feature('lngr1').set('xdatadescr', 'r [mm]');
model.result('pg11').feature('lngr1').set('linecolor', 'black');
model.result('pg11').feature('lngr1').set('linewidth', 3);
model.result('pg11').feature('lngr1').set('legend', true);
model.result('pg11').feature('lngr1').set('resolution', 'normal');
model.result('pg13').label('Residual');
model.result('pg13').set('looplevel', [1]);
model.result('pg13').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg13').feature('surf1').set('resolution', 'normal');
model.result('pg14').label('Residual Pressure');
model.result('pg14').set('looplevelinput', {'last'});
model.result('pg14').set('xlabel', 'r-coordinate (mm)');
model.result('pg14').set('ylabel', 'residual(tffs.p)');
model.result('pg14').set('xlabelactive', false);
model.result('pg14').set('ylabelactive', false);
model.result('pg14').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg14').feature('lngr1').set('xdataexpr', 'r');
model.result('pg14').feature('lngr1').set('xdatadescr', 'r-coordinate');
model.result('pg14').feature('lngr1').set('resolution', 'normal');
model.result('pg15').label('Radius-time');
model.result('pg15').set('data', 'none');
model.result('pg15').set('titletype', 'manual');
model.result('pg15').set('title', 'Radius of trapped air as a function of time');
model.result('pg15').set('xlabel', 'Time (s)');
model.result('pg15').set('ylabel', 'x (mm)');
model.result('pg15').set('xlabelactive', false);
model.result('pg15').set('ylabelactive', false);
model.result('pg15').feature('tblp1').set('xaxisdata', 1);
model.result('pg15').feature('tblp1').set('plotcolumninput', 'manual');
model.result('pg15').feature('tblp1').set('linewidth', 3);
model.result('pg16').label('Tip profile color');
model.result('pg16').set('looplevelinput', {'interp'});
model.result('pg16').set('titletype', 'manual');
model.result('pg16').set('title', ['Profile during impact\n']);
model.result('pg16').set('xlabel', 'r [m]');
model.result('pg16').set('ylabel', 'z-coordinate (mm)');
model.result('pg16').set('ylog', true);
model.result('pg16').set('xlabelactive', false);
model.result('pg16').set('ylabelactive', false);
model.result('pg16').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg16').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('pg16').feature('lngr1').set('xdataunit', '');
model.result('pg16').feature('lngr1').set('xdatadescractive', true);
model.result('pg16').feature('lngr1').set('xdatadescr', 'r [m]');
model.result('pg16').feature('lngr1').set('linewidth', 3);
model.result('pg16').feature('lngr1').set('legend', true);
model.result('pg16').feature('lngr1').set('resolution', 'normal');
model.result('pg17').label('Pressure_profile color');
model.result('pg17').set('looplevelinput', {'interp'});
model.result('pg17').set('titletype', 'manual');
model.result('pg17').set('titlecolor', 'black');
model.result('pg17').set('title', 'Pressure during impact');
model.result('pg17').set('xlabel', 'r [m]');
model.result('pg17').set('ylabel', 'Physical pressure (Pa)');
model.result('pg17').set('xlabelactive', false);
model.result('pg17').set('ylabelactive', false);
model.result('pg17').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg17').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('pg17').feature('lngr1').set('xdataunit', '');
model.result('pg17').feature('lngr1').set('xdatadescractive', true);
model.result('pg17').feature('lngr1').set('xdatadescr', 'r [m]');
model.result('pg17').feature('lngr1').set('linewidth', 3);
model.result('pg17').feature('lngr1').set('legend', true);
model.result('pg17').feature('lngr1').set('resolution', 'normal');
model.result('pg18').set('looplevel', [1]);
model.result('pg18').set('symmetryaxis', true);
model.result('pg18').set('frametype', 'spatial');
model.result('pg18').feature('surf1').active(false);
model.result('pg18').feature('surf1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.result('pg18').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg18').feature('surf1').set('resolution', 'normal');
model.result('pg18').feature('con1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.result('pg18').feature('con1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg18').feature('con1').set('levelmethod', 'levels');
model.result('pg18').feature('con1').set('levels', '10^{range(0,.1,2)}');
model.result('pg18').feature('con1').set('colortabletrans', 'nonlinear');
model.result('pg18').feature('con1').set('resolution', 'normal');
model.result('pg18').feature('surf2').active(false);
model.result('pg18').feature('surf2').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg18').feature('surf2').set('resolution', 'normal');
model.result('pg18').feature('surf3').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg18').feature('surf3').set('resolution', 'normal');
model.result('pg19').set('looplevel', [1]);
model.result('pg19').set('frametype', 'spatial');
model.result('pg19').feature('vol1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg19').feature('vol1').set('resolution', 'normal');


model.result.export('anim1').label('Pressure_gif');
model.result.export('anim1').set('plotgroup', 'pg17');
model.result.export('anim1').set('giffilename', strcat(path, 'Automated_data', 'tip_profile_', stiffness, '_', velocity, 'ms.gif'));
model.result.export('anim1').set('solnumtype', 'inner');
model.result.export('anim1').set('timeinterp', true);
model.result.export('anim1').set('t', 'range(0, 7.7e-6/20,7.7e-6)');
model.result.export('anim1').set('framesel', 'all');
model.result.export('anim1').set('frametime', 0.3);
model.result.export('anim1').set('fontsize', '20');
model.result.export('anim1').set('colortheme', 'globaltheme');
model.result.export('anim1').set('customcolor', [1 1 1]);
model.result.export('anim1').set('background', 'color');
model.result.export('anim1').set('gltfincludelines', 'on');
model.result.export('anim1').set('title1d', 'on');
model.result.export('anim1').set('legend1d', 'on');
model.result.export('anim1').set('logo1d', 'on');
model.result.export('anim1').set('options1d', 'on');
model.result.export('anim1').set('title2d', 'on');
model.result.export('anim1').set('legend2d', 'on');
model.result.export('anim1').set('logo2d', 'on');
model.result.export('anim1').set('options2d', 'off');
model.result.export('anim1').set('title3d', 'on');
model.result.export('anim1').set('legend3d', 'on');
model.result.export('anim1').set('logo3d', 'on');
model.result.export('anim1').set('options3d', 'off');
model.result.export('anim1').set('axisorientation', 'on');
model.result.export('anim1').set('grid', 'on');
model.result.export('anim1').set('axes1d', 'on');
model.result.export('anim1').set('axes2d', 'on');
model.result.export('anim1').set('showgrid', 'on');
model.result.export('anim2').label('Tip_gif');
model.result.export('anim2').set('plotgroup', 'pg16');
model.result.export('anim2').set('giffilename',strcat(path, 'Automated_data', 'tip_profile_', stiffness, '_', velocity, 'ms.gif'));
model.result.export('anim2').set('solnumtype', 'inner');
model.result.export('anim2').set('timeinterp', true);
model.result.export('anim2').set('t', 'range(0, 3.5e-6/20,3.5e-6)');
model.result.export('anim2').set('framesel', 'all');
model.result.export('anim2').set('frametime', 0.3);
model.result.export('anim2').set('synchronize', false);
model.result.export('anim2').set('fontsize', '20');
model.result.export('anim2').set('colortheme', 'globaltheme');
model.result.export('anim2').set('customcolor', [1 1 1]);
model.result.export('anim2').set('background', 'color');
model.result.export('anim2').set('gltfincludelines', 'on');
model.result.export('anim2').set('title1d', 'on');
model.result.export('anim2').set('legend1d', 'on');
model.result.export('anim2').set('logo1d', 'on');
model.result.export('anim2').set('options1d', 'on');
model.result.export('anim2').set('title2d', 'on');
model.result.export('anim2').set('legend2d', 'on');
model.result.export('anim2').set('logo2d', 'on');
model.result.export('anim2').set('options2d', 'off');
model.result.export('anim2').set('title3d', 'on');
model.result.export('anim2').set('legend3d', 'on');
model.result.export('anim2').set('logo3d', 'on');
model.result.export('anim2').set('options3d', 'off');
model.result.export('anim2').set('axisorientation', 'on');
model.result.export('anim2').set('grid', 'on');
model.result.export('anim2').set('axes1d', 'on');
model.result.export('anim2').set('axes2d', 'on');
model.result.export('anim2').set('showgrid', 'on');
model.result.export('anim3').set('plotgroup', 'pg19');
model.result.export('anim3').set('target', 'player');
model.result.export('anim3').set('solnumtype', 'inner');
model.result.export('anim3').set('solnum', [1]);
model.result.export('anim3').set('maxframes', 7);
model.result.export('anim3').set('showframe', 7);
model.result.export('anim3').set('shownparameter', '1.1827E-5');
model.result.export('anim3').set('fontsize', '20');
model.result.export('anim3').set('colortheme', 'globaltheme');
model.result.export('anim3').set('customcolor', [1 1 1]);
model.result.export('anim3').set('background', 'color');
model.result.export('anim3').set('gltfincludelines', 'on');
model.result.export('anim3').set('title1d', 'on');
model.result.export('anim3').set('legend1d', 'on');
model.result.export('anim3').set('logo1d', 'on');
model.result.export('anim3').set('options1d', 'on');
model.result.export('anim3').set('title2d', 'on');
model.result.export('anim3').set('legend2d', 'on');
model.result.export('anim3').set('logo2d', 'on');
model.result.export('anim3').set('options2d', 'off');
model.result.export('anim3').set('title3d', 'on');
model.result.export('anim3').set('legend3d', 'on');
model.result.export('anim3').set('logo3d', 'on');
model.result.export('anim3').set('options3d', 'off');
model.result.export('anim3').set('axisorientation', 'on');
model.result.export('anim3').set('grid', 'on');
model.result.export('anim3').set('axes1d', 'on');
model.result.export('anim3').set('axes2d', 'on');
model.result.export('anim3').set('showgrid', 'on');
model.result.export('anim4').set('plotgroup', 'pg18');
model.result.export('anim4').set('giffilename', strcat(path, 'Automated_data', 'deformation_isolines_', stiffness, '_', velocity, 'ms.gif'));
model.result.export('anim4').set('solnumtype', 'inner');
model.result.export('anim4').set('solnum', [1]);
model.result.export('anim4').set('frametime', 0.3);
model.result.export('anim4').set('synchronize', false);
model.result.export('anim4').set('fontsize', '20');
model.result.export('anim4').set('colortheme', 'globaltheme');
model.result.export('anim4').set('customcolor', [1 1 1]);
model.result.export('anim4').set('background', 'color');
model.result.export('anim4').set('gltfincludelines', 'on');
model.result.export('anim4').set('title1d', 'on');
model.result.export('anim4').set('legend1d', 'on');
model.result.export('anim4').set('logo1d', 'on');
model.result.export('anim4').set('options1d', 'on');
model.result.export('anim4').set('title2d', 'on');
model.result.export('anim4').set('legend2d', 'on');
model.result.export('anim4').set('logo2d', 'on');
model.result.export('anim4').set('options2d', 'off');
model.result.export('anim4').set('title3d', 'on');
model.result.export('anim4').set('legend3d', 'on');
model.result.export('anim4').set('logo3d', 'on');
model.result.export('anim4').set('options3d', 'off');
model.result.export('anim4').set('axisorientation', 'on');
model.result.export('anim4').set('grid', 'on');
model.result.export('anim4').set('axes1d', 'on');
model.result.export('anim4').set('axes2d', 'on');
model.result.export('anim4').set('showgrid', 'on');

%Save file
mphsave(save,strcat(path,'\Automated_simulatons\', stiffness, '_', velocity))

out = model;
