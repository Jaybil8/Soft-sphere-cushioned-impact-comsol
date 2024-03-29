function out = model
%
% Soft_0_1_ms.m
%
% Model exported on Mar 6 2024, 15:38 by COMSOL 6.0.0.318.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\matlab_codes\Soft-sphere-cushioned-impact-comsol\Automated_simulations');

model.label('Soft_0_1_ms.mph');

model.title('Soft_sphere_impact_mediated_by_a_fluid');

model.description('We try to capture the physics and the different regimes of a solid sphere impacting on a rigid surface with the presence of a fluid all around');

model.param.set('rho_m', '1200 [kg/m^3]', 'density of the fluid');
model.param.set('rho_i', '1140 [kg/m^3]', 'density of the sphere');
model.param.set('velocity', '0.1[m/s]', 'initial velocity of impact');
model.param.set('radius', '0.0071[m]', 'radius of the ball');
model.param.set('sphere_E', '250000[Pa]', 'young modulus of the ball');
model.param.set('nu_i', '0.47', 'poisson ratio of the ball (almost incompressible)');
model.param.set('mu_m', '1.81e-05[Pa*s]', 'fluid viscosity');
model.param.set('mesh_size', '2e-05[m]', 'minimum size of edge');
model.param.set('h0', '2e-05 [m]', 'initial height at the start of the smulation (check 5 times is good enough)');
model.param.set('time_to_impact', 'h0/velocity', 'time it would take the ball to impact if there was no fluid cushioning');
model.param.set('total_time', '20*time_to_impact', 'time of the simulation to capture dynamics till contact');
model.param.set('time_step', 'time_to_impact/20', 'time step as a function of time to contact');
model.param.set('vertical_ref', '0.05', 'minimum size of edge');
model.param.set('radial_ref', '0.85', 'minimum size of edge');
model.param.set('G_i', 'sphere_E/(2*(1+nu_i))');
model.param.set('c_s', 'sqrt(G_i/rho_i)*0.95', 'Shear wave velocity');
model.param.set('c_p', 'sqrt(sphere_E/(3*(1-nu_i)*rho_i))', 'P-wave velocity in the impactor');
model.param.set('ratio_impact_wave', 'velocity/c_s');
model.param.set('Phi', 'ratio_impact_wave/delta_in', 'transition parameter elastic to inertial');
model.param.set('psi', 'velocity/(delta_in*c_p)', 'transition parameter inertial to solid compressibility');
model.param.set('delta_in', '(12*mu_m/(rho_i*velocity*radius))^(1/3)', 'small parameter inertial regime');
model.param.set('delta_el', '(velocity*12*mu_m/(G_i*radius))^(1/5)', 'small parameter elastic regime');
model.param.set('l_inertial', 'radius*delta_in', 'horizontal scale in the radial direction inertial regime');
model.param.set('l_elastic', 'radius*delta_el', 'horizontal scale in the radial direction inertial regime');
model.param.set('h_inertial', 'radius*delta_in^2', 'Height of dimple in the inertial regime');
model.param.set('h_elastic', 'radius*delta_el^2', 'Height of dimple in the elastic regime');
model.param.set('p_inertial', 'rho_i*velocity^2/delta_in', 'inertial pressure scale');
model.param.set('p_elastic', 'G_i*h_elastic/l_elastic', 'elastic pressure scale');
model.param.set('tau_inertial', 'h_inertial/velocity', 'inertial time scale');
model.param.set('tau_elastic', 'h_elastic/velocity', 'elastic time scale');
model.param.set('elasticity_parameter', '4*(1-nu_i^2)/sphere_E*mu_m*velocity*radius^(3/2)/h0^(5/2)');
model.param.set('disp_scale', '0.02[mm]');
model.param.set('pressure_scale', 'rho_i*velocity^2/delta_in');
model.param.set('p_hertz', '2*sphere_E/(1-nu_i^2)*radius^(1/2)');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('evl2', 'Table');

model.component('comp1').geom('geom1').axisymmetric(true);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').label('Ball');
model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('pos', {'0' 'h0+radius'});
model.component('comp1').geom('geom1').feature('c1').set('rot', 270);
model.component('comp1').geom('geom1').feature('c1').set('r', 'radius');
model.component('comp1').geom('geom1').feature('c1').set('angle', 180);
model.component('comp1').geom('geom1').create('pare1', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare1').setIndex('param', '0.85', 0);
model.component('comp1').geom('geom1').feature('pare1').selection('edge').set('c1(1)', 1);
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').create('pare2', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare2').setIndex('param', '0.05', 0);
model.component('comp1').geom('geom1').feature('pare2').selection('edge').set('fin(1)', 1);
model.component('comp1').geom('geom1').run;

model.view.create('view2', 3);
model.view.create('view3', 3);

model.component('comp1').cpl.create('minop1', 'Minimum');
model.component('comp1').cpl('minop1').selection.geom('geom1', 1);
model.component('comp1').cpl('minop1').selection.set([4]);

model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
model.component('comp1').physics('solid').create('bndl1', 'BoundaryLoad', 1);
model.component('comp1').physics('solid').feature('bndl1').selection.set([4 6]);
model.component('comp1').physics('solid').create('hmm1', 'HyperelasticModel', 2);
model.component('comp1').physics('solid').feature('hmm1').selection.set([1]);
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

model.component('comp1').view('view1').axis.set('xmin', -11.669486999511719);
model.component('comp1').view('view1').axis.set('xmax', 18.769487380981445);
model.component('comp1').view('view1').axis.set('ymin', -0.6979999542236328);
model.component('comp1').view('view1').axis.set('ymax', 14.92199993133545);

model.component('comp1').physics('solid').prop('ShapeProperty').set('order_displacement', 2);
model.component('comp1').physics('solid').prop('EquationForm').set('form', 'Transient');
model.component('comp1').physics('solid').feature('lemm1').set('E_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('E', 'sphere_E');
model.component('comp1').physics('solid').feature('lemm1').set('nu_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('nu', 'nu_i');
model.component('comp1').physics('solid').feature('lemm1').set('rho_mat', 'userdef');
model.component('comp1').physics('solid').feature('lemm1').set('rho', 'rho_i');
model.component('comp1').physics('solid').feature('init1').set('ut', {'0'; '0'; '-velocity'});
model.component('comp1').physics('solid').feature('init1').label('Initial_velocity');
model.component('comp1').physics('solid').feature('bndl1').set('FperArea_src', 'root.comp1.tffs.fwallr');
model.component('comp1').physics('solid').feature('bndl1').label('Fluid_pressure_acting_on_ball');
model.component('comp1').physics('solid').feature('hmm1').set('IsotropicOption', 'Enu');
model.component('comp1').physics('solid').feature('hmm1').set('E_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('E', 'sphere_E');
model.component('comp1').physics('solid').feature('hmm1').set('nu_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('nu', 'nu_i');
model.component('comp1').physics('solid').feature('hmm1').set('rho_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('rho', 'rho_i');
model.component('comp1').physics('solid').feature('hmm1').label('Neo_hokean');
model.component('comp1').physics('tffs').prop('EquationForm').set('form', 'Transient');
model.component('comp1').physics('tffs').prop('ReferencePressure').set('pref', '0[atm]');
model.component('comp1').physics('tffs').feature('ffp1').set('hw1', 'radius + h0 -sqrt(radius^2- r^2)');
model.component('comp1').physics('tffs').feature('ffp1').set('TangentialWallVelocity', 'FromDeformation');
model.component('comp1').physics('tffs').feature('ffp1').set('uw_src', 'root.comp1.u');
model.component('comp1').physics('tffs').feature('ffp1').set('mure_mat', 'userdef');
model.component('comp1').physics('tffs').feature('ffp1').set('mure', 'mu_m');
model.component('comp1').physics('tffs').feature('ffp1').set('rho_mat', 'userdef');
model.component('comp1').physics('tffs').feature('ffp1').set('rho', 'rho_m');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 4);
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis1').label('vertical_refinement');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis1').set('numelem', 'floor(radius*vertical_ref/mesh_size)');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis2').label('radial_refinement');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis2').set('numelem', 'floor(3.14*radius*(1-radial_ref)/(2*mesh_size))');
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

model.result.dataset.create('rev1', 'Revolve2D');
model.result.numerical.create('min1', 'MinLine');
model.result.numerical('min1').selection.set([4]);
model.result.numerical('min1').set('probetag', 'none');
model.result.create('tip_bw', 'PlotGroup1D');
model.result.create('pressure_bw', 'PlotGroup1D');
model.result.create('residual_displacement', 'PlotGroup2D');
model.result.create('residual_pressure', 'PlotGroup1D');
model.result.create('Radius_time', 'PlotGroup1D');
model.result.create('surf_plot', 'PlotGroup2D');
model.result.create('tip_profile19', 'PlotGroup3D');
model.result.create('tip_profile_color', 'PlotGroup1D');
model.result.create('pressure_profile_color', 'PlotGroup1D');
model.result('tip_bw').create('lngr1', 'LineGraph');
model.result('tip_bw').feature('lngr1').set('xdata', 'expr');
model.result('tip_bw').feature('lngr1').selection.set([4]);
model.result('tip_bw').feature('lngr1').set('expr', 'z');
model.result('pressure_bw').create('lngr1', 'LineGraph');
model.result('pressure_bw').feature('lngr1').set('xdata', 'expr');
model.result('pressure_bw').feature('lngr1').selection.set([4]);
model.result('pressure_bw').feature('lngr1').set('expr', 'tffs.p/p_elastic');
model.result('residual_displacement').create('surf1', 'Surface');
model.result('residual_displacement').feature('surf1').set('expr', 'residual(solid.disp)');
model.result('residual_pressure').create('lngr1', 'LineGraph');
model.result('residual_pressure').feature('lngr1').set('xdata', 'expr');
model.result('residual_pressure').feature('lngr1').selection.set([4]);
model.result('residual_pressure').feature('lngr1').set('expr', 'residual(tffs.p)');
model.result('Radius_time').create('tblp1', 'Table');
model.result('surf_plot').create('surf1', 'Surface');
model.result('surf_plot').create('con1', 'Contour');
model.result('surf_plot').create('surf2', 'Surface');
model.result('surf_plot').create('surf3', 'Surface');
model.result('surf_plot').feature('surf2').set('expr', 'solid.sz');
model.result('surf_plot').feature('surf3').set('expr', 'solid.el33');
model.result('tip_profile_color').create('lngr1', 'LineGraph');
model.result('tip_profile_color').feature('lngr1').set('xdata', 'expr');
model.result('tip_profile_color').feature('lngr1').selection.set([4]);
model.result('tip_profile_color').feature('lngr1').set('expr', 'z');
model.result('pressure_profile_color').create('lngr1', 'LineGraph');
model.result('pressure_profile_color').feature('lngr1').set('xdata', 'expr');
model.result('pressure_profile_color').feature('lngr1').selection.set([4]);
model.result('pressure_profile_color').feature('lngr1').set('expr', 'tffs.p/p_elastic');
model.result.export.create('anim1', 'Animation');
model.result.export.create('anim2', 'Animation');
model.result.export.create('anim3', 'Animation');
model.result.export.create('anim4', 'Animation');

model.study('std1').feature('time').set('tlist', 'range(0,total_time/400,total_time)');
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', '1e-05');
model.study('std1').feature('time').set('plot', true);
model.study('std1').feature('time').set('plotgroup', 'pressure_bw');
model.study('std1').feature('time').set('plotfreq', 'tsteps');
model.study('std1').feature('time').set('probesel', 'none');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol1').feature('st1').set('keeplog', true);
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clist', {'range(0,total_time/400,total_time)' '1.0E-5[s]'});
model.sol('sol1').feature('v1').set('keeplog', true);
model.sol('sol1').feature('v1').feature('comp1_pfilm').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_pfilm').set('scaleval', 'pressure_scale');
model.sol('sol1').feature('v1').feature('comp1_pfilm').set('resscalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_pfilm').set('resscaleval', 'pressure_scale');
model.sol('sol1').feature('v1').feature('comp1_u').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u').set('scaleval', 'disp_scale');
model.sol('sol1').feature('v1').feature('comp1_u').set('resscalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_u').set('resscaleval', 'disp_scale');
model.sol('sol1').feature('t1').label('Time-Dependent Solver 1.1');
model.sol('sol1').feature('t1').set('tlist', 'range(0,total_time/400,total_time)');
model.sol('sol1').feature('t1').set('rtol', '1e-05');
model.sol('sol1').feature('t1').set('atolglobalfactor', '.05');
model.sol('sol1').feature('t1').set('atolfactor', {'comp1_pfilm' '0.01' 'comp1_u' '0.01'});
model.sol('sol1').feature('t1').set('tstepsbdf', 'manual');
model.sol('sol1').feature('t1').set('timestepbdf', 'time_step');
model.sol('sol1').feature('t1').set('eventtol', 1);
model.sol('sol1').feature('t1').set('stabcntrl', true);
model.sol('sol1').feature('t1').set('rescaleafterinitbw', true);
model.sol('sol1').feature('t1').set('plot', true);
model.sol('sol1').feature('t1').set('plotgroup', 'pressure_bw');
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
model.sol('sol1').feature('t1').feature('st1').label('Stop Condition 1.1');
model.sol('sol1').feature('t1').feature('st1').set('stopcondterminateon', {'true' 'true'});
model.sol('sol1').feature('t1').feature('st1').set('stopcondActive', {'on' 'on'});
model.sol('sol1').feature('t1').feature('st1').set('stopconddesc', {'Stop if times step is too small' 'Stop if gap is smaller than 100nm'});
model.sol('sol1').feature('t1').feature('st1').set('stopcondarr', {'1/timestep > 1e12' 'comp1.minop1(root.z) < 1e-7 [m]'});
model.sol('sol1').feature('t1').feature('st1').set('storestopcondsol', 'stepafter');
model.sol('sol1').runAll;

model.result.dataset('rev1').set('revangle', 90);
model.result.numerical('min1').set('table', 'tbl1');
model.result.numerical('min1').set('expr', {'z'});
model.result.numerical('min1').set('unit', {'mm'});
model.result.numerical('min1').set('descr', {'z-coordinate'});
model.result.numerical('min1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result.numerical('min1').set('includepos', true);
model.result.numerical('min1').setResult;
model.result('tip_bw').label('Tip profile');
model.result('tip_bw').set('titletype', 'manual');
model.result('tip_bw').set('title', 'Profile during impact');
model.result('tip_bw').set('xlabel', 'r [m]');
model.result('tip_bw').set('ylabel', 'z-coordinate (mm)');
model.result('tip_bw').set('ylog', true);
model.result('tip_bw').set('xlabelactive', false);
model.result('tip_bw').set('ylabelactive', false);
model.result('tip_bw').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('tip_bw').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('tip_bw').feature('lngr1').set('xdataunit', '');
model.result('tip_bw').feature('lngr1').set('xdatadescractive', true);
model.result('tip_bw').feature('lngr1').set('xdatadescr', 'r [m]');
model.result('tip_bw').feature('lngr1').set('linestyle', 'cycle');
model.result('tip_bw').feature('lngr1').set('linecolor', 'cyclereset');
model.result('tip_bw').feature('lngr1').set('linewidth', 3);
model.result('tip_bw').feature('lngr1').set('linemarker', 'cycle');
model.result('tip_bw').feature('lngr1').set('legend', true);
model.result('tip_bw').feature('lngr1').set('resolution', 'normal');
model.result('pressure_bw').label('Pressure_profile');
model.result('pressure_bw').set('looplevelinput', {'last'});
model.result('pressure_bw').set('titletype', 'manual');
model.result('pressure_bw').set('titlecolor', 'black');
model.result('pressure_bw').set('title', 'Pressure during impact');
model.result('pressure_bw').set('xlabel', 'r [mm]');
model.result('pressure_bw').set('ylabel', 'tffs.p/p_elastic (1)');
model.result('pressure_bw').set('xlabelactive', false);
model.result('pressure_bw').set('ylabelactive', false);
model.result('pressure_bw').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pressure_bw').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('pressure_bw').feature('lngr1').set('xdataunit', '');
model.result('pressure_bw').feature('lngr1').set('xdatadescractive', true);
model.result('pressure_bw').feature('lngr1').set('xdatadescr', 'r [mm]');
model.result('pressure_bw').feature('lngr1').set('linestyle', 'cycle');
model.result('pressure_bw').feature('lngr1').set('linecolor', 'black');
model.result('pressure_bw').feature('lngr1').set('linewidth', 3);
model.result('pressure_bw').feature('lngr1').set('legend', true);
model.result('pressure_bw').feature('lngr1').set('resolution', 'normal');
model.result('residual_displacement').label('Residual');
model.result('residual_displacement').set('looplevel', [1]);
model.result('residual_displacement').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('residual_displacement').feature('surf1').set('resolution', 'normal');
model.result('residual_pressure').label('Residual Pressure');
model.result('residual_pressure').set('looplevelinput', {'last'});
model.result('residual_pressure').set('xlabel', 'r-coordinate (mm)');
model.result('residual_pressure').set('ylabel', 'residual(tffs.p)');
model.result('residual_pressure').set('xlabelactive', false);
model.result('residual_pressure').set('ylabelactive', false);
model.result('residual_pressure').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('residual_pressure').feature('lngr1').set('xdataexpr', 'r');
model.result('residual_pressure').feature('lngr1').set('xdatadescr', 'r-coordinate');
model.result('residual_pressure').feature('lngr1').set('resolution', 'normal');
model.result('Radius_time').label('Table_plot radius-time');
model.result('Radius_time').set('data', 'none');
model.result('Radius_time').set('titletype', 'manual');
model.result('Radius_time').set('title', 'Radius of trapped air as a function of time');
model.result('Radius_time').set('xlabel', 'Time (s)');
model.result('Radius_time').set('ylabel', 'x (mm)');
model.result('Radius_time').set('xlabelactive', false);
model.result('Radius_time').set('ylabelactive', false);
model.result('Radius_time').feature('tblp1').set('xaxisdata', 1);
model.result('Radius_time').feature('tblp1').set('plotcolumninput', 'manual');
model.result('Radius_time').feature('tblp1').set('linewidth', 3);
model.result('surf_plot').label('2D Plot Group');
model.result('surf_plot').set('looplevel', [1]);
model.result('surf_plot').set('symmetryaxis', true);
model.result('surf_plot').set('frametype', 'spatial');
model.result('surf_plot').feature('surf1').active(false);
model.result('surf_plot').feature('surf1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.result('surf_plot').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf1').set('resolution', 'normal');
model.result('surf_plot').feature('con1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.result('surf_plot').feature('con1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('con1').set('levelmethod', 'levels');
model.result('surf_plot').feature('con1').set('levels', '10^{range(0,.1,2)}');
model.result('surf_plot').feature('con1').set('colortabletrans', 'nonlinear');
model.result('surf_plot').feature('con1').set('resolution', 'normal');
model.result('surf_plot').feature('surf2').active(false);
model.result('surf_plot').feature('surf2').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf2').set('resolution', 'normal');
model.result('surf_plot').feature('surf3').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf3').set('resolution', 'normal');
model.result('tip_profile_color').label('Tip profile color');
model.result('tip_profile_color').set('looplevelinput', {'manual'});
model.result('tip_profile_color').set('titletype', 'manual');
model.result('tip_profile_color').set('title', 'Profile during impact');
model.result('tip_profile_color').set('xlabel', 'r [m]');
model.result('tip_profile_color').set('ylabel', 'z-coordinate (mm)');
model.result('tip_profile_color').set('ylog', true);
model.result('tip_profile_color').set('xlabelactive', false);
model.result('tip_profile_color').set('ylabelactive', false);
model.result('tip_profile_color').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('tip_profile_color').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('tip_profile_color').feature('lngr1').set('xdataunit', '');
model.result('tip_profile_color').feature('lngr1').set('xdatadescractive', true);
model.result('tip_profile_color').feature('lngr1').set('xdatadescr', 'r [m]');
model.result('tip_profile_color').feature('lngr1').set('linewidth', 3);
model.result('tip_profile_color').feature('lngr1').set('legend', true);
model.result('tip_profile_color').feature('lngr1').set('resolution', 'normal');
model.result('pressure_profile_color').label('Pressure_profile color');
model.result('pressure_profile_color').set('looplevelinput', {'interp'});
model.result('pressure_profile_color').set('interp', {'range(-time_to_impact, tau_elastic,  0.0006286*1.3)'});
model.result('pressure_profile_color').set('titletype', 'manual');
model.result('pressure_profile_color').set('titlecolor', 'black');
model.result('pressure_profile_color').set('title', 'Pressure during impact');
model.result('pressure_profile_color').set('xlabel', 'r [m]');
model.result('pressure_profile_color').set('ylabel', 'tffs.p/p_elastic (1)');
model.result('pressure_profile_color').set('xlabelactive', false);
model.result('pressure_profile_color').set('ylabelactive', false);
model.result('pressure_profile_color').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pressure_profile_color').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r/l_elastic,none)');
model.result('pressure_profile_color').feature('lngr1').set('xdataunit', '');
model.result('pressure_profile_color').feature('lngr1').set('xdatadescractive', true);
model.result('pressure_profile_color').feature('lngr1').set('xdatadescr', 'r [m]');
model.result('pressure_profile_color').feature('lngr1').set('linewidth', 3);
model.result('pressure_profile_color').feature('lngr1').set('legend', true);
model.result('pressure_profile_color').feature('lngr1').set('resolution', 'normal');
model.result.export('anim1').label('Pressure_gif');
model.result.export('anim1').set('plotgroup', 'pressure_profile_color');
model.result.export('anim1').set('giffilename', '\Automated_data_lower_drop\Soft_0_1_mstip.gif');
model.result.export('anim1').set('solnumtype', 'inner');
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
model.result.export('anim1').set('options2d', 'on');
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
model.result.export('anim2').set('plotgroup', 'tip_profile_color');
model.result.export('anim2').set('solnumtype', 'inner');
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
model.result.export('anim2').set('options2d', 'on');
model.result.export('anim2').set('title3d', 'on');
model.result.export('anim2').set('legend3d', 'on');
model.result.export('anim2').set('logo3d', 'on');
model.result.export('anim2').set('options3d', 'off');
model.result.export('anim2').set('axisorientation', 'on');
model.result.export('anim2').set('grid', 'on');
model.result.export('anim2').set('axes1d', 'on');
model.result.export('anim2').set('axes2d', 'on');
model.result.export('anim2').set('showgrid', 'on');
model.result.export('anim3').set('plotgroup', 'tip_profile19');
model.result.export('anim3').set('target', 'player');
model.result.export('anim3').set('solnumtype', 'inner');
model.result.export('anim3').set('solnum', [1]);
model.result.export('anim3').set('maxframes', 7);
model.result.export('anim3').set('showframe', 7);
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
model.result.export('anim3').set('options2d', 'on');
model.result.export('anim3').set('title3d', 'on');
model.result.export('anim3').set('legend3d', 'on');
model.result.export('anim3').set('logo3d', 'on');
model.result.export('anim3').set('options3d', 'off');
model.result.export('anim3').set('axisorientation', 'on');
model.result.export('anim3').set('grid', 'on');
model.result.export('anim3').set('axes1d', 'on');
model.result.export('anim3').set('axes2d', 'on');
model.result.export('anim3').set('showgrid', 'on');
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
model.result.export('anim4').set('options2d', 'on');
model.result.export('anim4').set('title3d', 'on');
model.result.export('anim4').set('legend3d', 'on');
model.result.export('anim4').set('logo3d', 'on');
model.result.export('anim4').set('options3d', 'off');
model.result.export('anim4').set('axisorientation', 'on');
model.result.export('anim4').set('grid', 'on');
model.result.export('anim4').set('axes1d', 'on');
model.result.export('anim4').set('axes2d', 'on');
model.result.export('anim4').set('showgrid', 'on');

out = model;
