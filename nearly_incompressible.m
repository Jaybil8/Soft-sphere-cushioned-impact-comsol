function out = model
%
% nearly_incompressible.m
%
% Model exported on Jun 19 2023, 08:59 by COMSOL 6.0.0.318.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Comsol_simulations\Automated_simulations_nh');

model.label('Soft_5_4242_ms_incompressible.mph');

model.title('Ball_impact_mediated_by_a_fluid');

model.description('We try to capture the physics and the different regimes of a solid ball impacting on a rigid surface with the presence of a fluid all around');

model.param.set('fluid_density', '1.2 [kg/m^3]', 'density of the fluid');
model.param.set('ball_density', '1140 [kg/m^3]', 'density of the ball');
model.param.set('velocity', '3[m/s]', 'initial velocity of impact');
model.param.set('ball_radius', '0.0071[m]', 'radius of the ball');
model.param.set('ball_E', '250 [kPa]', 'young modulus of the ball');
model.param.set('ball_nu', '0.47', 'poisson ratio of the ball (almost incompressible)');
model.param.set('fluid_mu', '1.81e-05[Pa*s]', 'fluid viscosity');
model.param.set('mesh_size', '2e-05[m]', 'minimum size of edge');
model.param.set('delta', '(fluid_mu/(ball_density*velocity*ball_radius))^(1/3)', 'small parameter');
model.param.set('height', 'ball_radius*delta^2', 'height at which the pressure makes a leading order contribution to the deformation');
model.param.set('h0', '1.2E-5 [m]', 'initial height at the start of the smulation (check 5 times is good enough)');
model.param.set('length', 'ball_radius*delta', 'initial length of deformation in the radial direction');
model.param.set('time_to_impact', 'h0/velocity', 'time it would take the ball to impact if there was no fluid cushioning');
model.param.set('total_time', '10*time_to_impact', 'time of the simulation to capture dynamics till contact');
model.param.set('CFL', '0.1', 'CFL number to ensure time accuracy (stability is already provided by the implicit solver)');
model.param.set('time_step', 'time_to_impact/15', 'time_to_impact/10');
model.param.set('shear_modulus', 'ball_E/(2*(1+ball_nu))');
model.param.set('shear_wave_speed', 'sqrt(shear_modulus/ball_density)*0.95');
model.param.set('ratio_impact_wave', 'velocity/shear_wave_speed');
model.param.set('Phi', 'ratio_impact_wave/delta');
model.param.set('pressure_scale', 'ball_density*velocity^2/delta');
model.param.set('disp_scale', '0.02[mm]');
model.param.set('h_center', 'ball_radius*delta^2');
model.param.set('time_scale', 'delta^2');
model.param.set('c_p', 'sqrt(ball_E/(3*(1-ball_nu)*ball_density))');
model.param.set('phi2', 'velocity/(delta*c_p)');
model.param.set('vertical_ref', '0.05', 'minimum size of edge');
model.param.set('radial_ref', '0.85', 'minimum size of edge');

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

model.component('comp1').view('view1').axis.set('xmin', -0.21034100651741028);
model.component('comp1').view('view1').axis.set('xmax', 0.6175994873046875);
model.component('comp1').view('view1').axis.set('ymin', -0.04186049848794937);
model.component('comp1').view('view1').axis.set('ymax', 0.3212556838989258);

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
model.component('comp1').physics('solid').feature('bndl1').set('FperArea_src', 'root.comp1.tffs.fwallr');
model.component('comp1').physics('solid').feature('bndl1').label('Fluid_pressure_acting_on_ball');
model.component('comp1').physics('solid').feature('hmm1').set('kappa', '200*solid.Gequ');
model.component('comp1').physics('solid').feature('hmm1').set('IsotropicOption', 'Enu');
model.component('comp1').physics('solid').feature('hmm1').set('Compressibility_NeoHookean', 'NearlyIncompressibleQuadratic');
model.component('comp1').physics('solid').feature('hmm1').set('E_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('E', 'ball_E');
model.component('comp1').physics('solid').feature('hmm1').set('nu_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('nu', 'ball_nu');
model.component('comp1').physics('solid').feature('hmm1').set('rho_mat', 'userdef');
model.component('comp1').physics('solid').feature('hmm1').set('rho', 'ball_density');
model.component('comp1').physics('solid').feature('hmm1').set('CalculateDissipatedEnergy', true);
model.component('comp1').physics('solid').feature('hmm1').label('Neo_hokean');
model.component('comp1').physics('tffs').prop('ShapeProperty').set('order_pressure', 4);
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
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis1').set('numelem', 'floor(ball_radius*vertical_ref/mesh_size)');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis2').label('radial_refinement');
model.component('comp1').mesh('mesh1').feature('fq1').feature('dis2').set('numelem', 'floor(3.14*ball_radius*(1-radial_ref)/(2*mesh_size))');
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
model.result.create('surf_plot1', 'PlotGroup2D');
model.result.create('surf_plot2', 'PlotGroup2D');
model.result.create('surf_plot3', 'PlotGroup2D');
model.result('tip_bw').create('lngr1', 'LineGraph');
model.result('tip_bw').feature('lngr1').set('xdata', 'expr');
model.result('tip_bw').feature('lngr1').selection.set([4]);
model.result('tip_bw').feature('lngr1').set('expr', 'z');
model.result('pressure_bw').create('lngr1', 'LineGraph');
model.result('pressure_bw').feature('lngr1').set('xdata', 'expr');
model.result('pressure_bw').feature('lngr1').selection.set([4]);
model.result('pressure_bw').feature('lngr1').set('expr', 'tffs.p');
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
model.result('surf_plot').create('surf4', 'Surface');
model.result('surf_plot').create('surf5', 'Surface');
model.result('surf_plot').create('surf6', 'Surface');
model.result('surf_plot').feature('surf2').set('expr', 'solid.sz');
model.result('surf_plot').feature('surf3').set('data', 'dset1');
model.result('surf_plot').feature('surf3').set('expr', '(solid.sGpr+solid.sGpz+solid.sGpphi)/3');
model.result('surf_plot').feature('surf4').set('data', 'dset1');
model.result('surf_plot').feature('surf4').set('expr', '-solid.Gequ*d(solid.curlUPHI,z)');
model.result('surf_plot').feature('surf5').set('data', 'dset1');
model.result('surf_plot').feature('surf5').set('expr', 'solid.Gequ*(1/r*d(r*solid.curlUPHI,r))');
model.result('surf_plot').feature('surf6').set('data', 'dset1');
model.result('surf_plot').feature('surf6').set('expr', 'solid.K*(1/r*d(r*u, r) + d(w, z))');
model.result('tip_profile_color').create('lngr1', 'LineGraph');
model.result('tip_profile_color').feature('lngr1').set('xdata', 'expr');
model.result('tip_profile_color').feature('lngr1').selection.set([4]);
model.result('tip_profile_color').feature('lngr1').set('expr', 'z');
model.result('pressure_profile_color').create('lngr1', 'LineGraph');
model.result('pressure_profile_color').feature('lngr1').set('xdata', 'expr');
model.result('pressure_profile_color').feature('lngr1').selection.set([4]);
model.result('pressure_profile_color').feature('lngr1').set('expr', 'tffs.p');
model.result('surf_plot1').create('surf3', 'Surface');
model.result('surf_plot1').feature('surf3').set('data', 'dset1');
model.result('surf_plot1').feature('surf3').set('expr', 'solid.J');
model.result('surf_plot2').create('surf3', 'Surface');
model.result('surf_plot2').feature('surf3').set('data', 'dset1');
model.result('surf_plot2').feature('surf3').set('expr', '');
model.result('surf_plot3').create('surf3', 'Surface');
model.result('surf_plot3').feature('surf3').set('data', 'dset1');
model.result('surf_plot3').feature('surf3').set('expr', '');
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
model.sol('sol1').feature('v1').set('clist', {'range(0,total_time/400,total_time)' '2.6666666666666667E-7[s]'});
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
model.sol('sol1').feature('t1').set('atolfactor', {'comp1_pfilm' '0.01' 'comp1_u' '0.01' 'comp1_solid_hmm1_pw' '0.1'});
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
model.sol('sol1').feature('t1').feature('st1').set('stopcondarr', {'1/timestep > 1e12' 'comp1.minop1(root.z) < 0 [m]'});
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
model.result('pressure_bw').set('titletype', 'manual');
model.result('pressure_bw').set('titlecolor', 'black');
model.result('pressure_bw').set('title', 'Pressure during impact');
model.result('pressure_bw').set('xlabel', 'r [mm]');
model.result('pressure_bw').set('ylabel', 'Physical pressure (Pa)');
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
model.result('residual_displacement').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('residual_displacement').feature('surf1').set('resolution', 'normal');
model.result('residual_pressure').label('Residual Pressure');
model.result('residual_pressure').set('innerinput', 'last');
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
model.result('Radius_time').set('solrepresentation', 'solutioninfo');
model.result('Radius_time').set('titletype', 'manual');
model.result('Radius_time').set('title', 'Radius of trapped air as a function of time');
model.result('Radius_time').set('xlabel', 'Time (s)');
model.result('Radius_time').set('ylabel', 'x (mm)');
model.result('Radius_time').set('xlabelactive', false);
model.result('Radius_time').set('ylabelactive', false);
model.result('Radius_time').feature('tblp1').set('xaxisdata', 1);
model.result('Radius_time').feature('tblp1').set('plotcolumninput', 'manual');
model.result('Radius_time').feature('tblp1').set('linewidth', 3);
model.result('surf_plot').label('Pressure');
model.result('surf_plot').set('symmetryaxis', true);
model.result('surf_plot').set('frametype', 'spatial');
model.result('surf_plot').feature('surf1').active(false);
model.result('surf_plot').feature('surf1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.result('surf_plot').feature('surf1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf1').set('resolution', 'normal');
model.result('surf_plot').feature('con1').active(false);
model.result('surf_plot').feature('con1').set('unit', [native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.result('surf_plot').feature('con1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('con1').set('levelmethod', 'levels');
model.result('surf_plot').feature('con1').set('levels', '10^{range(0,.1,2)}');
model.result('surf_plot').feature('con1').set('colortabletrans', 'nonlinear');
model.result('surf_plot').feature('con1').set('resolution', 'normal');
model.result('surf_plot').feature('surf2').active(false);
model.result('surf_plot').feature('surf2').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf2').set('resolution', 'normal');
model.result('surf_plot').feature('surf3').active(false);
model.result('surf_plot').feature('surf3').label('Pressure');
model.result('surf_plot').feature('surf3').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf3').set('resolution', 'normal');
model.result('surf_plot').feature('surf4').label('CurlCurlR');
model.result('surf_plot').feature('surf4').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf4').set('resolution', 'normal');
model.result('surf_plot').feature('surf5').active(false);
model.result('surf_plot').feature('surf5').label('CurlCurlZ');
model.result('surf_plot').feature('surf5').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf5').set('resolution', 'normal');
model.result('surf_plot').feature('surf6').active(false);
model.result('surf_plot').feature('surf6').label('DIvergence');
model.result('surf_plot').feature('surf6').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot').feature('surf6').set('resolution', 'normal');
model.result('tip_profile_color').label('Tip profile color');
model.result('tip_profile_color').set('innerinput', 'manualindices');
model.result('tip_profile_color').set('solnumindices', 'range(3,20,87)');
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
model.result('pressure_profile_color').set('innerinput', 'manualindices');
model.result('pressure_profile_color').set('solnumindices', 'range(3,20,143)');
model.result('pressure_profile_color').set('titletype', 'manual');
model.result('pressure_profile_color').set('titlecolor', 'black');
model.result('pressure_profile_color').set('title', 'Pressure during impact');
model.result('pressure_profile_color').set('xlabel', 'r [m]');
model.result('pressure_profile_color').set('ylabel', 'Physical pressure (Pa)');
model.result('pressure_profile_color').set('xlabelactive', false);
model.result('pressure_profile_color').set('ylabelactive', false);
model.result('pressure_profile_color').feature('lngr1').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pressure_profile_color').feature('lngr1').set('xdataexpr', 'if(r<1.5[mm],r,none)');
model.result('pressure_profile_color').feature('lngr1').set('xdataunit', '');
model.result('pressure_profile_color').feature('lngr1').set('xdatadescractive', true);
model.result('pressure_profile_color').feature('lngr1').set('xdatadescr', 'r [m]');
model.result('pressure_profile_color').feature('lngr1').set('linewidth', 3);
model.result('pressure_profile_color').feature('lngr1').set('legend', true);
model.result('pressure_profile_color').feature('lngr1').set('resolution', 'normal');
model.result('surf_plot1').label('Volume_change');
model.result('surf_plot1').set('symmetryaxis', true);
model.result('surf_plot1').set('frametype', 'spatial');
model.result('surf_plot1').feature('surf3').label('J');
model.result('surf_plot1').feature('surf3').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot1').feature('surf3').set('resolution', 'normal');
model.result('surf_plot2').label('Pressure_from_bulk');
model.result('surf_plot2').set('symmetryaxis', true);
model.result('surf_plot2').set('frametype', 'spatial');
model.result('surf_plot2').feature('surf3').label('J');
model.result('surf_plot2').feature('surf3').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot2').feature('surf3').set('resolution', 'normal');
model.result('surf_plot3').label('Shear_term');
model.result('surf_plot3').set('symmetryaxis', true);
model.result('surf_plot3').set('frametype', 'spatial');
model.result('surf_plot3').feature('surf3').label('J');
model.result('surf_plot3').feature('surf3').set('const', {'solid.refpntr' '0' 'Reference point for moment computation, r coordinate'; 'solid.refpntphi' '0' 'Reference point for moment computation, phi coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('surf_plot3').feature('surf3').set('resolution', 'normal');
model.result.export('anim1').label('Pressure_gif');
model.result.export('anim1').set('plotgroup', 'pressure_profile_color');
model.result.export('anim1').set('giffilename', '\Automated_data_nh\Soft_5_4242_mstip.gif');
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
model.result.export('anim2').set('plotgroup', 'tip_profile_color');
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
model.result.export('anim3').set('plotgroup', 'tip_profile19');
model.result.export('anim3').set('target', 'player');
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

out = model;