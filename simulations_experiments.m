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

for i=1:length(vel_list)
    % soft impactor
    simulation_comsol(vel_list(i), 20e-6, 7.1e-3, 1140, .47, 1.81e-5, 1200,12e-6,250e3,30)
    % hard impactor
    simulation_comsol(vel_list(i), 20e-6, 7.1e-3, 1140, .47, 1.81e-5, 1200,12e-6,1100e3,30)
end

for i=1:length(vel_list)
    % soft impactor
    simulation_comsol(vel_list(i), 10e-6, 7.1e-3, 1140, .47, 1.81e-5, 1200,12e-6,250e3,60)
    % hard impactor
    simulation_comsol(vel_list(i), 10e-6, 7.1e-3, 1140, .47, 1.81e-5, 1200,12e-6,1100e3,60)
end


