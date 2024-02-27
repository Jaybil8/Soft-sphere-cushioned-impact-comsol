%%
%launch server

Current_dir = pwd;
cd ('C:\Program Files\COMSOL\COMSOL60\Multiphysics\bin\win64')
system('comsolmphserver.exe &');
cd (Current_dir)

%%
%Establish connection
Current_dir = pwd;
cd ('C:\Program Files\COMSOL\COMSOL60\Multiphysics\mli');
mphstart(2036);
cd (Current_dir)