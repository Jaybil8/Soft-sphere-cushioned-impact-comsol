cd '\Users\jbilotto\Documents\EPFL E3 program 2021\research project\Experimental_data\'
data_soft = csvread("soft_impact.csv");
data_hard = csvread("rigid_impact.csv");

data_elastic_soft = data_soft(1:5, :);
b1 = log(data_elastic_soft(:,1))\log(data_elastic_soft(:,2));
%{
figure()
scatter(data_elastic_soft(:,1), data_elastic_soft(:,2))
hold on
plot(data_elastic_soft(:,1), b1*data_elastic_soft(:,1))
%}
data_elastic_hard = data_hard(1:5, :);
b2 = log(data_elastic_hard(:,1))\log(data_elastic_hard(:,2));

data_inertial_soft = data_soft(7:end, :);
b3 = log(data_inertial_soft(:,1))\log(data_inertial_soft(:,2))

data_inertial_hard = data_hard(7:end, :);
b4 = log(data_inertial_hard(:,1))\log(data_inertial_hard(:,2))