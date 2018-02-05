% Plot the diffrence vector over time,
% For a keplerian orbit with diffrent integrators and propagators.
clear all;
clc;

%set value to plot in same plot
%Which value to plot in same graps
 n = 0; % Stepsizes
% n = 1; % Propagators
% n = 2; % Integrators

% define directory locations.
%input
jsondirectory = 'BenchmarkSettings.json';
inputDirectory = 'SimulationOutput\TUDATBenchmarkskeplerapprox\';
%output
outputDirectory = 'MatlabOutput\';
TableOutput = 'TableOutput\';

%make directories if not present.
mkdir(outputDirectory);
mkdir(TableOutput);

% read json
settings = jsondecode(fileread(jsondirectory));

propagatorName = settings.propagators';
integratorsFixed = settings.integrators.fixed';
integratorsVariable = settings.integrators.variable';
integrators = [integratorsFixed,integratorsVariable];

fixedSteps = settings.integratorsettings.fixed.step_size_exp_begin:...
             settings.integratorsettings.fixed.step_size_exp_step:...
             settings.integratorsettings.fixed.step_size_exp_end;
fixedSteps = strcat('E',cellstr(num2str(fixedSteps')));

variableSteps = settings.integratorsettings.variable.rel_error_tol_exp_begin:...
                settings.integratorsettings.variable.rel_error_exp_step:...
                settings.integratorsettings.variable.rel_error_tol_exp_end;
variableSteps = strcat('E',cellstr(num2str(variableSteps')));


% Chance order of the for loops with variable n
loopmaster = {integrators;propagatorName;num2cell(1:3)};
if n>0 
    loopmaster = circshift(loopmaster,n,1);
end

close all        
for simulation = 1:size(settings.satellite.states,1)% loop over all cases
    for l1 = 1:size(loopmaster{1},2)            % loop 1
        for l2 = 1:size(loopmaster{2},2)        % loop 2
            figure('visible', 'off') % make invisible figures
            for l3 = 1:size(loopmaster{3},2)    % loop 3
                
                property1 = loopmaster{1}{l1};
                property2 = loopmaster{2}{l2};
                property3 = loopmaster{3}{l3};
                
                % find integrator and propagator names.
                integrator = eval(strcat('property',num2str(mod(n+0,3)+1)));
                propagator = eval(strcat('property',num2str(mod(n+1,3)+1)));
                varname =         strcat('property',num2str(mod(n+2,3)+1));
                
                % check if current integrator is a fixed step integrator
                if  any(strcmp(integrator,integratorsFixed))
                    step = fixedSteps{eval(varname)};
                else
                    step = variableSteps{eval(varname)};
                end
                assignin('base',varname,step)
                
                % find the number of the current integrator and propagator
                intnum   = eval(strcat('l',num2str(mod(n+0,3)+1)));
                propnum  = eval(strcat('l',num2str(mod(n+1,3)+1)));
                stepnum  = eval(strcat('l',num2str(mod(n+2,3)+1)));
                
                
                % Import results from TUDAT
                tit = strcat('Case_',num2str(simulation-1),'_',...
                              propagator,'_',integrator);
                file = strcat(tit,'_',step,'.dat');
                path = strcat(inputDirectory,file);
                results = importdata(path);
                t = results(:,1);

                % Calculate Analytical Solution
                R = AnalyticSolution(t,settings,simulation);

                % Calculate diffrence vector
                O = results(:,2:4)-R(1:size(results(:,2:4),1),:);
                error = sqrt(sum(O.^2,2));
                
                rsslim = find(t>=0.9*t(end),1);
                rss = sum(error(rsslim:end).^2)/(size(error,1)-rsslim+1);
                rsstable(stepnum,propnum,intnum) = rss;
                
                leg = strcat(property3,', RSS = ',num2str(rss,3));
                
                % dont show 1st entry, error(1)=0 does not work on log plot
                semilogy(t(2:end),error(2:end),'DisplayName',leg)
                hold on
                legend('show')
            end
            tit = {strcat('Kepler Orbit');...
                   strcat('Case_',num2str(simulation-1),...
                   '_',property1,'_',property2)};
            
            title(tit,'Interpreter', 'none'); %disable latex interpreter
            xlabel('Time [seconds since J2000]');
            ylabel('Error in position [m]');
        end
    end
    % Make latex tables
    for i = 1:length(integrators)
        tablegen(TableOutput,rsstable(:,:,i),settings,simulation,integrators{i})
    end
end


% save all figures
figHandles = get(0,'Children');
for f = figHandles'
    % The title of the figure will be the title of the image.
    tit = get(get(get(f,'CurrentAxes'),'title'),'String');
    outputDir = strcat(pwd,'\',outputDirectory,tit);
    saveas(f,sprintf('%s.png',outputDir{2}))
end
close all %close invisible figures
%% Functions
function [ E ] = KeplerEq(manom,ecc)
%solves inverse kepler equation with a newtons-raphson approximation
    tol=1e-12;
    E=zeros(length(manom),1);
    for i=1:length(manom)
        Mi = manom(i);
        diff=inf;
        E(i,:)=Mi;
        while diff>tol
            f = E(i,:) - ecc*sin(E(i,:));
            df = 1-ecc*cos(E(i,:));
            E(i,:) = E(i,:)-((f-Mi)/(df));
            diff=abs(Mi-f);
        end
    end
end
function [ R ] = Orbital2State2(a,e,i,w,RA,v)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% This function computes the position vector r, from the
% classical orbital elements.
%
% a = Semi major axis [m]
% e = eccentricity [-]
% i = inclination of the orbit [rad] 
% w = argument of perigee [rad]
% RA = right ascension of the ascending node [rad]
% v = true anomaly [rad]
% r - Position vector magnitude [m]
% ------------------------------------------------------------
r=(a*(1-e^2))./(1+e*cos(v));

l1 = cos(w)*cos(RA)-sin(w)*sin(RA)*cos(i);
m1 = cos(w)*sin(RA)+sin(w)*cos(RA)*cos(i);
n1 = sin(w)*sin(i);
l2 = -sin(w)*cos(RA)-cos(w)*sin(RA)*cos(i);
m2 = -sin(w)*sin(RA)+cos(w)*cos(RA)*cos(i);
n2 = cos(w)*sin(i);
X = r.*(l1*cos(v)+l2*sin(v));
Y = r.*(m1*cos(v)+m2*sin(v));
Z = r.*(n1*cos(v)+n2*sin(v));

R=[X,Y,Z];
end
function [ R ] = AnalyticSolution(t,json,sim)
% This function calculates the analytical solution of an orbit.
%
% t = time vector [s]
% json = the json input file containing information about the simulation
% sim = the case number.
%
%
% This function requires two other functions:
%
% KeplerEq - to reverse the kepler equation
% Orbital2State2 - to calculate position from keplerian elements.
    

    set = json.satellite.states(sim);
    
    % Earth’s gravitational parameter [m^3/s^2]
    mu = json.gravparam.Earth;

    %  Inital state
    a     =  set.semiMajorAxis;           % [m] Semimajor axis
    e     =  set.eccentricity;            % [-] Eccentricity
    i     =  set.inclination;             % [rad] Inclination
    RAAN  =  set.longitudeOfAscendingNode;% [rad] longitude of ascending node
    omega =  set.argumentOfPeriapsis;     % [rad] Argument of perigee
    theta =  set.trueAnomaly;             % [rad] True anomaly at epoch

    T = 2*pi*sqrt(a^3/mu);                % Period [s]
    
    % Vectorized position calclation
    E0 = 2*atan(tan(theta/2)*sqrt((1-e)/(1+e)));
    M0 = E0  - e*sin(E0);
    t0 = (M0/2*pi)*T;
    tp = t0 + t;
    M = 2*pi*(tp/T);
    E = KeplerEq(M,e);
    theta = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
    R = Orbital2State2(a,e,i,omega,RAAN,theta);
end
function [   ] = tablegen(outputdir,array,json,sim,intname)
% This function generates a latex table
%
% For this function the latexTable add-on from Eli Duenisch has to be
% installed.
% Version 1.21 used by creation of this script.
%
% outputdir = [string], gives the path to the output directory
% array = [2 dimensional matrix], contains values to put in table
% json = [struct], json input file containing information about the simulation
% sim = [number], the case number
% intname = [string], name of the current integrator

input.data = array;

propagatorName = json.propagators';
integratorsFixed = json.integrators.fixed';
fixedSteps = json.integratorsettings.fixed.step_size_exp_begin:...
             json.integratorsettings.fixed.step_size_exp_step:...
             json.integratorsettings.fixed.step_size_exp_end;
fixedSteps = strcat('1e',cellstr(num2str(fixedSteps')));

variableSteps = json.integratorsettings.variable.rel_error_tol_exp_begin:...
                json.integratorsettings.variable.rel_error_exp_step:...
                json.integratorsettings.variable.rel_error_tol_exp_end;
variableSteps = strcat('1e',cellstr(num2str(variableSteps')));

% Collumn labels
input.tableColLabels = propagatorName;

% Row labels
if any(strcmp(intname,integratorsFixed))
    input.tableRowLabels = fixedSteps;
else
    input.tableRowLabels = variableSteps;
end

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'.

input.dataFormat = {'%.3e',4}; % three digits precision

% LaTex table caption:
input.tableCaption = ['RSS value for Case ',num2str(sim-1),...
                      ' and a ',intname,' integrator'];

% LaTex table label:
input.tableLabel = ['tab:keplerapprox_rss_C',num2str(sim-1),'_',intname];

% call latexTable:
latex = latexTable(input);

% save LaTex code as file
filename = strcat(outputdir,'Keplerapprox_C',num2str(sim-1),...
                  '_int_',intname,'_RSS_Table.tex');
fid=fopen(filename,'w');
[nrows,~] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
end