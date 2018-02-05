% Input the keplerian elements AnalyticSolution
clear all;
close all
clc;
%format shortE

settings = jsondecode(fileread('BenchmarkSettings.json'));
inputDirectory = 'SimulationOutput\TUDATBenchmarkskeplerapprox\';
outputDirectory = 'MatlabOutput\';
TableOutput = 'TableOutput\';


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
            
loopmaster = {integrators;propagatorName;num2cell(1:3)};
%set value to plot in same plot
%Which value to plot in same graps
 n = 0; % Stepsizes
% n = 1; % Propagators
% n = 2; % Integrators

if n>0 
    loopmaster = circshift(loopmaster,n,1);
end
            
for simulation = 1:size(settings.satellite.states,1)% loop over all cases
    for l1 = 1:size(loopmaster{1},2)            % loop 1
        for l2 = 1:size(loopmaster{2},2)        % loop 2
            figure('visible', 'off') % make invisible figures
            for l3 = 1:size(loopmaster{3},2)    % loop 3
                
                property1 = loopmaster{1}{l1};
                property2 = loopmaster{2}{l2};
                property3 = loopmaster{3}{l3};
                
                integrator = eval(strcat('property',num2str(mod(n+0,3)+1)));
                propagator = eval(strcat('property',num2str(mod(n+1,3)+1)));
                varname =         strcat('property',num2str(mod(n+2,3)+1));
                
                intnum   = eval(strcat('l',num2str(mod(n+0,3)+1)));
                propnum  = eval(strcat('l',num2str(mod(n+1,3)+1)));
                stepnum  = eval(strcat('l',num2str(mod(n+2,3)+1)));
                
                % check if current integrator is a fixed step integrator
                if  any(strcmp(integrator,integratorsFixed))
                    step = fixedSteps{eval(varname)};
                else
                    step = variableSteps{eval(varname)};
                end
                assignin('base',varname,step)
            
                % Import results from TUDAT
                tit = strcat('Case_',num2str(simulation-1),'_',...
                              propagator,'_',integrator);
                file = strcat(tit,'_',step,'.dat');
                path = strcat(inputDirectory,file);
                results = importdata(path);
                t = results(:,1);

                % Calculate Analytical Solution
                R = AnalyticSolution(t,settings,simulation);

                % Error function
                O = results(:,2:4)-R(1:size(results(:,2:4),1),:);
                error = sqrt(sum(O.^2,2));
                rsslim = find(t>=0.9*t(end),1);
                rss = sum(error(rsslim:end).^2)/(size(error,1)-rsslim+1);
                rsstable(stepnum,propnum,intnum) = rss;
                
                leg = strcat(property3,', RSS = ',num2str(rss,3));
                semilogy(t(2:end),error(2:end),'DisplayName',leg) % do not show 1st entry, error(1)=0
                hold on
                legend('show')
            end
            tit = {strcat('Kepler Orbit');strcat('Case_',num2str(simulation-1),'_',property1,'_',property2)};
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
    outputDir = strcat(pwd,'\',outputDirectory,get(get(get(f,'CurrentAxes'),'title'),'String'));
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
% This function computes the state vector (r,v) from the
% classical orbital elements.
%
% a = Semi major axis [m]
% e = eccentricity [-]
% i = inclination of the orbit [rad] 
% w = argument of perigee [rad]
% RA = right ascension of the ascending node [rad]
% v = true anomaly [rad]
% r - position vector in the geocentric equatorial frame [m]
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
    
    
    set = json.satellite.states(sim);
    % constants
    mu = 3.986004418e14;                 % Earth’s gravitational parameter [m^3/s^2]

    %  Inital state
    a     =  set.semiMajorAxis;           % [m] Semimajor axis
    e     =  set.eccentricity;            % [-] Eccentricity
    i     =  set.inclination;             % [rad] True anomaly at epoch
    omega =  set.argumentOfPeriapsis;     % [rad] Right ascension of the ascending node
    RAAN  =  set.longitudeOfAscendingNode;% [rad] Inclination
    theta =  set.trueAnomaly;             % [rad] Argument of perigee

    T = 2*pi*sqrt(a^3/mu);          % Period
    
    E0 = 2*atan(tan(theta/2)*sqrt((1-e)/(1+e)));
    M0 = E0  - e*sin(E0);
    t0 = (M0/2*pi)*T;
    tp = t0 + t;
    M = 2*pi*(tp/T);
    E = KeplerEq(M,e);
    theta = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
    R = Orbital2State2(a,e,i,omega,RAAN,theta);
end
function [   ] = tablegen(ourputdir,array,json,sim,intname)
% numeric values you want to tabulate:
% this field has to be an array or a MATLAB table
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

input.tableColLabels = propagatorName;
if any(strcmp(intname,integratorsFixed))
    input.tableRowLabels = fixedSteps;
else
    input.tableRowLabels = variableSteps;
end

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'.

input.dataFormat = {'%.3e',4}; % three digits precision

% LaTex table caption:
input.tableCaption = ['RSS value for Case ',num2str(sim-1),' and a ',intname,' integrator'];

% LaTex table label:
input.tableLabel = ['tab:keplerapprox_rss_C',num2str(sim-1),'_',intname];

% call latexTable:
latex = latexTable(input);

% save LaTex code as file
filename = strcat(ourputdir,'Keplerapprox_C',num2str(sim-1),'_int_',intname,'_RSS_Table.tex');
fid=fopen(filename,'w');
[nrows,~] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
end