clear all
% Import data
% store data in big 3D array
files = dir('ISS_Output\*.dat');
json = jsondecode(fileread('ISSD\MultiCase_ISSD.json'));
outputDirectory = 'MatlabOutput\';
TableOutput = 'TableOutput\';


integrators = [json.integrators.fixed; json.integrators.variable];
propagators = json.propagators;

fixedSteps = json.integratorsettings.fixed.step_size_exp_begin:...
             json.integratorsettings.fixed.step_size_exp_step:...
             json.integratorsettings.fixed.step_size_exp_end;
fixedSteps = strcat('E',cellstr(num2str(fixedSteps')));

variableSteps = json.integratorsettings.variable.rel_error_tol_exp_begin:...
                json.integratorsettings.variable.rel_error_exp_step:...
                json.integratorsettings.variable.rel_error_tol_exp_end;
variableSteps = strcat('E',cellstr(num2str(variableSteps')));

loopmaster = {integrators';propagators';num2cell(1:3)};
%set value to plot in same plot
%Which value to plot in same graps
% n = 0; % Stepsizes
 n = 1; % Propagators
% n = 2; % Integrators

if n>0 
    loopmaster = circshift(loopmaster,n,1);
end


for l1 = 1:size(loopmaster{1},2) %loop 1
    for l2 = 1:size(loopmaster{2},2)% loop 2
        figure('visible', 'off') % make invisible figures
        for l3 = 1:size(loopmaster{3},2)% loop 3
            property1 = loopmaster{1}{l1};
            property2 = loopmaster{2}{l2};
            property3 = loopmaster{3}{l3};
            
            intnum   = eval(strcat('l',num2str(mod(n+0,3)+1)));
            propnum  = eval(strcat('l',num2str(mod(n+1,3)+1)));
            stepnum  = eval(strcat('l',num2str(mod(n+2,3)+1)));
            
            integrator = eval(strcat('property',num2str(mod(n+0,3)+1)));
            propagator = eval(strcat('property',num2str(mod(n+1,3)+1)));
            varname =         strcat('property',num2str(mod(n+2,3)+1));
            % check if current integrator is a fixed step integrator
            if  any(strcmp(integrator,json.integrators.fixed))
                stepsize = fixedSteps{eval(varname)};
            else
                stepsize = variableSteps{eval(varname)};
            end
            assignin('base',varname,stepsize)
            
            % run actual calculation  
            ii=1;
            clear data
            for i = 1:size(files,1)
                if ~isempty(strfind(files(i).name,property1)) && ...
                   ~isempty(strfind(files(i).name,property2)) && ...
                   ~isempty(strfind(files(i).name,property3))
                    % download data from files
                    path = strcat(files(i).folder,'\',files(i).name);
    
                    data(:,:,ii) = importdata(path);

                    % get mu
                    index = strfind(files(i).name,'_');
                    name = files(i).name(1:index(1)-1);

                    mu(ii) = json.gravparam.(name);

                    ii = ii + 1;
                end
            end
            
            e = zeros(size(data,1),size(data,3));
            for i=1:size(data,3)
                k = 1/2*mu(i)*sum(data(:,5:7,i).^2,2);
                pot=-mu(i)*json.gravparam.('Sun')./sqrt(sum(data(:,2:4,i).^2,2));
                for j=1:size(data,3)
                    if i==j
                        % no gravitational influnce from itsself
                    else
                        p = -1/2*mu(i)*mu(j)./sqrt(sum((data(:,2:4,i)-data(:,2:4,j)).^2,2));
                        pot = pot+p;
                    end
                end
                e(:,i) = k+pot;
            end

            eTotal = sum(e,2);
            t = data(:,1,1);
            rsslim = find(t>=0.9*t(end),1);
            rss = sum(eTotal(rsslim:end).^2)/(size(eTotal,1)-rsslim+1);
            rsstable(stepnum,propnum,intnum) = rss;
            
            leg = strcat(property3,' RSS = ',num2str(rss,3));
            plot(data(:,1,1),eTotal/eTotal(1),'DisplayName',leg)
            hold on
            legend('show')
        end
        tit = {strcat('Total energy over time');strcat(property1,'_',property2)};
        title(tit,'Interpreter', 'none')
        xlabel('time [seconds since J2000]')
        ylabel('E(t)/E(0)')
    end
end
% make tables
for i = 1:length(integrators)
    tablegen(TableOutput,rsstable(:,:,i),json,integrators{i})
end

% save all figures
figHandles = get(0,'Children');
for f = figHandles'
    outputDir = strcat(pwd,'\',outputDirectory,get(get(get(f,'CurrentAxes'),'title'),'String'));
    saveas(f,sprintf('%s.png',outputDir{2}))
end
close all %close invisible figures

function [   ] = tablegen(ourputdir,array,json,intname)
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
input.tableCaption = ['RSS values for the ',intname,' integrator'];

% LaTex table label:
input.tableLabel = ['tab:ISSD_rss_',intname];

% call latexTable:
latex = latexTable(input);

% save LaTex code as file
filename = strcat(ourputdir,'ISSD_int_',intname,'_RSS_Table.tex');
fid=fopen(filename,'w');
[nrows,~] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
end