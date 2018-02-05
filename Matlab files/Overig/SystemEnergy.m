%% script to calculate energy balance of system
clear all
close all
clc
% Thanks Emiel for helping make this less messy :)
%% reading in output files from TuDat and extracting data
JsonInput = jsondecode(fileread('C:\Users\alks_\Documents\Matlab\Spaceflight_Assignment\Json_input\BenchmarkSettings.json'));

Graphs ='C:\Users\alks_\Documents\Matlab\Spaceflight_Assignment\Graphs\';
File = 'C:\Users\alks_\Documents\Matlab\Spaceflight_Assignment\Keppler_orbit_output\';
Dropbox = 'C:\Users\alks_\Dropbox\Benchmarking Tudat\Energy Graphs UnpKepler\';
NOC = 3;            % Number of cases
NOS = 3;            % Number of stepsizes
Data = dir(strcat(File,'*.dat'));
L = length (Data);
mu = 3.986004418e14;
count = 0;
%% Looping
while count <= NOC-1
  C0 = JsonInput.satellite.states(count+1);
  a = C0.semiMajorAxis;
  P = 2*pi*sqrt(a^3/mu);
for l = 1: L/NOC/NOS
for i = (count*(L/NOC))+l*NOS-2:(count*(L/NOC))+NOS*l
    s = strcat(File,Data(i).name);
    Import = importdata(s);
    
    %Importing parameters
Time = Import(:,1); 
rx = Import(:,2);
ry = Import(:,3);
rz = Import(:,4);
vx = Import(:,5);
vy = Import(:,6);
vz = Import(:,7);

ax = diff(vx)./diff(Time);
ay = diff(vy)./diff(Time);
az = diff(vz)./diff(Time);

Time = Time(1:end-1); 
rx = rx(1:end-1);
ry = ry(1:end-1);
rz = rz(1:end-1);
vx = vx(1:end-1);
vy = vy(1:end-1);
vz = vz(1:end-1);
hc = ones(size(Time))*-mu/(2*a);

% Calculations
rabs = sqrt(rx.^2+ry.^2+rz.^2);
h = 0.5*(vx.^2+vy.^2+vz.^2) - mu*(rabs.^-1);
hdot =(vx.*ax+vy.*ay+vz.*az)+...
      (mu./(rabs.^3)).*(rx.*vx+ry.*vy+rz.*vz);

% Plotting
name = Data(i).name(end-8:end-4);
semilogy (Time,abs(abs(h)-abs(hc)), 'DisplayName', name)
p = legend ('show');
set(p,'interpreter', 'none');
if i == (count*(L/NOC))+NOS*l       %Saving at last (depending on NOS) plot, closing and opening new fig
xlabel ('Time');
ylabel ('Specific energy');
set(gca,'yscale','log')
title(Data(i).name(1:end-6),'Interpreter', 'none' );
saveas(gcf,strcat(Graphs,Data(i).name(1:end-4),'.png'))
saveas(gcf,strcat(Dropbox,Data(i).name(1:end-4),'.png'))
 close
 figure('visible','off')
 hold on
end
end
end
count=count+1;
end
