%%Script for LabHack
clc;
clear all;
close all;

[filename, pathname] = uigetfile('*.csv','Select CSV file');
callname = strcat(pathname,'/',filename)
data = csvread(callname,1,0);

%data = csvread("/home/kevin/LabHack/S2_P1_R2_data.csv",1,0);
%data = csvread("/home/kevin/LabHack/filtered_output.csv",1,0);

% column 2 is vertical
% column 3 is horizontal

disp("done loading data")


polar_coord(:,2) = sqrt( data(:,2).*data(:,2) + data(:,2).*data(:,2) );
polar_coord(:,1) = atan2( data(2), data(2) );


%% Plotting script
maxx = max(abs(data(:,2)));
maxy = max(abs(data(:,2)));

supermax = 1.0 %max(maxx,maxy);
textloc = supermax * 0.9;

figure;
h = plot(data(1,2),data(1,2),"color",'g');
axis ([-supermax,supermax,-supermax,supermax], "square");
grid on;
grid minor;


timestamp = gmtime(time());
step = 5;

    h2 = text(-textloc,-textloc,strcat('Timestamp= ',num2str(data(1,1))));
    set(h2,'FontSize',16)

for i=10:step:length(data(:,1))
    yvector = [data(i,2),data(i-1,2),data(i-2,2),data(i-2,2),data(i-4,2),data(i-5,2),data(i-6,2),data(i-7,2),data(i-8,2),data(i-9,2)];
    xvector = [data(i,3),data(i-1,3),data(i-2,3),data(i-3,3),data(i-4,3),data(i-5,3),data(i-6,3),data(i-7,3),data(i-8,3),data(i-9,3)];
    
    set(h, 'YData', yvector,"marker",'o',"markersize", 15,"markerfacecolor",'g');
    set(h, 'XData', xvector,"marker",'o',"markersize", 15,"markerfacecolor",'g');
    
    set(h2,'String',strcat('Timestamp= ',num2str(data(i,1))));
    
    pause(data(i,1)-data(i-step,1));
    %pause(0.01);

end
