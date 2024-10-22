%% GRASPING OBJECTS WITH A ACOUSTIC SENSORY SUBSTITUTION GLOVE

% Program developed by Carlos de Paz (carlos.paz@uam.es,
% carlosdepazrios94@gmail.com).

clc; clear all; close all;
pp=[];
while isempty(pp)  %%This while is to create a folder for each participant
    pp = input('Please, introduce participant number: \n');
    fprintf('Participant number = %i \n', pp)
    reply = input('Is this correct? Y/N [Y]: \n','s');
    if isempty(reply)|| reply=='Y' || reply =='y'
    else
        pp=[];
    end
    if pp < 10
        Pnumber=['PP00',num2str(pp)];
    elseif pp>=10 && pp<100
        Pnumber=['PP0',num2str(pp)];
    else
        Pnumber=['PP',num2str(pp)];
    end
%     cd('xxxxx')
    DirID=Pnumber;
    mkdir(DirID);
    if isempty(mess)
    else
        fprintf(mess)
        pp=[];
    end
end
cd(DirID); MainCD = cd;
fprintf('Press a key to start the Training \n')
pause()
RepeatTrials = []; %Matrix for the recovery trials

% --------------------------Experiment----------------------------------
clc; 

mkdir('Experiment');cd('Experiment')

training = 0;
    
 for trial=1:81 % 27 conditions per 3 repetitions
     load InitialConditions.mat 

    % Display Information
    fprintf('Trial = %i \n', trial);
    %Call Qualysis and get the conf for each trial
    [conf, condition] = GetIniHvA(pp,trial,training);% We set a priori the experimental conditions regarding the participant, trial and the training
    conf.q = QMC;  %Call Qualysis and get the conf for each trial
    
    
    %Calculate variables to draw the circle (object)
    alfa = 0:0.01:2*pi;
    xt = conf.Sizet.*cos(alfa)+conf.Xt; 
    yt = conf.Sizet.*sin(alfa)+conf.Yt;
    
    if conf.Yt == 30 || conf.Yt == 25.9807621135332 || conf.Yt == 15
        distance_print = 30;
    elseif conf.Yt == 25 || conf.Yt == 21.6506350946110 || conf.Yt == 12.5
        distance_print = 25;
    elseif conf.Yt == 20 || conf.Yt == 17.3205080756888 || conf.Yt == 10
        distance_print = 20;
    end
    fprintf('Angle = %i \n', conf.Angle);
    fprintf('Distance = %i \n', distance_print);
    fprintf('Size = %i \n', conf.Sizet);
    
    %Pure tone (sound) Parameters 
    fs = 20500;  % sampling frequency
    duration = 0.1;  % .025; %.15;
    values=0:1/fs:duration;
    s = zeros(2,2051); %513 %1026   
    a = .009;  b = - .5;   c = 10;
    
    %Define the empty vectors
    xvectori = NaN(1,100); yvectori = xvectori;
    xvectort = xvectori;   yvectort = xvectori;
    
    %Define the empty scalars
    xintercept_i = NaN;          yintercept_i = xintercept_i;
    xintercept_t = xintercept_i; yintercept_t = xintercept_i;
    
    %Define the object initial position (X,Y,Z)

    [data, obj_ini] = QMC(conf.q,1); % get data from Qualisys
    
    % Check if the object is well detected (stop the script until the
    % object is detected)    
    Error = find(obj_ini(:,1)==0);
    flag_print=0;
      while Error ~= 0
        [data, obj_ini] = QMC(conf.q,1);
        if flag_print == 0
           fprintf('Looking for the object posini \n')
           flag_print = 1;
        end
           Error = find(obj_ini(:,1)==0);
      end
     clear Error
     obj_ini = obj_ini./10;
    
    % Dynamical subplot (conditions vs actual position of the object) 
    obj_pos = zeros(3,1);

    figure('units', 'normalized', 'outerposition', [0 0 1 1])
    subplot(1,2,1)
    hold on
    axis([10 60 0 32.5])
    for i = 1:4
        plot(linspace(50,M30(1,i),100), linspace(0,M30(2,i),100),'k--','LineWidth', 6)
        plot(M20(1,i), M20(2,i),'ko', 'LineWidth', 6)
        plot(M30(1,i), M30(2,i),'ko', 'LineWidth', 6)
    end
    for i = 1:3
        plot(M22(1,i), M22(2,i),'kx', 'LineWidth', 6)
        plot(M28(1,i), M28(2,i),'kx', 'LineWidth', 6) 
    end
    plot(conf.Xt,conf.Yt, 'rx', 'LineWidth', 6)
    title(['Angle = ', num2str(conf.Angle), ', Dist = ' num2str(distance_printer), 'Size = ',num2str(conf.Sizet)], 'FontSize', 32)
    hold off

    subplot(1,2,2)
    axis([10 60 0 40])
    hold on
    plot(xt,yt, 'k', 'LineWidth', 6)
    plot(conf.Xt, conf.Yt,'kx', 'LineWidth', 6)
    Pos_Y = plot([obj_pos(1,1), obj_pos(1,1)],[obj_pos(2,1)-conf.Sizet, obj_pos(2,1)+conf.Sizet],'r--', 'LineWidth', 6);
    Pos_X = plot([obj_pos(1,1)- conf.Sizet, obj_pos(1,1)+ conf.Sizet],[obj_pos(2,1), obj_pos(2,1)],'r--', 'LineWidth', 6);
    Centre = plot(obj_pos(1,1),obj_pos(2,1),'rx', 'Linewidth', 6); 
    
    while (abs(obj_pos(1,1) - conf.Xt) > .5 || abs(obj_pos(2,1) - conf.Yt) > .5)
    
        [data,obj_pos] = QMC(conf.q,1); 
        obj_pos = obj_pos./10;
        set(Pos_Y, 'Xdata', [obj_pos(1,1), obj_pos(1,1)], 'Ydata', [obj_pos(2,1)-conf.Sizet, obj_pos(2,1)+conf.Sizet]);
        set(Pos_X, 'Xdata', [obj_pos(1,1) - conf.Sizet, obj_pos(1,1)+ conf.Sizet], 'Ydata', [obj_pos(2,1), obj_pos(2,1)]);
        set(Centre, 'Xdata', obj_pos(1,1), 'Ydata', obj_pos(2,1));
        drawnow

    end
    hold off
    pause(1.5)
    close all
    [data,obj_ini] = QMC(conf.q,1); 
    obj_ini = obj_ini./10; data = data./10; %actualize the obj_ini pos
    conf.Xobj = obj_ini(1,1);
    conf.Yobj = obj_ini(2,1);

    
    t = 0; z = 0;%Just a varaible to control (On/Off) the while loop
    MatrixData = NaN(1,1); n = 0; %number of rows 
    fprintf('Press a key to Start the Trial \n')
    pause()
    
    tic %Calculate the time of each trial
    
    figure ('units', 'normalized', 'position', [0.05 0.35 0.35 .5])
    title(['Trial = ', num2str(trial)], 'FontSize', 14)
    hold on
    axis([10 60 0 40])

    H1 = uicontrol('Position', [5 10 50 20],'Style', 'PushButton', 'BackgroundColor','g', 'FontWeight', 'bold','String', 'NEXT', 'Callback', 't=NaN;'); %Next button, if it's pushed the next trial will begin
    H2 = uicontrol('Position', [5 60 50 20],'Style', 'PushButton', 'BackgroundColor','r', 'FontWeight', 'bold','String', 'STOP', 'Callback', 'z=1;'); %Stop button, if it's pushed the experiment will end
    H3 = uicontrol('Position', [5 35 50 20],'Style', 'PushButton', 'BackgroundColor',[.4 .4 1], 'FontWeight', 'bold','String', 'REPEAT', 'Callback', 't=NaN; z=2;'); %Repeat button, if it's pushed this trial will be repeat  
    plot(conf.Xt,conf.Yt,'kx',xt,yt,'k-') %Center of the circle and the perimeter
    
    %Define the variables (empty) that we will plot
    hVI = plot(xvectori,yvectori,'r--');
    hVT = plot(xvectort,yvectort,'b--');
    hIntI = plot(xintercept_i,yintercept_i, 'ro');
    hIntT = plot(xintercept_t,yintercept_t, 'bo');
        
    while ~isnan(t)
        
        [data,obj_pos] = QMC(conf.q,1); 

        % "While" solution for the critical error (data = 0) 
        Error = find(data(1:2,1:4)==0);
        flag_print=0;
        while Error ~= 0
            [data,obj_pos] = QMC(conf.q,1);
            if flag_print == 0
               fprintf('Looking for data \n')
               flag_print = 1;
            end
            Error = find(data(1:2,1:4)==0);
        end
        clear Error
        data = data./10; obj_pos = obj_pos./10; %Transform the data to cm
        
       
        %Calculate whether the participant is pointing to objects with each finger. In that
        %case, calculate the coordinates of the intercept and the distance
        [data,m_i,b_i,xintercept_i,yintercept_i,distance_i,m_t,b_t,xintercept_t,yintercept_t,distance_t] = intercept(obj_ini(1,1),obj_ini(2,1),conf.Sizet,data);
         
        % Index finger
        if isnan(distance_i) %If there is not an intercept, just plot a line from xdistal_i to xfinali (X position for conf.Yt) 
            xfinali = (max(yt)-b_i)/m_i;
            xvectori = linspace(data(1,1),xfinali,100);
            yvectori = linspace(data(2,1),max(yt),100);
            A_i = 0;
            
        else
            xvectori = linspace(data(1,1),xintercept_i,100);
            yvectori = linspace(data(2,1),yintercept_i,100);
            A_i = a.*distance_i.^2 + b.*distance_i + c;

        end
        
        s(1,:) = A_i.*sin(2.*pi.* 200.*values); %pure tone for the left channel (~ index finger)

        % Thumb finger
        if isnan(distance_t)
            xfinalt = (max(yt)-b_t)/m_t;
            xvectort = linspace(data(1,3),xfinalt,100);
            yvectort = linspace(data(2,3),max(yt),100);
            A_t = 0;
            
        else
            xvectort = linspace(data(1,3),xintercept_t,100);
            yvectort = linspace(data(2,3),yintercept_t,100);
            A_t = a.*distance_t.^2 + b.*distance_t + c;

        end
        
        s(2,:) = A_t.*sin(2.*pi.* 200.*values); %pure tone for the rigth channel (~ thumb)

            
        %After calculate the value if each data, the plot is going to update 
        set(hVI, 'Xdata', xvectori, 'Ydata', yvectori);
        set(hVT, 'Xdata', xvectort, 'Ydata', yvectort);
        set(hIntI, 'Xdata', xintercept_i, 'Ydata', yintercept_i);
        set(hIntT, 'Xdata', xintercept_t, 'Ydata', yintercept_t);
        drawnow
        hold off
        
        sound(s,fs)
        
        %Saving the data in a matrix 
        data_save = data(1:2,1:5); data_save=data_save(:)'; %4 finger markers + wrist
        n = n+1; %add a row for the MatrixData
        time = toc; %Calculate the time of each trial
        
        MatrixData(n,1:21) = [data_save, m_i, b_i, distance_i, xintercept_i,yintercept_i, m_t, b_t, distance_t, xintercept_t,yintercept_t, time]; %21
        
        if z==1 %If Stop button is pressed, stop the experiment
            fprintf('Last Trial = %i \n', trial) 
            return
            close all
        end
        
        if z==2 %If Repeat button is pressed, it shows which trial needs to be repeated and start the next one
           RepeatTrials (trial,1:3) = [trial, conf.Xt, conf.Yt];
           RepeatTrials=nonzeros(RepeatTrials(:,1));
           fprintf('Repeat Trial = %i \n', trial)
           z=0;
        end
        
%         if abs(data(1,6) - conf.Xt) > 1 || abs(data(2,6) - conf.Yt) > 1 || round(data (3,6)) ~= 14 %If the participant move the object, the trial will finish (Grasped = 1; Hit = 0)
        if sum(abs(obj_pos - obj_ini))>.5           
        daq.reset()
           close all
           conf.Answer = input ('Has the participant grasped the target? (N = 0; Y = 1; Repeat = 3) \n');
           fprintf('Next trial')
           pause()
               if conf.Answer == 3
                 RepeatTrials (trial,1:3) = [trial, conf.Xt, conf.Yt];
                 RepeatTrials=nonzeros(RepeatTrials(:,1));
                 fprintf('Repeat Trial = %i \n', trial)
               end
           t=NaN; %this value will stop the loop and start the next trial
        end
        
    end
    
    MatrixDataTable = array2table(MatrixData, 'VariableNames', {'Xdistal_i','Ydistal_i', 'Xprox_i', 'Yprox_i','Xdista_t', 'Ydistal_t', 'Xprox_t','Yprox_t','Xwrist','Ywrist','m_i', 'b_i', 'distance_i', 'xintercept_i','yintercept_i', 'm_t', 'b_t', 'distance_t', 'xintercept_t','yintercept_t', 'time'});
    close all

% %     ClosePort(PS)
    QMC(conf.q, 'quit')
%     clear PS
    
    if trial<10
        Tnumber = ['TR00', num2str(trial)];
    else
        Tnumber = ['TR0', num2str(trial)];
    end
    FileID = [Pnumber,Tnumber];
    save(FileID, 'conf', 'MatrixData', 'MatrixDataTable');
    fprintf('End of the trial \n')
    clearvars -except training alpha trial pp Pnumber Tnumber a b c RepeatTrials MainCD
 end
cd(MainCD)
fprintf('-----------------End of the Experiment------------------------ \n')

%% 1. GetIniHvA: this functions set the coordinates of the to-be-grasped object

% Orientation (90-120-150) --> Distance (30-25-20) --> Size (4-3-2) = 27 x 3 reps

function [ conf, condition ] = GetIniHvA(pp,trial,training)
    load('InitialConditions_HvA.mat') 
    if training == 0
        load HvAOrder.mat
        condition = HvAOrder(pp,trial);
        switch condition
            case 1
                conf.Xt = M30(1,1);
                conf.Yt = M30(2,1);
                conf.Sizet = 4;
                conf.Angle = 90;
            case 2
                conf.Xt = M30(1,1);
                conf.Yt = M30(2,1);
                conf.Sizet = 3;
                conf.Angle = 90;
            case 3 
                conf.Xt = M30(1,1);
                conf.Yt = M30(2,1);
                conf.Sizet = 2;
                conf.Angle = 90;
            case 4
                conf.Xt = M25(1,1);
                conf.Yt = M25(2,1);
                conf.Sizet = 4;
                conf.Angle = 90;
            case 5
                conf.Xt = M25(1,1);
                conf.Yt = M25(2,1);
                conf.Sizet = 3;
                conf.Angle = 90;
            case 6
                conf.Xt = M25(1,1);
                conf.Yt = M25(2,1);
                conf.Sizet = 2;
                conf.Angle = 90;
            case 7 
                conf.Xt = M20(1,1);
                conf.Yt = M20(2,1);
                conf.Sizet = 4;
                conf.Angle = 90;
            case 8
                conf.Xt = M20(1,1);
                conf.Yt = M20(2,1);
                conf.Sizet = 3;
                conf.Angle = 90;   
            case 9
                conf.Xt = M20(1,1);
                conf.Yt = M20(2,1);
                conf.Sizet = 2;
                conf.Angle = 90; 
            case 10
                conf.Xt = M30(1,2);
                conf.Yt = M30(2,2);
                conf.Sizet = 4;
                conf.Angle = 120;
            case 11
                conf.Xt = M30(1,2);
                conf.Yt = M30(2,2);
                conf.Sizet = 3;
                conf.Angle = 120;
            case 12 
                conf.Xt = M30(1,2);
                conf.Yt = M30(2,2);
                conf.Sizet = 2;
                conf.Angle = 120;
            case 13
                conf.Xt = M25(1,2);
                conf.Yt = M25(2,2);
                conf.Sizet = 4;
                conf.Angle = 120;
            case 14
                conf.Xt = M25(1,2);
                conf.Yt = M25(2,2);
                conf.Sizet = 3;
                conf.Angle = 120;
            case 15
                conf.Xt = M25(1,2);
                conf.Yt = M25(2,2);
                conf.Sizet = 2;
                conf.Angle = 120;
            case 16 
                conf.Xt = M20(1,2);
                conf.Yt = M20(2,2);
                conf.Sizet = 4;
                conf.Angle = 120;
            case 17
                conf.Xt = M20(1,2);
                conf.Yt = M20(2,2);
                conf.Sizet = 3;
                conf.Angle = 120;   
            case 18
                conf.Xt = M20(1,2);
                conf.Yt = M20(2,2);
                conf.Sizet = 2;
                conf.Angle = 120;
            case 19
                conf.Xt = M30(1,3);
                conf.Yt = M30(2,3);
                conf.Sizet = 4;
                conf.Angle = 150;
            case 20
                conf.Xt = M30(1,3);
                conf.Yt = M30(2,3);
                conf.Sizet = 3;
                conf.Angle = 150;
            case 21 
                conf.Xt = M30(1,3);
                conf.Yt = M30(2,3);
                conf.Sizet = 2;
                conf.Angle = 150;
            case 22
                conf.Xt = M25(1,3);
                conf.Yt = M25(2,3);
                conf.Sizet = 4;
                conf.Angle = 150;
            case 23
                conf.Xt = M25(1,3);
                conf.Yt = M25(2,3);
                conf.Sizet = 3;
                conf.Angle = 150;
            case 24
                conf.Xt = M25(1,3);
                conf.Yt = M25(2,3);
                conf.Sizet = 2;
                conf.Angle = 150;
            case 25 
                conf.Xt = M20(1,3);
                conf.Yt = M20(2,3);
                conf.Sizet = 4;
                conf.Angle = 150;
            case 26
                conf.Xt = M20(1,3);
                conf.Yt = M20(2,3);
                conf.Sizet = 3;
                conf.Angle = 150;   
            case 27
                conf.Xt = M20(1,3);
                conf.Yt = M20(2,3);
                conf.Sizet = 2;
                conf.Angle = 150;                 
        end
        
        % Orientation (105-135) --> Distance (28-22) --> Size (3.5-2.5) = 8
    elseif training == 1
        load HvAOrder.mat
        condition = HvAOrderTR(pp,trial);
        switch condition
            case 1
                conf.Xt = M28(1,1);
                conf.Yt = M28(2,1);
                conf.Sizet = 3.5;
                conf.Angle = 105;
            case 2
                conf.Xt = M28(1,1);
                conf.Yt = M28(2,1);
                conf.Sizet = 2.5;
                conf.Angle = 105;
            case 3 
                conf.Xt = M22(1,1);
                conf.Yt = M22(2,1);
                conf.Sizet = 3.5;
                conf.Angle = 105;
            case 4
                conf.Xt = M22(1,2);
                conf.Yt = M22(2,2);
                conf.Sizet = 2;
                conf.Angle = 105;
            case 5
                conf.Xt = M28(1,1);
                conf.Yt = M28(2,1);
                conf.Sizet = 3.5;
                conf.Angle = 135;
            case 6
                conf.Xt = M28(1,1);
                conf.Yt = M28(2,1);
                conf.Sizet = 2.5;
                conf.Angle = 135;
            case 7 
                conf.Xt = M22(1,1);
                conf.Yt = M22(2,1);
                conf.Sizet = 3.5;
                conf.Angle = 135;
            case 8
                conf.Xt = M22(1,2);
                conf.Yt = M22(2,2);
                conf.Sizet = 2;
                conf.Angle = 135;
        end
    end
end

%% 2. Intercept: this functions calculate whenever the index finger and the thumb are pointing to the object. In that case, the function calculates the intercept coordinates and the distance between the fingers and the intercept
function [m_i,b_i,xintercept_i,yintercept_i,distance_i,m_t,b_t,xintercept_t,yintercept_t,distance_t] = intercept(xcenter,ycenter,r,data)
%% 1. Define the position for the index finger
distali = data(:,1); % Coordinates of the distal index
proximali = data(:,2); %Coordinates of the proximal index

%If there is not an increment of x, there is no m (slope), so we add some
%increment
if distali(1,1) == proximali(1,1)
   proximali(1,1) = proximali(1,1)+0.001;
end

%Define the line Y = m*X + b
m_i = (distali(2,1)-proximali(2,1))/(distali(1,1)-proximali(1,1)); %The slope of the function
b_i = distali(2,1)-m_i*distali(1,1); %The y-intercept

%Then we calculate the position (x,y) of the intercept. 
[xouti,youti] = linecirc(m_i,b_i,xcenter,ycenter,r);

%Calculate the distance between the distal marker and to the intercep
%point

d_i = sqrt ((youti(1,:)-distali(2,1)).^2+((xouti(1,:)-distali(1,1)).^2));

% In secant lines there are 2 values but we only want the smallest one
if d_i(1,1) <= d_i(1,2)
    distance_i = d_i(1,1);
    xintercept_i = xouti(1,1);
    yintercept_i = youti(1,1);
else
    distance_i = d_i(1,2);
    xintercept_i = xouti(1,2);
    yintercept_i = youti(1,2);
end

%We repeat the process for the thumb finger

%% 2. Define the position for the thumb finger
distalt = data(:,3); % Coordinates of the distal index
proximalt = data(:,4); %Coordinates of the proximal index

%If there is not an increment of x, there is no m (slope), so we add some
%increment
if distalt(1,1) == proximalt(1,1)
   proximalt(1,1) = proximalt(1,1)+0.001;
end

%Define the line Y = m*X + b
m_t = (distalt(2,1)-proximalt(2,1))/(distalt(1,1)-proximalt(1,1)); %The slope of the function
b_t = distalt(2,1)-m_t*distalt(1,1); %The y-intercept

%Then we calculate the position (x,y) of the intercept. 
[xoutt,youtt] = linecirc(m_t,b_t,xcenter,ycenter,r);

%Calculate the distance between the distal marker and to the intercep
%point

d_t = sqrt ((youtt(1,:)-distalt(2,1)).^2+((xoutt(1,:)-distalt(1,1)).^2));

% In secant lines there are 2 values but we only want the smallest one
if d_t(1,1) <= d_t(1,2)
    distance_t = d_t(1,1);
    xintercept_t = xoutt(1,1);
    yintercept_t = youtt(1,1);
else
    distance_t = d_t(1,2);
    xintercept_t = xoutt(1,2);
    yintercept_t = youtt(1,2);
end

end

