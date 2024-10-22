%% GRASPING OBJECTS WITH A VIBROTACTILE SENSORY SUBSTITUTION GLOVE

% Program developed by Carlos de Paz (carlos.paz@uam.es,
% carlosdepazrios94@gmail.com).

clc; clear all; close all;
pp=[];
while isempty(pp)  % Create a folder for each participant
    pp = input('Please, introduce participant number: \n');
    fprintf('Participant number = %i \n', pp)
    reply = input('Is this correct? Y/N [Y]: \n','s');
    if isempty(reply)|| reply=='Y' || reply =='y'
    else
        pp=[];
    end
    if pp<10
        Pnumber=['PP00',num2str(pp)];
    elseif pp>=10 && pp<100
        Pnumber=['PP0',num2str(pp)];
    else
        Pnumber=['PP',num2str(pp)];
    end
    cd('C:\Users\Usuario\Desktop\DAVIDT\Exp3\Programs\Data')
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
% -------------------------Training-------------------------------------
mkdir('Training'); cd('Training')
 for trial=1:6 
     training = 1;
    % Display Information
    fprintf('Trial = %i \n', trial);
    load InitialConditions.mat % Condition for the subplot(1,2,1)
    
    [conf, condition] = GetIniExp3(pp,trial,training); % We set a priori the experimental conditions regarding the participant, trial and the training
    conf.q = QMC;  %Call Qualysis and get the conf for each trial

    % Print the conditions
    fprintf('Angle = %i \n', conf.Angle)
    if condition <= 6
        fprintf('Distance = 27 cm \n')
    else
        fprintf('Distance = 22 cm \n')
    end
    fprintf('Size = %i \n', conf.Sizet);
    
    
    %Calculate variables to draw the circle (object)
    alfa = 0:0.01:2*pi;
    xt = conf.Sizet.*cos(alfa)+conf.Xt; 
    yt = conf.Sizet.*sin(alfa)+conf.Yt;

    %Call the vibrotactile device (a02 = index; a03 = thumb)
    devices = daq.getDevices;
    s = daq.createSession('ni');
    s.addAnalogOutputChannel('Dev1',2:3,'Voltage'); 
    VoltageT = zeros(1,2); %index(1,1), thumb(1,2)
    
    %Define the empty vectors
    xvectori = NaN(1,100); yvectori = xvectori;
    xvectort = xvectori;   yvectort = xvectori;
    
    %Define the empty scalars
    xintercept_i = NaN;          yintercept_i = xintercept_i;
    xintercept_t = xintercept_i; yintercept_t = xintercept_i;

    [data, obj_ini] = QMC(conf.q,1); %get data from Qualisys
    
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
    
    % Dynamical plot to set the object position (match the position of the
    % virtual object with the real object)
    obj_pos = zeros(3,1);
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1])
    subplot(1,2,1)
    hold on
    axis([10 60 0 32.5])
    for i = 1:4
        plot(linspace(50,M30(1,i),100), linspace(0,M30(2,i),100),'k--', 'LineWidth', 8)
        plot(M20(1,i), M20(2,i),'ko','LineWidth', 8)
        plot(M30(1,i), M30(2,i),'ko', 'LineWidth', 8)
    end
    for i = 1:3
        plot(M22(1,i), M22(2,i),'kx', 'LineWidth', 8)
        plot(M27(1,i), M27(2,i),'kx', 'LineWidth', 8) 
    end
    plot(conf.Xt,conf.Yt, 'rx', 'LineWidth', 10)
    hold off
    if condition <= 6
        title(['Angle = ', num2str(conf.Angle), ', Dist = 27, Size = ',num2str(conf.Sizet)], 'FontSize', 32)
    else
        title(['Angle = ', num2str(conf.Angle), ', Dist = 22, Size = ',num2str(conf.Sizet)], 'FontSize', 32)
    end

    subplot(1,2,2)
    axis([10 60 0 40])
    hold on
    plot(xt,yt, 'k')
    plot(conf.Xt, conf.Yt,'kx')
    Pos_Y = plot([obj_pos(1,1), obj_pos(1,1)],[obj_pos(2,1)-conf.Sizet, obj_pos(2,1)+conf.Sizet],'r--', 'LineWidth', 10);
    Pos_X = plot([obj_pos(1,1)- conf.Sizet, obj_pos(1,1)+ conf.Sizet],[obj_pos(2,1), obj_pos(2,1)],'r--', 'LineWidth', 10);

    Centre = plot(obj_pos(1,1),obj_pos(2,1),'rx', 'Linewidth', 8); 
    while (abs(obj_pos(1,1) - conf.Xt) > .5 || abs(obj_pos(2,1) - conf.Yt) > .5)
    
        [data,obj_pos] = QMC(conf.q,1); 
        obj_pos = obj_pos./10;
        set(Pos_Y, 'Xdata', [obj_pos(1,1), obj_pos(1,1)], 'Ydata', [obj_pos(2,1)-conf.Sizet, obj_pos(2,1)+conf.Sizet]);
        set(Pos_X, 'Xdata', [obj_pos(1,1) - conf.Sizet, obj_pos(1,1)+ conf.Sizet], 'Ydata', [obj_pos(2,1), obj_pos(2,1)]);

        set(Centre, 'Xdata', obj_pos(1,1), 'Ydata', obj_pos(2,1));
        drawnow

    end
    pause(1.5)
    hold off
    close all
   
    %Parameters of the voltage's function (V = a*d^2 + b*d +c) (dmax = 30)
    a = .009;  b = - .5;   c = 10;
    
    t = 0; z = 0;%Just a varaible to control (On/Off) the while loop
    MatrixData = NaN(1,1); n = 0; %number of rows 
    fprintf('Press a key to Start the Trial \n')
    pause()
    [data,obj_ini] = QMC(conf.q,1); 
    obj_ini = obj_ini./10; data = data./10; %actualize the obj_ini
    
    tic %Calculate the time of each trial
    
    figure ('units', 'normalized', 'position', [0.05 0.35 0.35 .5])
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
    
    while ~isnan(t) %Start the trial
        
        [data,obj_pos] = QMC(conf.q,1); %get data from Qualisys

        %Check whether the data is correctly detected
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
            VoltageT(1,1) = 0;
           
        else
            xvectori = linspace(data(1,1),xintercept_i,100);
            yvectori = linspace(data(2,1),yintercept_i,100);
            VoltageT(1,1)= a.*distance_i.^2 + b.*distance_i + c;  % Quadratic function

        end
        
        % Thumb finger
        if isnan(distance_t)
            xfinalt = (max(yt)-b_t)/m_t;
            xvectort = linspace(data(1,3),xfinalt,100);
            yvectort = linspace(data(2,3),max(yt),100);
            VoltageT(1,2) = 0;
            
        else
            xvectort = linspace(data(1,3),xintercept_t,100);
            yvectort = linspace(data(2,3),yintercept_t,100);
            VoltageT(1,2)= a.*distance_t.^2 + b.*distance_t + c;  % Quadratic function
        end
            
        %After calculate the value if each data, the plot is going to update 
        set(hVI, 'Xdata', xvectori, 'Ydata', yvectori);
        set(hVT, 'Xdata', xvectort, 'Ydata', yvectort);
        set(hIntI, 'Xdata', xintercept_i, 'Ydata', yintercept_i);
        set(hIntT, 'Xdata', xintercept_t, 'Ydata', yintercept_t);
        drawnow
        hold off
        
        s.outputSingleScan(VoltageT); %Voltage for each motor
        
        %Saving the data in a matrix 
        data_save = data(1:2,1:5); data_save=data_save(:)'; %4
        n = n+1; %add a row for the MatrixData
        time = toc; %Calculate the time of each trial
        %size([data_save, m_i, b_i, distance_i, xintercept_i,yintercept_i, m_t, b_t, distance_t, xintercept_t,yintercept_t, VoltageT, time])
        MatrixData(n,1:23) = [data_save, m_i, b_i, distance_i, xintercept_i,yintercept_i, m_t, b_t, distance_t, xintercept_t,yintercept_t, VoltageT, time]; %21
        if z==1 %If Stop button is pressed, stop the experiment
            fprintf('Last Trial = %i \n', trial) 
            return
            close all
        end
        
        if z==2 %If Repeat button is pressed, save which trial needs to be repeated and start the next one
           RepeatTrials (trial,1:3) = [trial, conf.Xt, conf.Yt];
           RepeatTrials=nonzeros(RepeatTrials(:,1));
           fprintf('Repeat Trial = %i \n', trial)
           z=0;
        end
        
%         if abs(data(1,6) - conf.Xt) > 1 || abs(data(2,6) - conf.Yt) > 1 || round(data (3,6)) ~= 14 %If the participant move the object, the trial will finish (Grasped = 1; Hit = 0)
          if sum(abs(obj_pos - obj_ini))>.5 
         close all
            daq.reset()
          fprintf('Next trial')
           pause()
%          conf.Answer = input ('Has the participant grasped the target? (N = 0; Y = 1; Repeat = 3) \n');
           t=NaN; %this value will stop the loop and start the next trial
          end
        
    end
    MatrixDataTable = array2table(MatrixData, 'VariableNames', {'Xdistal_i','Ydistal_i', 'Xprox_i', 'Yprox_i','Xdista_t', 'Ydistal_t', 'Xprox_t','Yprox_t','Xwrist','Ywrist','m_i', 'b_i', 'distance_i', 'xintercept_i','yintercept_i', 'm_t', 'b_t', 'distance_t', 'xintercept_t','yintercept_t', 'VoltageT_i', 'VoltageT_t', 'time'});
    close all
    daq.reset()
%     ClosePort(PS)
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
 end
fprintf('--------------------End of the training------------------------ \n')
cd(MainCD)
% Then we repeat the code for the experiment and for the repeat trials


% --------------------------Experiment----------------------------------
clc; 
clearvars -except training alpha pp Pnumber Tnumber a b c RepeatTrials MainCD

mkdir('Experiment');cd('Experiment')

training = 0;
    
 for trial=1:80 % 16 conditions per 5 repetitions
     load InitialConditions.mat %Load the conditions for the subplot(1,2,1)

    % Display Information
    fprintf('Trial = %i \n', trial);
    %Call Qualysis and get the conf for each trial    
    [conf, condition] = GetIniExp3(pp,trial,training);
    conf.q = QMC;
    
    fprintf('Angle = %i \n', conf.Angle)
    if condition <= 8
        fprintf('Distance = 30 cm \n')
    else
        fprintf('Distance = 20 cm \n')
    end
    fprintf('Size = %i \n', conf.Sizet);
    alfa = 0:0.01:2*pi;
    
    %Parameters of the voltage's function (V = a*d^2 + b*d +c) (dmax = 28.5)
    a = .009;  b = - .5;   c = 10; %10V is too much

    %Calculate variables to draw the circle
    xt = conf.Sizet.*cos(alfa)+conf.Xt; 
    yt = conf.Sizet.*sin(alfa)+conf.Yt;

    %Call NI (a02 = index; a03 = thumb)
    devices = daq.getDevices;
    s = daq.createSession('ni');
    s.addAnalogOutputChannel('Dev1',2:3,'Voltage'); 
    VoltageT = zeros(1,2); %index(1,1), thumb(1,2)
    
    
    
    %Define the empty vectors
    xvectori = NaN(1,100); yvectori = xvectori;
    xvectort = xvectori;   yvectort = xvectori;
    
    %Define the empty scalars
    xintercept_i = NaN;          yintercept_i = xintercept_i;
    xintercept_t = xintercept_i; yintercept_t = xintercept_i;
    
    %Define the object initial position (X,Y,Z)

%     data = QMC(conf.q,1)./10;
%     obj_ini = data(:,6);

    [data, obj_ini] = QMC(conf.q,1);
    
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
        plot(M27(1,i), M27(2,i),'kx', 'LineWidth', 6) 
    end
    plot(conf.Xt,conf.Yt, 'rx', 'LineWidth', 6)
    if condition <= 8
        title(['Angle = ', num2str(conf.Angle), ', Dist = 30, Size = ',num2str(conf.Sizet)], 'FontSize', 32)
    else
        title(['Angle = ', num2str(conf.Angle), ', Dist = 20, Size = ',num2str(conf.Sizet)], 'FontSize', 32)
    end
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
    
%     figure (1)
    figure ('units', 'normalized', 'position', [0.05 0.35 0.35 .5])
    title(['Trial = ', num2str(trial)], 'FontSize', 14)
    hold on
%     axis([10 30 60 85])
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
        
       
        %Intercept function
        [data,m_i,b_i,xintercept_i,yintercept_i,distance_i,m_t,b_t,xintercept_t,yintercept_t,distance_t] = intercept(obj_ini(1,1),obj_ini(2,1),conf.Sizet,data);
         
        % Index finger
        if isnan(distance_i) %If there is not an intercept, just plot a line from xdistal_i to xfinali (X position for conf.Yt) 
            xfinali = (max(yt)-b_i)/m_i;
            xvectori = linspace(data(1,1),xfinali,100);
            yvectori = linspace(data(2,1),max(yt),100);
            VoltageT(1,1) = 0;
           
        else
            xvectori = linspace(data(1,1),xintercept_i,100);
            yvectori = linspace(data(2,1),yintercept_i,100);
           %VoltageT(1,1) = 5-((distance_i*7)/(conf.Yt-50))+3; %Lineal function
            VoltageT(1,1)= a.*distance_i.^2 + b.*distance_i + c;  % Quadratic function
        end
        
        % Thumb finger
        if isnan(distance_t)
            xfinalt = (max(yt)-b_t)/m_t;
            xvectort = linspace(data(1,3),xfinalt,100);
            yvectort = linspace(data(2,3),max(yt),100);
            VoltageT(1,2) = 0;
            
        else
            xvectort = linspace(data(1,3),xintercept_t,100);
            yvectort = linspace(data(2,3),yintercept_t,100);
           %VoltageT(1,2) = 5-((distance_t*7)/(conf.Yt-50))+3; %Lineal function
            VoltageT(1,2)= a.*distance_t.^2 + b.*distance_t + c; % Quadratic function
        end
            
        %After calculate the value if each data, the plot is going to update 
        set(hVI, 'Xdata', xvectori, 'Ydata', yvectori);
        set(hVT, 'Xdata', xvectort, 'Ydata', yvectort);
        set(hIntI, 'Xdata', xintercept_i, 'Ydata', yintercept_i);
        set(hIntT, 'Xdata', xintercept_t, 'Ydata', yintercept_t);
        drawnow
        hold off
        
        s.outputSingleScan(VoltageT); %Voltage for each motor
        
        %Saving the data in a matrix 
        data_save = data(1:2,1:5); data_save=data_save(:)'; %4 finger markers + wrist
        n = n+1; %add a row for the MatrixData
        time = toc; %Calculate the time of each trial
        MatrixData(n,1:23) = [data_save, m_i, b_i, distance_i, xintercept_i,yintercept_i, m_t, b_t, distance_t, xintercept_t,yintercept_t, VoltageT, time]; %21
        
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
    MatrixDataTable = array2table(MatrixData, 'VariableNames', {'Xdistal_i','Ydistal_i', 'Xprox_i', 'Yprox_i','Xdista_t', 'Ydistal_t', 'Xprox_t','Yprox_t','Xwrist','Ywrist','m_i', 'b_i', 'distance_i', 'xintercept_i','yintercept_i', 'm_t', 'b_t', 'distance_t', 'xintercept_t','yintercept_t', 'VoltageT_i', 'VoltageT_t', 'time'});
    close all
    daq.reset()
%     ClosePort(PS)
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
% ------------------------ Repeat the trials----------------------------
    clearvars -except training alpha pp Pnumber Tnumber a b c RepeatTrials MainCD

if ~isempty(RepeatTrials)
i=0;    
mkdir('Recovery'); cd('Recovery')
    
    for i = 1:length(RepeatTrials)
        load InitialConditions.mat 

        fprintf('Repeated Trial = %i \n', RepeatTrials(i,1));
        training = 0;
%         PS = OpenPort;
%         Mode = 1;
        [conf, condition] = GetIniExp3(pp, RepeatTrials(i,1),training);
        conf.q = QMC;
        
        fprintf('Angle = %i \n', conf.Angle)
        if condition <= 8
            fprintf('Distance = 30 cm')
        else
            fprintf('Distance = 20 cm')
        end
        fprintf('Size = %i \n', conf.Sizet);
        
    
        %Calculate variables to draw the circle
        trial = RepeatTrials(i,1);
        alfa = 0:0.01:2*pi;
        xt = conf.Sizet.*cos(alfa)+conf.Xt; 
        yt = conf.Sizet.*sin(alfa)+conf.Yt;

        %Call NI (a02 = index; a03 = thumb)
        devices = daq.getDevices;
        s = daq.createSession('ni');
        s.addAnalogOutputChannel('Dev1',2:3,'Voltage'); 
        VoltageT = zeros(1,2); %index(1,1), thumb(1,2)

        

        %Define the empty vectors
        xvectori = NaN(1,100); yvectori = xvectori;
        xvectort = xvectori;   yvectort = xvectori;

        %Define the empty scalars
        xintercept_i = NaN;          yintercept_i = xintercept_i;
        xintercept_t = xintercept_i; yintercept_t = xintercept_i;

        %Define the object initial position (X,Y,Z)
        [data, obj_ini] = QMC(conf.q,1);

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
         obj_pos = zeros(3,1);

        figure('units', 'normalized', 'outerposition', [0 0 1 1])
        subplot(1,2,1)
        hold on
        axis([10 60 0 32.5])
        for j = 1:4
            plot(linspace(50,M30(1,j),100), linspace(0,M30(2,j),100),'k--', 'LineWidth', 6)
            plot(M20(1,j), M20(2,j),'ko', 'LineWidth', 6)
            plot(M30(1,j), M30(2,j),'ko', 'LineWidth', 6)
        end
        for j = 1:3
            plot(M22(1,j), M22(2,j),'kx', 'LineWidth', 6)
            plot(M27(1,j), M27(2,j),'kx', 'LineWidth', 6) 
        end
        plot(conf.Xt,conf.Yt, 'rx', 'LineWidth', 6)
        if condition <= 8
            title(['Angle = ', num2str(conf.Angle), ', Dist = 30, Size = ',num2str(conf.Sizet)], 'FontSize', 32)
        else
            title(['Angle = ', num2str(conf.Angle), ', Dist = 20, Size = ',num2str(conf.Sizet)], 'FontSize', 32)
        end

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
        close all
        pause(1.5)
        [data,obj_ini] = QMC(conf.q,1); 
        obj_ini = obj_ini./10; data = data./10;
        conf.Xobj = obj_ini(1,1);
        conf.Yobj = obj_ini(2,1);
        %Parameters of the voltage's function (V = a*d^2 + b*d +c) (dmax = 30)
        a = .009;  b = - .5; c = 10;
        t = 0; z = 0;%Just a varaible to control (On/Off) the while loop
        MatrixData = NaN(1,1); n = 0; %number of rows 
        fprintf('Press a key to Start the Trial \n')
        pause()
        
        tic %Calculate the time of each trial
        figure (1)
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
            
            %Intercept function
%             [data,m_i,b_i,xintercept_i,yintercept_i,distance_i,m_t,b_t,xintercept_t,yintercept_t,distance_t] = intercept(conf.Xt,conf.Yt,conf.Sizet,data);
            [data,m_i,b_i,xintercept_i,yintercept_i,distance_i,m_t,b_t,xintercept_t,yintercept_t,distance_t] = intercept(obj_ini(1,1),obj_ini(2,1),conf.Sizet,data);

            % Index finger
            if isnan(distance_i) %If there is not an intercept, just plot a line from xdistal_i to xfinali (X position for conf.Yt) 
                xfinali = (max(yt)-b_i)/m_i;
                xvectori = linspace(data(1,1),xfinali,100);
                yvectori = linspace(data(2,1),max(yt),100);
                VoltageT(1,1) = 0;
           
            else
                xvectori = linspace(data(1,1),xintercept_i,100);
                yvectori = linspace(data(2,1),yintercept_i,100);
               %VoltageT(1,1) = 5-((distance_i*7)/(conf.Yt-50))+3; %Lineal function
                VoltageT(1,1)= a.*distance_i.^2 + b.*distance_i + c;  % Quadratic function
            end

            % Thumb finger
            if isnan(distance_t)
                xfinalt = (max(yt)-b_t)/m_t;
                xvectort = linspace(data(1,3),xfinalt,100);
                yvectort = linspace(data(2,3),max(yt),100);
                VoltageT(1,2) = 0;

            else
                xvectort = linspace(data(1,3),xintercept_t,100);
                yvectort = linspace(data(2,3),yintercept_t,100);
               %VoltageT(1,2) = 5-((distance_t*7)/(conf.Yt-50))+3; %Lineal function
                VoltageT(1,2)= a.*distance_t.^2 + b.*distance_t + c; % Cuadratic function
            end

            %After calculate the value if each data, the plot is going to update 
            set(hVI, 'Xdata', xvectori, 'Ydata', yvectori);
            set(hVT, 'Xdata', xvectort, 'Ydata', yvectort);
            set(hIntI, 'Xdata', xintercept_i, 'Ydata', yintercept_i);
            set(hIntT, 'Xdata', xintercept_t, 'Ydata', yintercept_t);
            drawnow
            hold off

            s.outputSingleScan(VoltageT); %Voltage for each motor

            %Saving the data in a matrix 
            data_save = data(1:2,1:5); data_save=data_save(:)';
            n = n+1; %add a row for the MatrixData
            time = toc; %Calculate the time of each trial
            MatrixData(n,1:23) = [data_save, m_i, b_i, distance_i, xintercept_i,yintercept_i, m_t, b_t, distance_t, xintercept_t,yintercept_t, VoltageT, time];

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
        MatrixDataTable = array2table(MatrixData, 'VariableNames', {'Xdistal_i','Ydistal_i', 'Xprox_i', 'Yprox_i','Xdista_t', 'Ydistal_t', 'Xprox_t','Yprox_t','Xwrist','Ywrist','m_i', 'b_i', 'distance_i', 'xintercept_i','yintercept_i', 'm_t', 'b_t', 'distance_t', 'xintercept_t','yintercept_t', 'VoltageT_i', 'VoltageT_t', 'time'});
        close all
        daq.reset()
    %     ClosePort(PS)
        QMC(conf.q, 'quit')
    %     clear PS

        if RepeatTrials(i,1)<10
            Tnumber = ['RepitedTR00', num2str(RepeatTrials(i,1))];
        else
            Tnumber = ['RepitedTR0', num2str(RepeatTrials(i,1))];
        end
        FileID = [Pnumber,Tnumber];

        save(FileID, 'conf', 'MatrixData', 'MatrixDataTable');

        fprintf('End of the trial \n')
    end
cd(MainCD)
fprintf('------------------End of the Repeated Trials---------------- \n')
end
%% Auxiliar functions
%% 1. GetIniExp3: this functions set the coordinates of the to-be-grasped object

function [ conf, condition ] = GetIniExp3(pp,trial,training)
    load InitialConditions.mat 
    if training == 0
        load Exp3Order.mat
        condition = Exp3Order(pp,trial);
        switch condition
            case 1
                conf.Xt = M30(1,1);
                conf.Yt = M30(2,1);
                conf.Sizet = 2.5;
                conf.Angle = 90;
            case 2
                conf.Xt = M30(1,1);
                conf.Yt = M30(2,1);
                conf.Sizet = 3.5;
                conf.Angle = 90;

            case 3 
                conf.Xt = M30(1,2);
                conf.Yt = M30(2,2);
                conf.Sizet = 2.5;
                conf.Angle = 110;
            case 4
                conf.Xt = M30(1,2);
                conf.Yt = M30(2,2);
                conf.Sizet = 3.5;
                conf.Angle = 110;
            case 5
                conf.Xt = M30(1,3);
                conf.Yt = M30(2,3);
                conf.Sizet = 2.5;
                conf.Angle = 130;
            case 6
                conf.Xt = M30(1,3);
                conf.Yt = M30(2,3);
                conf.Sizet = 3.5;
                conf.Angle = 130;
            case 7
                conf.Xt = M30(1,4);
                conf.Yt = M30(2,4);
                conf.Sizet = 2.5;
                conf.Angle = 150;
            case 8
                conf.Xt = M30(1,4);
                conf.Yt = M30(2,4);
                conf.Sizet = 3.5;
                conf.Angle = 150;
            case 9
                conf.Xt = M20(1,1);
                conf.Yt = M20(2,1);
                conf.Sizet = 2.5;
                conf.Angle = 90;
            case 10
                conf.Xt = M20(1,1);
                conf.Yt = M20(2,1);
                conf.Sizet = 3.5;
                conf.Angle = 90;
            case 11 
                conf.Xt = M20(1,2);
                conf.Yt = M20(2,2);
                conf.Sizet = 2.5;
                conf.Angle = 110;
            case 12
                conf.Xt = M20(1,2);
                conf.Yt = M20(2,2);
                conf.Sizet = 3.5;
                conf.Angle = 110;
            case 13
                conf.Xt = M20(1,3);
                conf.Yt = M20(2,3);
                conf.Sizet = 2.5;
                conf.Angle = 130;
            case 14
                conf.Xt = M20(1,3);
                conf.Yt = M20(2,3);
                conf.Sizet = 3.5;
                conf.Angle = 130;
            case 15
                conf.Xt = M20(1,4);
                conf.Yt = M20(2,4);
                conf.Sizet = 2.5;
                conf.Angle = 150;
            case 16
                conf.Xt = M20(1,4);
                conf.Yt = M20(2,4);
                conf.Sizet = 3.5;
                conf.Angle = 150;
        end
    elseif training == 1
        load Exp3OrderTr.mat
        condition = Exp3OrderTR(pp,trial);
        switch condition
            case 1
                conf.Xt = M27(1,1);
                conf.Yt = M27(2,1);
                conf.Sizet = 2;
                conf.Angle = 100;
            case 2
                conf.Xt = M27(1,1);
                conf.Yt = M27(2,1);
                conf.Sizet = 4;
                conf.Angle = 100;
            case 3 
                conf.Xt = M27(1,2);
                conf.Yt = M27(2,2);
                conf.Sizet = 2;
                conf.Angle = 120;
            case 4
                conf.Xt = M27(1,2);
                conf.Yt = M27(2,2);
                conf.Sizet = 4;
                conf.Angle = 120;
            case 5
                conf.Xt = M27(1,3);
                conf.Yt = M27(2,3);
                conf.Sizet = 2;
                conf.Angle = 140;
            case 6
                conf.Xt = M27(1,3);
                conf.Yt = M27(2,3);
                conf.Sizet = 4;
                conf.Angle = 140;
            case 7
                conf.Xt = M22(1,1);
                conf.Yt = M22(2,1);
                conf.Sizet = 2;
                conf.Angle = 100;
            case 8
                conf.Xt = M22(1,1);
                conf.Yt = M22(2,1);
                conf.Sizet = 4;
                conf.Angle = 100;
            case 9
                conf.Xt = M22(1,2);
                conf.Yt = M22(2,2);
                conf.Sizet = 2;
                conf.Angle = 120;
            case 10
                conf.Xt = M22(1,2);
                conf.Yt = M22(2,2);
                conf.Sizet = 4;
                conf.Angle = 120;
            case 11
                conf.Xt = M22(1,3);
                conf.Yt = M22(2,3);
                conf.Sizet = 2;
                conf.Angle = 140;
            case 12
                conf.Xt = M22(1,3);
                conf.Yt = M22(2,3);
                conf.Sizet = 4;
                conf.Angle = 140;
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


