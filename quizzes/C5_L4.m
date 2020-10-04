% Class 5. Clustering - Lesson 4. MATLAB Sensor Fusion Guided Walkthrough
% Carlos Gómez Huélamo

clc, clear all, close all;
%% Define the scenario

% Define an empty scenario

scenario = drivingScenario;
scenario.SampleTime = 0.01;

% Add a stretch of 500 m of typical highway road with two lanes

% The road is defined using a set of points, where each point defines the
% center of the road in 3D space, and a road width

roadCenters = [0 0; 50 0; 100 0; 250 20; 500 40];
roadWidth = 7.2; % Two lanes, each 3.6 meters
road(scenario, roadCenters, roadWidth);

% Create the ego vehicle and three cars around it: one that overtakes the ego vehicle 
% and passes it on the left, one that drives right in front of the ego vehicle and 
% one that drives right behind the ego vehicle. 

% All the cars follow the path defined by the road waypoints by using the path driving 
% policy. The passing car will start on the right lane, move to the left lane to pass, 
% and return to the right lane.

% Create the ego vehicle that travels at 25 m/s along the road.

egoCar = vehicle(scenario, 'ClassID', 1);
path(egoCar, roadCenters(2:end,:) - [0 1.8], 25); % On right lane

% Add a car in front of the ego vehicle.

leadCar = vehicle(scenario, 'ClassID', 1);
path(leadCar, [70 0; roadCenters(3:end,:)] - [0 1.8], 25); % On right lane

% Add a car that travels at 35 m/s along the road and passes the ego vehicle.

passingCar = vehicle(scenario, 'ClassID', 1);
waypoints = [0 -1.8; 50 1.8; 100 1.8; 250 21.8; 400 32.2; 500 38.2];
path(passingCar, waypoints, 35);

% Add a car behind the ego vehicle

chaseCar = vehicle(scenario, 'ClassID', 1);
path(chaseCar, [25 0; roadCenters(2:end,:)] - [0 1.8], 25); % On right lane

%% Define Radar

sensors = cell(6,1); % 6 radars

% Front-facing long-range radar sensor at the center of the front bumper of the car.
sensors{1} = radarDetectionGenerator('SensorIndex', 1, 'Height', 0.2, 'MaxRange', 174, ...
 'SensorLocation', [egoCar.Wheelbase + egoCar.FrontOverhang, 0], 'FieldOfView', [20, 5]);

% Rear-facing long-range radar sensor at the center of the rear bumper of the car.
sensors{2} = radarDetectionGenerator('SensorIndex', 2, 'Height', 0.2, 'Yaw', 180, 'MaxRange', 174, ...
 'SensorLocation', [-egoCar.RearOverhang, 0], 'FieldOfView', [20, 5]);

% Rear-left-facing short-range radar sensor at the left rear wheel of the car.
sensors{3} = radarDetectionGenerator('SensorIndex', 3, 'Height', 0.2, 'Yaw', 120, 'MaxRange', 30, ...
 'SensorLocation', [0, egoCar.Width/2], 'FieldOfView', [90, 5], ...
 'ReferenceRange', 50, 'AzimuthResolution', 10, 'RangeResolution', 1.25);

% Rear-right-facing short-range radar sensor at the right rear wheel of the car.
sensors{4} = radarDetectionGenerator('SensorIndex', 4, 'Height', 0.2, 'Yaw', -120, 'MaxRange', 30, ...
 'SensorLocation', [0, -egoCar.Width/2], 'FieldOfView', [90, 5], ...
 'ReferenceRange', 50, 'AzimuthResolution', 10, 'RangeResolution', 1.25);

% Front-left-facing short-range radar sensor at the left front wheel of the car.
sensors{5} = radarDetectionGenerator('SensorIndex', 5, 'Height', 0.2, 'Yaw', 60, 'MaxRange', 30, ...
 'SensorLocation', [egoCar.Wheelbase, egoCar.Width/2], 'FieldOfView', [90, 5], ...
 'ReferenceRange', 50, 'AzimuthResolution', 10, 'RangeResolution', 1.25);

% Front-left-facing short-range radar sensor at the left front wheel of the car.
sensors{6} = radarDetectionGenerator('SensorIndex', 6, 'Height', 0.2, 'Yaw', -60, 'MaxRange', 30, ...
 'SensorLocation', [egoCar.Wheelbase, -egoCar.Width/2], 'FieldOfView', [90, 5], ...
 'ReferenceRange', 50, 'AzimuthResolution', 10, 'RangeResolution', 1.25);

%% Create a MultiObjectTracker

tracker = multiObjectTracker('FilterInitializationFcn', @initSimDemoFilter, ...
 'AssignmentThreshold', 30, 'ConfirmationParameters', [4 5]);

positionSelector = [1 0 0 0; 0 0 1 0]; % Position selector
velocitySelector = [0 1 0 0; 0 0 0 1]; % Velocity selector

% Create the display and return a handle to the Bird Eye Plot

BEP = createDemoDisplay(egoCar, sensors);

%% Simulate the scenario

toSnap = true;
while advance(scenario) && ishghandle(BEP.Parent)    
    % Get the scenario time
    time = scenario.SimulationTime;

    % Get the position of the other vehicle in ego vehicle coordinates
    ta = targetPoses(egoCar);

    % Simulate the sensors
    detections = {};
    isValidTime = false(1,6);
    for i = 1:6
        [sensorDets,numValidDets,isValidTime(i)] = sensors{i}(ta, time);
        if numValidDets
            detections = [detections; sensorDets]; %#ok<AGROW>
        end
    end

    % Update the tracker if there are new detections
    if any(isValidTime)
        vehicleLength = sensors{1}.ActorProfiles.Length;
        detectionClusters = cluster_detections(detections, vehicleLength); % Function implemented in C5_L2
        confirmedTracks = updateTracks(tracker, detectionClusters, time);

        % Update bird's-eye plot
        updateBEP(BEP, egoCar, detections, confirmedTracks, positionSelector, velocitySelector);
    end

    % Snap a figure for the document when the car passes the ego vehicle
    if ta(1).Position(1) > 0 && toSnap
        toSnap = false;
        snapnow
    end
end

%% Define the Kalman Filter

function filter = initSimDemoFilter(detection)

% Use a 2-D constant velocity model to initialize a trackingKF filter.
% The state vector is [x;vx;y;vy]
% The detection measurement vector is [x;y;vx;vy]
% As a result, the measurement model is H = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

H = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
filter = trackingKF('MotionModel', '2D Constant Velocity', ...
 'State', H' * detection.Measurement, ...
 'MeasurementModel', H, ...
 'StateCovariance', H' * detection.MeasurementNoise * H, ...
 'MeasurementNoise', detection.MeasurementNoise);
end