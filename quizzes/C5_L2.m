% Class 5. Clustering - Lesson 2. Clustering
% Carlos Gómez Huélamo

function detection_clusters = cluster_detections(detections, vehicle_size)

    N = numel(detections); % Number of detections
    distances = zeros(N); 
    
    for i = 1:N
        for j = i+1:N
            if detections{i}.SensorIndex == detections{j}.SensorIndex % If the detections are from the same sensor
                distances(i,j) = norm(detections{i}.Measurement(1:2) - detections{j}.Measurement(1:2)); % Euclidean distance
            else
                distances(i,j) = inf;
            end
        end
    end
    
    left_to_check = 1:N;
    i = 0;
    detection_clusters = cell(N,1);
    
    while ~isempty(left_to_check)
        % Remove the detections that are in the same cluster as the one
        % under consideration
        
        under_consideration = left_to_check(1); % Take first element under consideration
        
        cluster_indices = (distances(under_consideration, left_to_check) < vehicle_size);
        
        detection_indices = left_to_check(clister_indices);
        cluster_detections = [detections{detection_indices}];
        cluster_measurements = [cluster_detections.Measurement];
        
        measurements = mean(cluster_measurements, 2);
        measurements_2D = [measurements(1:2); meas(4:5)];
        
        i = i + 1; % Create a new cluster id
        detection_clusters{i} = detections{detection_indices(1)};
        detection_clusters{i}.Measurement = measurements_2D;
        left_to_check(cluster_indices) = []; % Erase this index     
    end
    
    detection_clusters(i+1:end) = [];
end

