% Class 5. Clustering - Lesson 3. Kalman Tracking
% Carlos G�mez Hu�lamo

filter = trackingKF('MotionModel', model, ...
                    'State', state, ...
                    'MeasurementModel',measurementModel, ...
                    'StateCovariance', stateCovrariance, ...
                    'MeasurementNoise', measurementNoise);