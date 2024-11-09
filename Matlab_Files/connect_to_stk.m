% Initialize STK
app = actxserver('STK11.application');
app.Visible = 1;
root = app.Personality2;
root.NewScenario('AsteroidEarthInteraction');
scenario = root.CurrentScenario;
startTime = '1 Jul 2024 00:00:00.000';
stopTime = '2 Jul 2024 00:00:00.000';
scenario.SetTimePeriod(startTime, stopTime);
scenario.StartTime = startTime;
scenario.StopTime = stopTime;
root.Rewind;

% Create the asteroid as a satellite object
root.ExecuteCommand('New / */Satellite Asteroid');
root.ExecuteCommand(['SetState */Satellite/Asteroid Classical TwoBody J2000 "' startTime '" "' stopTime '" ' ...
                     '60 1.5e8 0.1 10 0 0 0']);  % Adjusted semi-major axis for far distance

% Use Astrogator to propagate the asteroid's trajectory
root.ExecuteCommand('Astrogator */Satellite/Asteroid ClearAll');
root.ExecuteCommand('Astrogator */Satellite/Asteroid SetPropagator Heliocentric');
root.ExecuteCommand('Astrogator */Satellite/Asteroid Propagate Astrogator End');

% Create the satellite around the asteroid
root.ExecuteCommand('New / */Satellite AsteroidSat');
root.ExecuteCommand(['SetState */Satellite/AsteroidSat Classical TwoBody J2000 "' startTime '" "' stopTime '" ' ...
                     '60 1000 0.01 30 0 0 180']);  % True anomaly set to 180 degrees for different start

% Set the asteroid mass
earthMass = 5.972e24; % Earth's mass in kg
asteroidMass = earthMass / 12000;
root.ExecuteCommand(['SetMass */Satellite/Asteroid Mass ' num2str(asteroidMass)]);

% Save and visualize the scenario
scenario.SaveAs('AsteroidEarthInteraction.sc');
root.ExecuteCommand('Animate * Reset');
root.ExecuteCommand('Animate * Start');
