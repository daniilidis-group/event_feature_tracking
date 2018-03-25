% This is the main script to run the event-based feature tracking methods
% outlined in:
%  Alex Zihao Zhu, Nikolay Atanasov and Kostas Daniilidis.
%  "Event-based Feature Tracking with Probabilistic Data Association", ICRA 2017.
%  and:
%  Alex Zihao Zhu, Nikolay Atanasov, and Kostas Daniilidis.
%  "Event-Based Visual Inertial Odometry." CVPR 2017.
% Please cite these works if you use this code in an academic publication.
%
% Author: Alex Zihao Zhu, University of Pennsylvania
% Email: alexzhu(at)seas.upenn.edu
% Copyright 2018 University of Pennsylvania 
% Alex Zihao Zhu, Nikolay Atanasov, Kostas Daniilidis
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS, CONTRIBUTORS, AND THE 
% TRUSTEES OF THE UNIVERSITY OF PENNSYLVANIA "AS IS" AND ANY EXPRESS OR 
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
% IN NO EVENT SHALL THE COPYRIGHT OWNER, CONTRIBUTORS OR THE TRUSTEES OF 
% THE UNIVERSITY OF PENNSYLVANIA BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

params = get_params();
if (exist('dataset_loaded', 'var') == 0) || ~strcmp(params.data_path, dataset_loaded)
    disp('Loading data. This may take some time.')
    clearvars -except params
    close all
    
    addpath(genpath('Tracker'))
    %% Load data
    load(params.data_path);
    load(params.undistort_map_path);
    
    [rows, cols] = size(undistort_map_x);
    
    % Read intrinsics
    fx = cinfo.K(1);
    fy = cinfo.K(5);
    px = cinfo.K(3);
    py = cinfo.K(6);
    
    pix_to_ideal = @(loc)(bsxfun(@rdivide, bsxfun(@minus, loc, [px;py]), [fx;fy]));    
    davis_t0 = events(3,1);
    
    events_zeroed = events;
    % Timestamps should start at 0.
    events_zeroed(3,:) = events_zeroed(3,:)-davis_t0;
    % Polarities should be -1 or 1 instead of 0 or 1.
    events_zeroed(4,:) = (events_zeroed(4,:)-0.5)*2;
    
    % Undistort events
    events_zeroed = undistort_events(...
        events_zeroed, ...
        undistort_map_x, ...
        undistort_map_y, ...
        [rows, cols]);
    dataset_loaded = params.data_path;
end

%% Initialization
if ~params.do_tracking
    params.debug = 0;
end

if params.debug
    params.do_parallel = 0;
end

if params.use_imu && (~exist('ang_vel', 'var') || isempty(ang_vel))
    disp('params.use_imu was set, but no angular velocity was found.')
    return
end

events_start = find(events_zeroed(3, :)> params.t_start, 1);
if params.t_end ~= -1
    events_end = find(events_zeroed(3, :) > params.t_end,1);
else
    events_end = size(events_zeroed, 2);
end
% Time since last event
event_t0 = events_zeroed(3,events_start);

imu_time_zeroed = imu_time-davis_t0;

% Find starting position in data
imu_start = find(imu_time_zeroed >= params.t_start, 1)-1;
imu_end = find(imu_time_zeroed >= params.t_end, 1);
imu_iter = imu_start;

R_curr = eye(3);

last_imu_time = imu_time_zeroed(imu_start);

time_taken = 0;

% Initialize counters
last_event_iter = events_start;
num_ima = 0;

% Count the reasons why each tracker died.
numicpmaxiters = 0;
numicpfailed = 0;
numicperr = 0;
numflow = 0;
numout = 0;
numransac = 0;

% Initialize figures so that they won't steal the focus when something is
% plotted. Use set(0, 'CurrFigure', fig) to plot to a fig.
if ~params.headless
    image_fig = figure(1);
    clf(image_fig)
    if params.debug
        debug_fig = figure(2);
        clf(debug_fig);
        set(debug_fig, 'Color',[1 1 1])
    else
        debug_fig = [];
    end
else
    image_fig = [];
    debug_fig = [];
end

id = 1;

valid_feature_positions = [];
valid_ids = [];

feature_positions = [];
% Boolean array, true if a feature is still actively tracked.
is_valid = [];
% Array of ids for each feature.
ids = [];
% Color of each feature in plots.
colors = [];

avg_int_time = params.def_int_time;
last_int_time = params.def_int_time;
% Feature tracker objects.
trackers = [];

iter = 1;
event_iter = events_start;
total_time = 0;

%% Main loop
while event_iter < events_end
    tic
    if params.constant_time
        int_time = params.def_int_time;
    else
        int_time = avg_int_time;
    end
    
    % Find end of next temporal window
    new_event_iter = find(int_time+event_t0 <= events_zeroed(3, event_iter:event_iter+params.max_events_per_window), 1);
    % Either the temporal window size is met, or max events have arrived
    if new_event_iter ~= -1
        event_iter = event_iter + new_event_iter - 1;
    else
        event_iter = event_iter + params.max_events_per_window;
    end
    
    fprintf('\nIteration %d, time %f, dt %f, num_events %d, total run time: %f\n----------------------\n', ...
        num_ima, ...
        events_zeroed(3, event_iter) - events_zeroed(3,events_start), ...
        events_zeroed(3, event_iter) - events_zeroed(3, last_event_iter), ...
        event_iter-last_event_iter, ...
        total_time);
    
    curr_events = events_zeroed(:, last_event_iter:event_iter);
      
    % Generate image
    image = accumarray(round([curr_events(2, :)' curr_events(1, :)']) + 1, ...
        1, [rows cols]);
    
    num_ima = num_ima+1;
    
    %% Tracking
    if params.do_tracking
        R_new = R_curr;
        if sum(is_valid)
            % Propagate rotation from the angular velocity.
            if params.do_ransac || params.use_imu
                [R_new, imu_iter] = propagate_rotation(...
                    imu_iter, ...
                    imu_time_zeroed, ...
                    events_zeroed(3, event_iter), ...
                    events_zeroed(3, event_iter) - events_zeroed(3, last_event_iter), ...
                    ang_vel, ...
                    R_new);
            end
                        
            fprintf('Updating %d trackers\n', sum(is_valid));
            errors = cell(1, length(trackers));
            old_feature_positions = feature_positions;

            % Update trackers
            tic
            if params.do_parallel 
                parfor i=1:length(trackers)
                    errors{i} = TrackingErrors.NotValid;
                    if trackers(i).isvalid
                        [trackers(i), ...
                            feature_positions(:, i), ...
                            is_valid(i), ...
                            error] = stepTracker(...
                            trackers(i), ...
                            curr_events, ...
                            R_new);
                        is_valid(i) = trackers(i).isvalid;
                        errors{i} = error;
                    end
                end
            else
                for i=1:length(trackers)
                    errors{i} = TrackingErrors.NotValid;
                    if trackers(i).isvalid
                        if ~is_valid(i)
                            keyboard
                        end
                        [trackers(i), ...
                            feature_positions(:, i), ...
                            is_valid(i), error] = stepTracker(...
                            trackers(i), ...
                            curr_events, ...
                            R_new);
                        is_valid(i) = trackers(i).isvalid;
                        errors{i} = error;
                    end
                end
            end
            time_taken = toc;
            fprintf('Tracking took %f seconds\n', time_taken);
            
            old_num_corners = sum(is_valid);
            
            % RANSAC
            if params.do_ransac && iter > 1
                dR = R_new * R_curr';
                
                inliers = two_point_ransac(...
                    pix_to_ideal(old_feature_positions(:, is_valid>0)),...
                    pix_to_ideal(feature_positions(:, is_valid>0)), ...
                    dR, params);
                
                is_valid(is_valid > 0) = inliers;
            end
            R_curr = R_new;
            
            % Compute the size of the next integration window and invalid
            % features that were outliers from RANSAC.
            integration_times = nan(length(trackers), 1);
            for i=1:length(integration_times)
                if is_valid(i)
                    integration_times(i) = params.integration_multiplier/norm(trackers(i).flow);
                else
                    trackers(i).isvalid = false;
                end
            end

            avg_int_time = prctile( ...
                integration_times(integration_times<1 & ~isnan(integration_times)), ...
                65);
            
            % Sanity checks.
            if isnan(avg_int_time)
                if last_int_time > 1e-4
                    avg_int_time = last_int_time;
                else
                    avg_int_time = params.def_int_time;
                end
            else
                last_int_time = avg_int_time;
            end
            
            % Count the errors for each tracker.
            for i=1:length(trackers)
                switch errors{i}
                    case TrackingErrors.Flow
                        numflow = numflow + 1;
                    case TrackingErrors.ICPMaxIters
                        numicpmaxiters = numicpmaxiters + 1;
                    case TrackingErrors.ICPError
                        numicperr= numicperr + 1;
                    case TrackingErrors.ICPFailed
                        numicpfailed= numicpfailed + 1;
                    case TrackingErrors.OutOfImg
                        numout = numout + 1;
                end
            end
            numransac = numransac + old_num_corners - sum(is_valid);
            
            fprintf('After update, %d trackers remaining. \n %d flow, %d ICP max, %d icp err, %d icp failed, %d out, %d ransac\n',...
                sum(is_valid), numflow, numicpmaxiters, numicperr, numicpfailed, numout, numransac);
        end
        
        % Detect new features if necessary.
        if sum(is_valid) < params.min_features
            feature_positions = feature_positions(:, is_valid > 0);
            ids = ids(is_valid > 0);
            detected_corners = detect_corners(image, ...
                feature_positions', ...
                params.num_features, ...
                params);
                        
            trackers = trackers(is_valid>0);
            colors = colors(is_valid>0, :);
            old_tracker_length = length(trackers);
            
            for i=1:size(detected_corners, 2)
                tracker = EventTracker(...
                    [rows, cols], ...
                    detected_corners(:, i), ...
                    id, ...
                    R_new, ...
                    cinfo, ...
                    params, ...
                    debug_fig);
                if ~tracker.isvalid
                    continue
                end
                feature_positions = [feature_positions, detected_corners(:, i)];
                trackers = [trackers, tracker];
                ids = [ids, id];
                id = id + 1;
            end
            colors = [colors; rand(length(trackers) - old_tracker_length, 3)];
            is_valid = ones(length(trackers), 1);
            
            if ~isempty(trackers)
                fprintf('Adding features, highest id is: %d\n', trackers(end).id)
            else
                fprintf('Tried to add features but none were valid.\n');
            end
        end
    end
    
    iter = iter + 1;
    valid_feature_positions = feature_positions(:, is_valid>0);
    valid_ids = ids(is_valid>0);
    %% Plotting
    if ~params.headless
        set(0, 'CurrentFigure', image_fig)
        imagesc(image)
        title(['Time: ' num2str(events_zeroed(3, event_iter)) ...
            ', Integration time: ' num2str(events_zeroed(3,event_iter)-event_t0)...
            ', Actual time: ' num2str(time_taken)])
        if params.do_tracking && sum(is_valid)
            hold on
            s = scatter(...
                feature_positions(1, is_valid>0), ...
                feature_positions(2, is_valid>0), ...
                50, colors(is_valid>0, :), 'filled', ...
                'MarkerEdgeColor', 'w');
            hold off
        end
        pause(0.001)
    end
   
    last_event_iter = event_iter+1;
    event_t0 = events_zeroed(3, event_iter);
    
    total_time = total_time + toc;
end

return
