function [flow_out, time_shifted_points_out] = em1_flow_with_prior(...
    events, ...
    feature_pos, ...
    prev_shifted_points, ...
    prev_shifted_weights, ...
    flow_init, ...
    params, ...
    fig)
%EM1_FLOW_WITH_PRIOR Estimates optical flow for an event feature.
%
% EM1_FLOW_WITH_PRIOR estimates the optical flow of a set of events
% combined to form a 'feature'. This EM step uses the time shifted points
% from the previous iteration as a template, as in:
% Alex Zihao Zhu, Nikolay Atanasov and Kostas Daniilidis.
% "Event-based Visual Inertial Odometry", 
% IEEE International Conference on Computer Vision and Pattern Recognition (CVPR), 2017.
%
% Syntax:  [flow_out, time_shifted_points_out] = EM1_FLOW_WITH_PRIOR(...
%            events, ... 
%            feature_pos, ...
%            prev_shifted_points, ...
%            prev_shifted_weights, ...
%            flow_init, ...
%            params, ...
%            fig)
%
% Inputs:
%    events               - 4xN, each column is (x,y,t,p).
%    feature_pos          - 2x1, pixel position of the feature.
%    prev_shifted_points  - 2xM, time shifted points from the previous
%                           iteration, used as a template.
%    prev_shifted_weights - 1xM, weights for each shifted point.
%    flow_init            - 2x1, initialization for the flow.
%    params               - parameters, defined in get_params().
%    fig                  - figure handle for plotting.
%
% Outputs:
%    flow_out                - 2x1, estimated flow.
%    time_shifted_points_out - 2xN, input events in the feature window
%                              shifted by the flow [x; y] + dt * flow.
%
% See also GET_PARAMS
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

%% Initialization
flow = flow_init;
prev_flow = flow;
% Estimated flow
flow_out = [nan; nan];
% Estimated canonical events
time_shifted_points_out = [];
delta_flows = zeros(params.em1_params.max_iters, 1);
num_iter = 0;
scatter_plot_handle = [];
event_window = [];

target_time = events(3, 1);

centered_events = events;
centered_events(1:2, :) = bsxfun(@minus, events(1:2, :), feature_pos);

kdtree = KDTreeSearcher(...
    prev_shifted_points' / (sqrt(2) * params.em1_params.sigma), ...
    'Distance', 'euclidean');

if params.debug
    set(0, 'CurrentFigure', fig)
    subplot(2,1,2)
    if ~isempty(prev_shifted_points)
        scatter(prev_shifted_points(1, :), prev_shifted_points(2, :), prev_shifted_weights * 10, 'b.')
    end
end

%% Main EM Loop
while true
    if num_iter > params.em1_params.max_iters
        return
    end
    
    time_shifted_points = centered_events(1:2,:) + ...
        bsxfun(@times, flow,(target_time-centered_events(3,:)));

    if isempty(event_window)
        event_window = ...
            time_shifted_points(1, :) >= -params.window_size/2 &...
            time_shifted_points(2, :) >= -params.window_size/2 & ...
            time_shifted_points(1, :) <= params.window_size/2 & ...
            time_shifted_points(2, :) <= params.window_size/2;
        
        n_in_window = sum(event_window);
        if n_in_window < params.min_events_for_em
            return
        end
        
        centered_events = centered_events(:, event_window);
        time_shifted_points = time_shifted_points(:, event_window);
        
        weights = zeros(size(prev_shifted_points, 2), size(time_shifted_points, 2));
    end
    
    [neighbors_cell, distances_cell] = rangesearch(...
        kdtree, ...
        time_shifted_points' / (sqrt(2) * params.em1_params.sigma), ...
        params.em1_params.max_distance);

    distancesstacked = cell2mat(cellfun(...
        @transpose, distances_cell, 'UniformOutput', false))';
    
    if isempty(distancesstacked)
        return
    end
    
    % NOTE: 'length' is not the same as @length
    num_neighbors = cellfun('length', neighbors_cell);
    prev_correspondences = cell2mat(cellfun(...
        @transpose, neighbors_cell, 'UniformOutput', false))';
    curr_correspondences = repelem(1:size(time_shifted_points, 2), num_neighbors);
    
    % It's cheaper to multiply to 0 to reset the weights matrix than to
    % reinitialize it using zeros.
    weightsstacked = exp(-distancesstacked); 
    weights = weights * 0;
    
    valid_inds = sub2indc(...
        curr_correspondences, ...
        prev_correspondences, ...
        size(weights));
    
    weights(valid_inds) = weightsstacked;
    weights = bsxfun(@times, prev_shifted_weights, weights);
    weights_sum = sum(weights, 1);
    valid_weights = weights_sum > 0;
    weights = bsxfun(@rdivide, weights, weights_sum + 1e-10);
        
    % Simultaneously minimize over flow and translation
    weighted_prior_points = prev_shifted_points * weights;
    weighted_prior_points = weighted_prior_points(:, valid_weights);
    
    valid_centered_events = centered_events(:, valid_weights);
    dx = weighted_prior_points(1:2, :) - valid_centered_events(1:2, :);
    dt = target_time - valid_centered_events(3, :);
    
    flow = (dx*dt') / (dt*dt');
    
    %% Calculate change in flow, plot debug information.
    if (norm(flow - prev_flow) < params.em1_params.min_err)
        break;
    end
    
    delta_flows(num_iter+1) = norm(flow - prev_flow);
    prev_flow = flow;
    
    if params.debug
        set(0, 'CurrentFigure', fig)
        subplot(2,1,1)
        plot(delta_flows(delta_flows > 0),'b')
        title('EM1 with prior change in flow (convergence criterion)')
        xlim([0 params.em1_params.max_iters])
        
        subplot(2,1,2)
        if (~isempty(scatter_plot_handle))
            delete(scatter_plot_handle)
        end
        
        if ~isempty(prev_shifted_points)
            hold on
        end
        
        scatter_plot_handle = scatter(time_shifted_points(1, :), time_shifted_points(2, :),'r.');
        hold off
        axis equal
        axis([-params.window_size/2-5 params.window_size/2+5 -params.window_size/2-5 params.window_size/2+5])
        axis ij
        
        title('EM1 with prior time shifted events')
        pause(0.01)
    end
    
    num_iter = num_iter + 1;
end

if params.debug
    pause(0.5)
end

flow_out = flow;

dt = events(3, end) - events(3, 1);
centered_events = bsxfun(@minus, events(1:2, :), feature_pos + flow * dt);
time_shifted_points = centered_events(1:2, :) + ...
    bsxfun(@times, flow, (events(3, end) - events(3,:)));


% Make the window a little bigger.
window_size = round(params.window_size * 1.5);

event_window = time_shifted_points(1, :) >= -window_size/2 & ... 
    time_shifted_points(2, :) >= -window_size/2 & ...
    time_shifted_points(1, :) <= window_size/2 & ...
    time_shifted_points(2, :) <= window_size/2;


time_shifted_points_out = time_shifted_points(:, event_window);
end