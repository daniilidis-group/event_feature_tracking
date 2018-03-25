function [flow_out, time_shifted_points_out] = em1_flow(...
    events, ...
    feature_pos, ...
    flow_init, ...
    params, ...
    fig)
%EM1_FLOW Estimates optical flow for an event feature.
%
% EM1_FLOW estimates the optical flow of a set of events % combined to form 
% a 'feature'. This EM step runs from the events alone, as in:
% Alex Zihao Zhu, Nikolay Atanasov and Kostas Daniilidis.
% "Event-based Feature Tracking with Probabilistic Data Associations", ICRA 2017.
%
% Syntax:  [flow_out, time_shifted_points_out] = EM1_FLOW(...
%            events, ... 
%            feature_pos, ...
%            flow_init, ...
%            params, ...
%            fig)
%
% Inputs:
%    events               - 4xN, each column is (x,y,t,p).
%    feature_pos          - 2x1, pixel position of the feature.
%    flow_init            - 2x1, initialization for the flow.
%    params               - parameters, defined in get_params().
%    fig                  - figure handle for plotting.
%
% Outputs:
%    flow_out                - 2x1, estimated flow.
%    time_shifted_points_out - 2xN, input events in the feature window
%                              shifted by the flow [x;y] + dt * flow.
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
% Estimated flow
flow_out = [nan; nan];
% Final time shifted events.
time_shifted_points_out = [];
% Store change in flow norm at each iteration for plotting.
delta_flows = zeros(params.em1_params.max_iters, 1);
num_iter = 0;
% Plot handle.
scatter_plot_handle = [];
% Boolean array of size equal to the number of events. Is true if an event
% is in the spatial window of the feature, and thus needs to be processed.
event_window = [];
n_in_window = 0;

target_time = events(3, 1);

% Translate the events so that they are centered at the last position of
% the feature.
centered_events = events;
centered_events(1:2, :) = bsxfun(@minus, events(1:2, :), feature_pos);
prev_flow = flow;

%% Main EM Loop
while true
    if num_iter > params.em1_params.max_iters
        return
    end
    % Shift the centered events in time by the current estimate of the
    % flow.
    time_shifted_points = centered_events(1:2,:) + ...
        bsxfun(@times, flow,(target_time - centered_events(3,:)));

    % Only compute the event window once.
    if isempty(event_window)
        event_window = time_shifted_points(1, :) >= -params.window_size/2 & ...
            time_shifted_points(2, :) >= -params.window_size/2 & ...
            time_shifted_points(1, :) <= params.window_size/2 & ...
            time_shifted_points(2, :) <= params.window_size/2;
        n_in_window = sum(event_window);
                
        % Don't bother optimizing if there aren't enough events.
        if n_in_window < params.min_events_for_em
            return
        end
    
        centered_events = centered_events(:, event_window);
        time_shifted_points = time_shifted_points(:, event_window);
        %         params.max_distance = n_in_window / 400;
    end
    
    % Scale the events so that the computed distances are equal to
    % (x1-x2)/(2*sigma^2). Saves us having to divide all the
    % correspondences later.
    normalized_events = time_shifted_points' / (sqrt(2) * params.em1_params.sigma);
    
    % The KD tree allows us to compute distances between neighboring
    % events, while also performing outlier rejection. Unfortunately, as
    % the time shifted events change every iteration, it must also be
    % reconstructed at every iteration.
    kdtree = KDTreeSearcher(...
        normalized_events, ...
        'Distance', ...
        'euclidean');
    
    % Threshold distances by the Malhalanobis distance.
    [neighbors_cell, distances_cell] = rangesearch(...
        kdtree, ...
        normalized_events, ...
        params.em1_params.max_distance);

    distancesstacked = cell2mat(cellfun(...
        @transpose, distances_cell, 'UniformOutput', false))';

    % Can't solve for flow without any correspondences :(
    if isempty(distancesstacked)
        return
    end
    
    % Number of neighbors for each event.
    % NOTE: 'length' is not the same as @length
    num_neighbors_per_event = cellfun('length', neighbors_cell);
    
    neighbor_inds = cell2mat(cellfun(...
        @transpose, neighbors_cell, 'UniformOutput', false))';
    
    event_inds = repelem(1:n_in_window, num_neighbors_per_event);
    
    % The event to neighbor graph is undirected, so no need to double count
    % the event-neighbor correspondences.
    valid_correspondences = neighbor_inds > event_inds;
    
    neighbor_inds = neighbor_inds(valid_correspondences);
    event_inds = event_inds(valid_correspondences);
    distancesstacked = distancesstacked(valid_correspondences);
    
    % Distances are already scaled by the variance. Note that the original
    % equation in the paper is the sum of the product of two weights. Here
    % we simplify it with a single weight for speed.
    weights = exp(-distancesstacked);
    
    %% Simultaneously minimize over flow and translation.
    neighbor_events = centered_events(:, neighbor_inds);
    original_events = centered_events(:, event_inds);
    
    % These are the X and D matrices specified in the original paper,
    % except the weighting is handled by multiplying the D with the full
    % weight.
    X = original_events(1:2, :) - neighbor_events(1:2, :);
    D = original_events(3, :)-neighbor_events(3, :);
    weighted_D = bsxfun(@times, D, weights);
    DDT = weighted_D*D';
    XDT = X*weighted_D';
    
    flow = XDT/DDT;
    
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
        title('EM1 change in flow (convergence criterion)')
        xlim([0 params.em1_params.max_iters])
        
        subplot(2,1,2)
        if (~isempty(scatter_plot_handle))
            delete(scatter_plot_handle)
        end
        scatter_plot_handle = scatter(...
            time_shifted_points(1, :), ...
            time_shifted_points(2, :), 'r.');
        
        axis equal
        axis([-params.window_size/2-5 params.window_size/2+5 -params.window_size/2-5 params.window_size/2+5])
        axis ij
        title('EM1 Time shifted events')
        pause(0.01)
    end
    
    num_iter = num_iter + 1;
end

if params.debug
    pause(0.5)
end

flow_out = flow;

% Calculate the final shifted events to be used in EM2 later.
dt = events(3, end) - events(3, 1);
centered_events = events;
centered_events(1:2, :) = bsxfun(@minus, events(1:2, :), feature_pos + flow * dt);
time_shifted_points = centered_events(1:2, :) + ...
    bsxfun(@times, flow,(events(3, end)-centered_events(3, :)));

% Make the window a little bigger.
window_size = round(params.window_size * 1.5);

event_window = time_shifted_points(1, :) >= -window_size/2 & ...
    time_shifted_points(2, :) >= -window_size/2 & ...
    time_shifted_points(1, :) <= window_size/2 & ...
    time_shifted_points(2, :) <= window_size/2;

time_shifted_points_out = time_shifted_points(:, event_window);
end