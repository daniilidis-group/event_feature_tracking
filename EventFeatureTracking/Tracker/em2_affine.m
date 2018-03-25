function [feature_pos_out, transformed_points_out, A_out, error] = em2_affine(...
    events, ...
    feature_pos, ...
    template_points, ...
    template_weights, ...
    A_in, ...
    flow, ...
    params, ...
    fig)
%EM2_AFFINE Estimates an affine transform between a set of events and a set
%of template points.
%
% EM2_AFFINE estimates an affine warp (A, b) between a set of events with 
% known optical flow, and a set of 2D template points. This EM step runs 
% from the events alone, as in:
% Alex Zihao Zhu, Nikolay Atanasov and Kostas Daniilidis.
% "Event-based Feature Tracking with Probabilistic Data Association", 
% IEEE International Conference on Robotics and Automation (ICRA), 2017.
%
% Syntax:  [feature_pos_out, transformed_points_out, A_out, error] = EM2_AFFINE(...
%            events, ...
%            feature_pos, ...
%            template_points, ...
%            template_weights, ...
%            A_in, ...
%            flow, ...
%            params, ...
%            fig)
%
% Inputs:
%    events           - 4xN, each column is (x,y,t,p).
%    feature_pos      - 2x1, pixel position of the feature.
%    template_points  - 2xM, set of points to align the events with.
%    template_weights - 1xM, weights for each template point.
%    A_in             - 2x2, initialization for the affine warp.
%    flow_init        - 2x1, flow for the event window.
%    params           - parameters, defined in get_params().
%    fig              - figure handle for plotting.
%
% Outputs:
%    feature_pos_out        - 2x1, feature pos after alignment, feature_pos - b.
%    transformed_points_out - 2xN, input events in the feature window
%                               shifted by the flow A[x; y] + dt * flow + b.
%    A_out                  - 2x2, estimate for the affine warp.
%    error                  - Error code, defined in TrackingErrors.
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
% Estimated canonical events
transformed_points_out = [];
% New estimated patch position if estimated
feature_pos_out = [nan; nan];
errs = zeros(params.em2_params.max_iters, 1);
num_iter = 0;
error = TrackingErrors.ICPFailed;
scatter_plot_handle = [];

target_time = events(3, 1);
A = A_in;
A_out = A;
pos_shift = [0;0];

event_window = [];

time_shifted_points = bsxfun(@minus, events(1:2, :), feature_pos) + ...
    bsxfun(@times, flow, (target_time - events(3,:)));

if isempty(template_points)
    return
end

kdtree = KDTreeSearcher(...
    template_points' / (sqrt(2) * params.em2_params.sigma), ...
    'Distance', 'euclidean');

while true
    if num_iter > params.em2_params.max_iters
        error = TrackingErrors.ICPMaxIters;
        return
    end
    
    transformed_points = bsxfun(@plus, A * time_shifted_points, pos_shift);
    
    if isempty(event_window)
        event_window = ...
            transformed_points(1, :) >= -params.window_size/2 & ...
            transformed_points(2, :) >= -params.window_size/2 & ...
            transformed_points(1, :) <= params.window_size/2 & ...
            transformed_points(2, :) <= params.window_size/2;
                
        n_in_window = sum(event_window);
        if n_in_window < params.min_events_for_em
            error = TrackingErrors.ICPFailed;
            return
        end
        
        time_shifted_points = time_shifted_points(:, event_window);
        transformed_points = transformed_points(:, event_window);
        %         if n_in_window > 200
        %             shifted_events = shifted_events(:, 1:200);
        %             transformed_events = transformed_events(:, 1:200);
        %             n_in_window = 200;
        %         end
        
        weights = zeros(size(template_points, 2), size(transformed_points, 2));
    end
    
    [neighbors_cell, distances_cell] = rangesearch(...
        kdtree, ...
        transformed_points' / (sqrt(2) * params.em2_params.sigma), ...
        params.em2_params.max_distance);
    
    distancesstacked = cell2mat(cellfun(...
        @transpose, distances_cell, 'UniformOutput', false))';
    
    if isempty(distancesstacked)
        error = TrackingErrors.ICPFailed;
        return
    end
    
    num_neighbors = cellfun('length', neighbors_cell);
    template_correspondences = cell2mat(cellfun(...
        @transpose, neighbors_cell, 'UniformOutput', false))';
    
    transformed_correspondences = repelem(1:n_in_window, num_neighbors);
    
    weightsstacked = exp(-distancesstacked);
    
    valid_inds = sub2indc(transformed_correspondences, template_correspondences, size(weights));
    
    % It's cheaper to multiply to 0 to reset the weights matrix than to
    % reinitialize it using zeros.
    weights = weights * 0;
    weights(valid_inds) = weightsstacked;
    
    if ~isempty(template_weights)
        weights = bsxfun(@times, template_weights, weights);
    end
    
    weight_sum = sum(weights, 1);
    valid_weights = (weight_sum > 0);
    
    weights = bsxfun(@rdivide, weights, weight_sum + 1e-10);
    
    weightsstacked = weights(valid_inds);
    
    %% Simultaneously minimize over affine warp and translation.
    weighted_template_points = template_points*weights;
    weighted_template_points = weighted_template_points(:, valid_weights);
    
    current_events = time_shifted_points(:, valid_weights);
    
    D = [current_events(1:2, :); ones(1, size(current_events, 2))];
    X = weighted_template_points;
    
    transform = X/D;
        
    A = transform(1:2, 1:2);
    scale = sqrt(abs(det(A)));
    pos_shift = (transform(:, 3));
    
    if scale < params.em2_params.min_scale
        error = TrackingErrors.ICPFailed;
        return
    end
    
    %% Calculate errors, plot debug information.
    err = weightsstacked*distancesstacked'/ size(current_events, 2);
    errs(num_iter+1) = err;
    
    if params.debug
        set(0, 'CurrentFigure', fig)
        subplot(2,1,1)
        plot(errs(errs > 0), 'b')
        title('Change in error (convergence criterion)')
        xlim([0 params.em2_params.max_iters])
        subplot(2,1,2)
        if (~isempty(scatter_plot_handle))
            delete(scatter_plot_handle)
        end        
        scatter_plot_handle = scatter(transformed_points(1, :), transformed_points(2, :),'r.');
        hold on
        scatter(template_points(1, :), template_points(2, :),template_weights*10,'b.');
        hold off
        axis equal
        axis([-params.window_size/2-5, ...
            params.window_size/2+5, ...
            -params.window_size/2-5, ...
            params.window_size/2+5])
        axis ij
        title('EM 2')

        pause(0.01)
    end
    
    if num_iter
        derr = err - errs(num_iter);
        if derr > -params.em2_params.min_err && derr <= 0
            break;
        end
    end
    num_iter = num_iter + 1;
end

A_out = A;
% Update the feature position based on the shift.
feature_pos_out = feature_pos - pos_shift;

% Compute the template events to be used in later iterations.
dt = events(3, end)-events(3, 1);
time_shifted_points = bsxfun(@minus, events(1:2, :), ...
    feature_pos + flow * dt) + ...
    bsxfun(@times, flow, (events(3, end) - events(3, :)));

transformed_points = A_out*time_shifted_points;

% Make the window size a bit bigger for the template.
window_size = round(params.window_size * 1.5);

event_window = ...
    transformed_points(1, :) >= -window_size/2 & ...
    transformed_points(2, :) >= -window_size/2 & ...
    transformed_points(1, :) < window_size/2 & ...
    transformed_points(2, :) < window_size/2;

transformed_points_out = transformed_points(:, event_window);

if params.debug
    pause(0.5)
end

% Measure how many transformed points are considered 'outliers' when
% matched against the templates. Outliers are defined as points without any
% template points within 2 pixels. Note that this function still requires
% more work to be properly robust.
[~, distances_cell] = knnsearch(...
    kdtree, ...
    transformed_points_out'/(sqrt(2)*params.em2_params.sigma), ...
    'K', 1);
distances_cell = distances_cell * 2 * params.em2_params.sigma^2;
percent_outliers = sum(distances_cell > 2.0) / length(distances_cell);

if percent_outliers > params.em2_params.max_outlier_thresh
    error = TrackingErrors.ICPError;
    feature_pos_out = [nan; nan];
    return;
end

end
