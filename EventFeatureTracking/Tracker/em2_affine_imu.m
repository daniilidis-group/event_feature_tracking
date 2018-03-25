function [feature_pos_out, transformed_points_out, scale_out, error] = em2_affine_imu(...
    events, ...
    feature_pos, ...
    template_points, ...
    template_weights, ...
    scale_in, ...
    flow, ...
    R_in, ...
    cinfo, ...
    params, ...
    fig)
%EM2_AFFINE_IMU Estimates a scale and translationbetween a set of events 
%and a set of template points, given camera rotation.
%
% EM2_AFFINE_IMU estimates an scale and translation (s, b) between a set of 
% events with known optical flow, and a set of 2D template points, given the
% camera rotation between the events and the template. This is an
% implementation of the method outlined in:
% Alex Zihao Zhu, Nikolay Atanasov, and Kostas Daniilidis. 
% "Event-Based Visual Inertial Odometry." CVPR 2017.
%
% Syntax:  [feature_pos_out, transformed_points_out, scale_out, error] = EM2_AFFINE_IMU(...
%            events, ...
%            feature_pos, ...
%            template_points, ...
%            template_weights, ...
%            scale_in, ...
%            flow, ...
%            R_in, ...
%            cinfo, ...
%            params, ...
%            fig)
%
% Inputs:
%    events           - 4xN, each column is (x,y,t,p).
%    feature_pos      - 2x1, pixel position of the feature.
%    template_points  - 2xM, set of points to align the events with.
%    template_weights - 1xM, weights for each template point.
%    scale_in         - 1x1, initialization for the scale factor.
%    flow_init        - 2x1, flow for the event window.
%    R_in             - 3x3, 3D camera rotation between the template and events.
%    cinfo            - camera info struct. Structure is defined by the ROS message:
%                       http://docs.ros.org/api/sensor_msgs/html/msg/CameraInfo.html
%    params           - parameters, defined in get_params().
%    fig              - figure handle for plotting.
%
% Outputs:
%    feature_pos_out        - 2x1, feature pos after alignment, feature_pos - b.
%    transformed_points_out - 2xN, input events in the feature window
%                               shifted by the flow A[x; y] + dt * flow + b.
%    scale_out              - 1x1, estimate for the scale factor.
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
scatter_plot_handle = [];
scale = scale_in;
scale_out = scale;
pos_shift = [0;0];

event_window = [];

error = TrackingErrors.ICPFailed;

% Transform the input events based on the 3D camera rotation between the
% input events and template points.
px = cinfo.K(3);
py = cinfo.K(6);
fx = cinfo.K(1);
fy = cinfo.K(5);

pos_u = [(feature_pos-[px;py])./[fx;fy]; 1];
pos_ur = R_in*pos_u;
pos_ur = pos_ur(1:2)/pos_ur(3);
pos_ur = pos_ur.*[fx;fy]+[px;py];

events_xy= events(1:2, :)+bsxfun(@times, flow, events(3, 1) - events(3, :));
events_xy = bsxfun(@rdivide, bsxfun(@minus, events_xy(1:2, :), [px;py]), [fx;fy]);
events_xy = [events_xy; ones(1, size(events_xy, 2))];
events_xy = R_in*events_xy;
events_xy = bsxfun(@rdivide, events_xy(1:2, :), events_xy(3, :));
events_xy = bsxfun(@plus, bsxfun(@times, events_xy, [fx;fy]), [px;py]);
rotated_points = bsxfun(@minus, events_xy, pos_ur);

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
    
    transformed_points= bsxfun(@plus, scale*rotated_points, pos_shift);
       
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
        
        rotated_points = rotated_points(:, event_window);
        transformed_points = transformed_points(:, event_window);
        
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
    
    % NOTE: 'length' is not the same as @length
    num_neighbors = cellfun('length', neighbors_cell);
    template_correspondences = cell2mat(cellfun(...
        @transpose, neighbors_cell, 'UniformOutput', false))';
    
    transformed_correspondences = repelem(1:size(transformed_points, 2), num_neighbors);
    
    weightsstacked = exp(-distancesstacked);
    
    valid_inds = sub2indc(...
        transformed_correspondences, ...
        template_correspondences, ...
        size(weights));
    
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
    
    %% Simultaneously minimize over scale and translation.
    weighted_template_barycenters = template_points*weights;
    weighted_template_barycenters = weighted_template_barycenters(:, valid_weights);

    valid_rotated_points = rotated_points(:, valid_weights);
    
    rotated_mean = mean(valid_rotated_points, 2);
    template_mean = mean(weighted_template_barycenters, 2);
    
    v_rotated_points_centered = bsxfun(@minus, valid_rotated_points, rotated_mean);
    w_template_b_centered = bsxfun(@minus, weighted_template_barycenters, template_mean);
    
    covariance = v_rotated_points_centered * w_template_b_centered';
    [U, S, V] = svd(covariance);
    if det(U*V') < 0
        S(2, 2) = -S(2, 2);
    end
    
    variance = sum(sum(v_rotated_points_centered.^2));
    new_scale = trace(S) / variance;
    
    pos_shift = template_mean - new_scale*rotated_mean;
    scale = new_scale;
    %
    %     scale = sqrt(sum(sum(...
    %         v_rotated_points_centered.*w_template_b_centered, 2)) / ...
    %         sum(sum(v_rotated_points_centered.^2)));
    
    if scale < 0.1
        error = TrackingErrors.ICPFailed;
        return
    end
    
    %% Calculate errors, plot debug information.
    err = weightsstacked*distancesstacked'/size(valid_rotated_points, 2);
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
        
        title('EM2 with IMU Rotation')
        pause(0.01)
    end
    
    if num_iter
        derr = err-errs(num_iter);
        if derr > -params.em2_params.min_err && derr <= 0
            break;
        end
    end
    
    num_iter = num_iter + 1;
end

scale_out = scale;
% Update the feature position based on the shift.
feature_pos_out = feature_pos-pos_shift;

% Compute the template events to be used in later iterations.
dt = events(3, end)-events(3, 1);
transformed_points = scale * ...
    (bsxfun(@minus, events(1:2, :), feature_pos_out + flow * dt) + ...
    bsxfun(@times, flow, (events(3, end) - events(3, :))));

% Make the window size a bit bigger for the template.
window_size = round(params.window_size * 1.5);

event_window = ...
    transformed_points(1, :) >= -window_size/2 & ...
    transformed_points(2, :) >= -window_size/2 & ...
    transformed_points(1, :) <= window_size/2 & ...
    transformed_points(2, :) <= window_size/2;

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
