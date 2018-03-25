classdef EventTracker
    %EventTracker implements feature tracking for event cameras.
    %  This class implements the methods described in:
    %  Alex Zihao Zhu, Nikolay Atanasov and Kostas Daniilidis.
    %  "Event-based Feature Tracking with Probabilistic Data Association", 
    %  IEEE International Conference on Robotics and Automation (ICRA), 2017.
    %  and:
    %  Alex Zihao Zhu, Nikolay Atanasov, and Kostas Daniilidis. 
    %  "Event-Based Visual Inertial Odometry." 
    %  IEEE International Conference on Computer Vision and Pattern Recognition (CVPR), 2017.
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
    
    properties
        feature_pos % Pixel position of the feature.
        num_frames % Number of frames this feature's been tracked.
        isvalid % If the feature is still valid.
        template_points % Set of template points for EM2.
        template_weights % Set of template weights for EM2.
        prev_shifted_points % Set of time shifted points from last iteration of EM1.
        prev_shifted_weights % Set of weights from last iteration of EM1.
        template_decimated % Whether the template points have been downsampled with sphere decimation.
        id % Identifier for the feature.
        flow % Latest estimate of the flow of the feature.
        debug_fig % Figure handle to plot debug information.
        image_size % Size of the image.
        pos_history % History of past positions for this feature.
        R_init % Rotation of the camera when the feature was initialized. Used for em2_affine_imu.
        cinfo % Struct containing all of the intrinsics. Defined in the ROS message: http://docs.ros.org/api/sensor_msgs/html/msg/CameraInfo.html.
        A % Most recent affine transform between the events received and the template events.
        scale % Most recent scale between the events received and the template events.
        params % General parameters specified in get_params().
    end
    
    methods
        function [ef] = EventTracker(...
                image_size, ...
                feature_pos, ...
                id, ...
                R_init, ...
                cinfo, ...
                params, ...
                debug_fig)
            % Initializes the class with the initial position, as well as intrinsics and other params.
            %
            % Inputs:
            %    image_size  - resolution of the image [rows, cols].
            %    feature_pos - 2x1 initial position of the feature.
            %    id          - ID of the feature.
            %    R_init      - 3x3 initial 3D rotation of the camera.
            %    cinfo       - struct containing intrinsics of the camera.
            %    params      - parameters defined in getParams().
            %    debug_fig   - figure handle for debug plots.
            if ~EventTracker.checkInImage(...
                    feature_pos, ...
                    (params.window_size-1)/2, ...
                    image_size)
                ef.isvalid = 0;
                return;
            end
            ef.image_size = image_size;
            ef.feature_pos = feature_pos;
            ef.pos_history = feature_pos;
            ef.prev_shifted_points = [];
            ef.prev_shifted_weights = [];
            ef.id = id;
            ef.R_init = R_init;
            ef.cinfo = cinfo;
            ef.num_frames = 0;
            ef.isvalid = 1;
            ef.flow = [0; 0];
            ef.debug_fig = debug_fig;
            ef.template_decimated = false;
            ef.A = eye(2);
            ef.scale = 1;
            ef.template_points = []; 
            ef.params = params;
        end
        
        function [ef, feature_pos, isvalid, error] = stepTracker(...
                ef, ...
                events, ...
                R_curr)
            % stepTracker Tracks events across a time window.
            %
            % Inputs:
            %    events - 4xN matrix of (x,y,t,p) events.
            %    R_curr - 3x3 3D rotation of the camera.
            error = TrackingErrors.None;
            isvalid = false;
            feature_pos = ef.feature_pos;
            if ~ef.isvalid
                error = TrackingErrors.NotValid;
                return
            end            
                      
            ef.isvalid = false;            
            dt = events(3,end)-events(3, 1);
                        
            icpreason = error;
            
            % Estimate flow.
            if ef.params.use_em1_template && ~isempty(ef.prev_shifted_points)
                [ef.flow, shifted_points] = em1_flow_with_prior(...
                    events, ...
                    ef.feature_pos, ...
                    ef.prev_shifted_points, ...
                    ef.prev_shifted_weights, ...
                    ef.flow, ...
                    ef.params, ...
                    ef.debug_fig);
            else
                [ef.flow, shifted_points] = em1_flow(...
                    events, ...
                    ef.feature_pos, ...
                    ef.flow, ...
                    ef.params, ...
                    ef.debug_fig);
            end
            % Save shifted points for next EM1 iteration.
            if ef.params.use_em1_template
                [ ef.prev_shifted_points, ef.prev_shifted_weights ] = ...
                    decimate_points(shifted_points, 1.0);
                ef.prev_shifted_weights = ef.prev_shifted_weights.^2;
            end
            
            if any(isnan(ef.flow(:)))
                error = TrackingErrors.Flow;
                return
            end

            % Align the shifted points with the template to correct for
            % errors in the flow.
            if ~isempty(ef.template_points)
                if ef.params.use_imu
                    dR = ef.R_init * R_curr';
                    [new_point, shifted_points, ef.scale, icpreason] = em2_affine_imu(...
                        events, ...
                        ef.feature_pos, ...
                        ef.template_points, ...
                        ef.template_weights, ...
                        ef.scale, ...
                        ef.flow, ...
                        dR, ...
                        ef.cinfo, ...
                        ef.params, ...
                        ef.debug_fig);
                else
                    [new_point, shifted_points, ef.A, icpreason] = em2_affine(...
                        events, ...
                        ef.feature_pos, ...
                        ef.template_points, ...
                        ef.template_weights, ...
                        ef.flow, ...
                        ef.params, ...
                        ef.debug_fig, ...
                        ef.A);
                end
                % Due to line 164 (ef.feature_pos=...), this implicitly
                % sets ef.feature_pos = new_point.
                ef.flow = ef.flow + (new_point-ef.feature_pos)/dt;
            else
                new_point = ef.feature_pos;
            end
            
            if any(isnan(ef.flow)) || any(isnan(new_point))
                error = icpreason;
                return
            end       
            
            % Generate the template for EM2.
            if ef.params.use_initial_template
                % Generate a template from the first ef.params.num_init
                % iterations.
                if ...
                        ef.num_frames < ef.params.num_init && ...
                        size(ef.template_points, 2)+size(shifted_points, 2) < ef.params.max_template_points && ...
                        ~ef.template_decimated
                    if ...
                            size(ef.template_points, 2) < ef.params.max_template_points ...
                            && ef.num_frames < ef.params.num_init
                        ef.template_points = [ef.template_points shifted_points];
                    end
                elseif ~ef.template_decimated
                    [ ef.template_points, ef.template_weights] = decimate_points([ef.template_points shifted_points], 1);
                    if size(ef.template_points, 1) > ef.params.window_size^2/2
                        error = TrackingErrors.NotValid;
                        keyboard;
                        return
                    end
                    % Square the weights to weigh pixels with more events
                    % higher.
                    ef.template_weights = ef.template_weights.^2;
                    ef.template_decimated = true;
                end
            else
                % Use the current shifted events for the next iteration.
                ef.template_points = shifted_points;
                ef.R_init = R_curr;
            end
                        
            % Update the feature position.
            ef.feature_pos = ef.feature_pos + ef.flow*dt;
            
            feature_pos = ef.feature_pos;
            
            if ~EventTracker.checkInImage(ef.feature_pos, (ef.params.window_size-1)/2, ef.image_size)
                error = TrackingErrors.OutOfImg;
                return
            end
            
            ef.pos_history = [ef.pos_history feature_pos];
            ef.isvalid = true;
            ef.num_frames = ef.num_frames + 1;
            isvalid = true;
        end
    end
    
    methods (Static=true, Access=private)
        function [valid] = checkInImage(point, window_size, image_size)
            % checkInImage checks if a point with a given window size is
            % inside an image with size image_size.
            %
            % Inputs:
            %    point       - 2x1 pixel position.
            %    window_size - size of the square spatial window.
            %    image_size  - size of the image [rows, cols].
            valid = 1;
            if (point(1) < window_size || ...
                    point(2) < window_size ||...
                    point(1) > image_size(2) - window_size || ...
                    point(2) > image_size(1) - window_size)
                valid = 0;
            end
        end
        
    end
    
end


