function [params] = get_params()
%GET_PARAMS sets a number of parameters for ebpda_main.
% You must edit this file to change the parameters.
%
% Syntax: [params] = GET_PARAMS()
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

%% Main parameters to set
% Path to load the data .mat file.
params.data_path = '../data/indoor_flying1_data.mat';
% Path to load the undistort_map .mat file.
params.undistort_map_path = '../data/left_undistort_map.mat';
% Time in data to start in seconds.
params.t_start = 5;
% Time to stop in seconds. Set to -1 to finish at the end of the bag.
params.t_end = -1;
% Number of features to track.
params.num_features = 30;
% Set to true to default to def_int_time every time. Can be used to have a
% fixed number of events in each window instead of an adaptive window size.
params.constant_time = false;
% If true, will plot each individual EM step for visualization. If things
% aren't working, turn this on, you should see points aligning nicely.
params.debug = false;
% Initial temporal window size in seconds. Dependent on optical flow in the 
% scene, needs to be larger than 3 pixels of motion roughly.
params.def_int_time = 1;
% Set to true to use the IMU rotations for RANSAC, as in EVIO.
params.do_ransac = false;
% If true, tracking will be done inside a parfor loop.
params.do_parallel = true;
% Do tracking.
params.do_tracking = true;
% Set to true for no plots at all
params.headless = false;
% Set to true to use the IMU rotations for the affine step (as in EVIO).
params.use_imu = true;
% Size of feature window.
params.window_size = 31;

%% Other params for fine tuning. Usually don't need to touch these.
% Maximum number of events allowed in a window, upper limit for when the
% lifetime estimate becomes too high.
params.max_events_per_window = 30000;
% Sets integration time time for flow to move multiplier pixels. Minimum is
% around 2, 2.5 to be safe.
params.integration_multiplier = 3;
% If true, will use the events from the first iteration of a tracker in
% EM2. Otherwise, will use the events from the previous iteration.
params.use_initial_template = false;
% If true, will use the events from the previous iteration of a tracker as
% a template for EM1 (as in EVIO). Otherwise, the full optimization will be
% performed from scratch. In MATLAB, the full optimization is actually
% faster as it requires fewer iterations to converge, and for loops are
% much worse than the other computations.
params.use_em1_template = true;
% Number of initial point sets to concatenate for canonical set of events.
params.num_init = 3;
% Minimum number of events in the feature window required to perform EM.
params.min_events_for_em = 50;
% Maximum number of points in the template for EM2 before sphere decimation
% is applied.
params.max_template_points = 800;
% Number of iterations for RANSAC.
params.n_ransac_iters = 70;
% Threshold for inliers in RANSAC.
params.ransac_thresh = 4e-4;
% Minimum number of tracked features before adding new features.
params.min_features = round(params.num_features*.8);

%% Params for EM. Usually don't need to touch these.
% Maximum iterations of EM
em1_params.max_iters = 50;
% Covariance of target point Gaussians
em1_params.sigma = 2;
% Min change in error before optimization stops
em1_params.min_err = 1.0;
% Max distance for rangesearch for em1 beyond which points are considered
% outliers
em1_params.max_distance = 1.0;
params.em1_params = em1_params;

% Maximum iterations of ICP
em2_params.max_iters = 50;
% Covariance of target point Gaussians
em2_params.sigma = 2;
% Min change in error before optimization stops.
em2_params.min_err = 1e-3;
% Max distance for rangesearch for em2 beyond which points are considered
% outliers.
em2_params.max_distance = 1.0;
% Set to true to estimate a full affine transform between the current
% events and the template. Otherwise only translation will be estimated.
em2_params.doaffinewarp = true;
% Minimum change in scale between template and current set of events before
% the feature is dropped.
em2_params.min_scale = 0.2;
% Maximum outlier percentage for EM2 before a feature is dropped.
em2_params.max_outlier_thresh = 0.8;
params.em2_params = em2_params;

end

