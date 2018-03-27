function [] = extract_ros_data(file_name, load_folder, save_folder, cam_name)
%EXTRACT_ROS_DATA reads a ROS bag, and saves the events, IMU messages and
%camera info in a mat file.
%
% Syntax: extract_ros_data(file_name, load_folder, save_folder)
%
% You will need the matlab_rosbag package:
% https://github.com/bcharrow/matlab_rosbag
%
% Inputs:
%    file_name   - name of the bag file (no extension).
%    load_folder - folder of the bag. No trailing slash. Set to . by
%                  default.
%    save_folder - folder the mat is saved in. Set to load_folder by
%                  default.
%    cam_name    - OPTIONAL name of the event camera topic if not davis or
%                  dvs.
%
% See also GENERATE_UNDISTORT_MAP
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

if nargin==1
    load_folder = '.';
end

if nargin < 3
    save_folder = load_folder;
end

if nargin == 4
    fprintf(['Warning: You have set the cam_name argument to ''%s''. ', ...
        'If this is intentional, please ignore this message.\n'], cam_name)
end

bag_name = strcat(load_folder,'/', file_name, '.bag');

fprintf('Reading ROS bag %s. This may take some time.\n', bag_name)

bag = ros.Bag.load(bag_name);

event_topics = bag.topics('.*events');

if isempty(event_topics)
    fprintf('Error, no topic with events found in bag %s\n', bag_name);
    return
end

if contains(event_topics{1}, 'davis')
    cam_name = 'davis';
elseif contains(event_topics{1}, 'dvs')
    cam_name = 'dvs';
elseif nargin < 4
    fprintf(['Could not find the default camera names ''davis'' or ''dvs'' in the bag.', ...
        'Please enter the camera name as a fourth argument to this function.\n']);
    return
end

if length(event_topics) == 1    
    event_topic = ['/' cam_name '/events'];
    cinfo_topic = ['/' cam_name '/camera_info'];
    imu_topic = ['/' cam_name '/imu'];
elseif length(event_topics) == 2
    if contains(event_topics{1}, 'left') || contains(event_topics{2}, 'left')
        event_topic = ['/' cam_name '/left/events'];
        cinfo_topic = ['/' cam_name '/left/camera_info'];
        imu_topic = ['/' cam_name '/left/imu'];
    else
        fprintf('No left event topic found despite two event topics found.\n')
    end
else
    fprintf('Can''t find monocular or stereo events.\n')
end

[event_msgs, ~] = bag.readAll(event_topic);
[imu_msgs, ~] = bag.readAll(imu_topic);
[cinfo_msgs, ~] = bag.readAll(cinfo_topic);

num_event = 0;

for i=1:length(event_msgs)
    num_event = num_event + size(event_msgs{i}.events, 2);
end

cinfo = cinfo_msgs{1};

events = zeros(4, num_event);

iter = 1;

for i=1:length(event_msgs)
    events(:, iter:iter+size(event_msgs{i}.events,2)-1) = event_msgs{i}.events;
    iter = iter+size(event_msgs{i}.events,2);
end

num_imu = length(imu_msgs);

lin_acc = zeros(3, num_imu);
ang_vel = zeros(3, num_imu);
imu_time = zeros(num_imu, 1);

for i=1:num_imu
    lin_acc(:,i) = imu_msgs{i}.linear_acceleration;
    ang_vel(:,i) = imu_msgs{i}.angular_velocity;
    imu_time(i) = imu_msgs{i}.header.stamp.time;
end

output_name = strcat(save_folder,'/',file_name);

fprintf('Saving data to %s. This may take some time.\n', output_name)

save(output_name,'file_name','events',...
    'lin_acc','ang_vel','imu_time', 'cinfo', '-v7.3');
end