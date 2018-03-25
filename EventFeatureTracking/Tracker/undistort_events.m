function [undistorted_events] = undistort_events(...
    events, ...
    undistort_map_x, ...
    undistort_map_y, ...
    image_size)
%UNDISTORT_EVENTS undistorts a set of events given the x and y undistortion maps.
%
%  Syntax: [undistorted_events] = UNDISTORT_EVENTS(...
%           events, ...
%           rectify_map_x, ...
%           rectify_map_y, ...
%           image_size)
%
%  Inputs:
%    events        - 4xN set of input (x,y,t,p) events.
%    undistort_map_x - RxC matrix where the distorted pixel (x,y) maps to
%      the x position of undistort_map_x(y-1, x-1)+1. Note the change from
%      zero index to one index.
%    undistort_map_y - RxC matrix where the distorted pixel (x,y) maps to
%      the y position of undistort_map_x(y-1, x-1)+1. Note the change from
%      zero index to one index.
%    image_size    - size of the image [rows, cols].
%
%  Outputs:
%    undistorted_events - 4xN set of undistorted events.
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

inds = sub2indc(events(1, :)+1, events(2, :)+1, size(undistort_map_x));
x_u = undistort_map_x(inds);
y_u = undistort_map_y(inds);
% Take only events that are in the original image.
valid_points = round(x_u) >= 0 & round(x_u) < image_size(2) & ...
    round(y_u) >= 0 & round(y_u) < image_size(1);
x_u = x_u(valid_points);
y_u = y_u(valid_points);
undistorted_events = [x_u; y_u; events(3:4, valid_points)];
end

