function [detected_corners] = detect_corners(...
    image, ...
    curr_feature_pos, ...
    num_features, ...
    params)
%DETECT_CORNERS detects spatially separated Harris corners in an image.
%  All Harris corners are detected, and then k corners are returned that
%  are spatially far apart from the current set of features, such that the
%  total number of features is equal to num_features.
%
%  Syntax: [detected_corners] = DETECT_CORNERS(...
%             image, ...
%             curr_feature_pos, ...
%             num_features, ...
%             params)
%
%  Inputs:
%    image            - RxC grayscale image.
%    curr_feature_pos - 2xN set of current feature points.
%    num_features     - desired total number of features.
%    params           - parameters defined in get_params().
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

% Detect all Harris corners.
corners = detectHarrisFeatures(image, 'ROI', ...
    [params.window_size, ...
    params.window_size, ...
    size(image, 2) - params.window_size, ...
    size(image, 1) - params.window_size]);

% Append current features with essentially infinite metric to ensure that
% they are kept.
if ~isempty(curr_feature_pos)
    curr_corners = cornerPoints(curr_feature_pos, ...
        'Metric', ones(size(curr_feature_pos, 1), 1) * realmax);
    all_corners = [corners; curr_corners];
else
    all_corners = corners;
end

% Find a subset of corners that are uniformly distributed in the image.
uniform_points = selectUniform(all_corners, num_features, size(image));
% Find the new corners.
detected_corners = uniform_points.Location(~isinf(uniform_points.Metric), :)';
end