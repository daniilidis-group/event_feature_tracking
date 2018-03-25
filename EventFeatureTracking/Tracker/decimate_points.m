function [ decimated_points, decimated_weights ] = decimate_points(points, radius)
% DECIMATE_POINTS Decimates points using sphere decimation.
% As in EM-ICP, the points are decimated by clustering points less than 
% radius apart with a heuristic.
%
% Syntax: [ decimated_points, decimated_weights ] = DECIMATE_POINTS(points, radius)
%
% See: S. Granger, and P. Xavier. 
%      "Multi-scale EM-ICP: A fast and robust approach for surface registration."
%      ECCV 2002.
%
% Inputs: points - 2xN matrix of 2D points
%         radius - max distance for clustering
%
% Outputs: decimated points  - 2xM, barycenter of each cluster
%          decimated_weights - Mx1, number of points in each cluster
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

num_points = size(points, 2);
decimated_points = zeros(2, num_points);
decimated_weights = zeros(num_points, 1);
num_dec_points = 0;
while size(points, 2)
    iter = randi(size(points, 2));
    converged = false;
    barycenter = points(:, iter);
    nearby_points = [];
    while ~converged
        [idx, ~] = rangesearch_bruteforce(barycenter, points, radius);
        barycenter = mean(points(:, idx), 2);
        converged = true;
        if all(size(nearby_points)==size(idx)) && ...
                all(nearby_points==new_nearby_points)
            converged = true;
        end
        nearby_points = idx;
    end
    decimated_points(:, num_dec_points+1) = barycenter;
        
    decimated_weights(num_dec_points + 1) = length(nearby_points);
    num_dec_points = num_dec_points + 1;
    points(:, nearby_points) = [];
end
decimated_points = decimated_points(:, 1:num_dec_points);
decimated_weights = decimated_weights(1:num_dec_points);
end

