function [ best_inliers ] = two_point_ransac( pointsA, pointsB, R, params )
%TWO_POINT_RANSAC performs outlier rejection on two sets of points given camera rotation.
% Given two sets of corresponding points between two iterations in time and
% the rotation of the camera between the iterations, finds the largest set
% of points that agrees on a translation direction using RANSAC.
%
% Syntax: [ best_inliers ] = TWO_POINT_RANSAC( pointsA, pointsB, R, params )
%
% Inputs: 
%    pointsA - 2xN or matrix of 2D undistorted points. Can be 3xN if
%              the last row is all 1s (homogeneous representation).
%    pointsB - 2xN or matrix of 2D undistorted points, as above.
%              Note: it is assumed that A(:, i) corresponds to B(:, i)
%    R       - 3x3 matrix corresponding to the rotation of the camera
%              between the iterations A and B.
%    params  - parameters defined in get_params().
% Outputs: 
%    best_inliers - indices into A and B for the largest inlier set.
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

num_pts = size(pointsA, 2);

% No point doing RANSAC on too few points.
if num_pts < 5
    best_inliers = ones(num_pts, 1);
    return
end

% Convert to homogeneous form if necessary.
if size(pointsA, 1) == 2
    pointsA = [pointsA; ones(1, num_pts)];
    pointsB = [pointsB; ones(1, num_pts)];
end

num_iters = 70;
thresh = params.ransac_thresh;

best_inliers = [];
most_inliers = -1;

for i=1:num_iters
    
    inds = randperm(num_pts, 2);
    t = solveTwoPoint(pointsA(:, inds), pointsB(:, inds), R);
    if isempty(t)
        continue;
    end
    
    errs = computeSampsonError(pointsA, pointsB, R, t);
    inliers = errs < thresh;
    num_inliers = sum(inliers);
    if num_inliers > most_inliers
        best_inliers = inliers;
        most_inliers = num_inliers;
    end
end

    function [err] = computeSampsonError(p1, p2, R, t)
        E = skewsymm(t) * R;
        Ex1 = E*p1;
        Ex2 = E'*p2;
        err = sum(p2.*Ex1, 1).^2;
        err = err ./ (Ex1(1, :).^2+Ex1(2,:).^2+Ex2(1,:).^2+Ex2(2,:).^2);
    end

    function [t] = solveTwoPoint(pointsA, pointsB, R)
        M = zeros(2, 3);
        M(1, :) = (R*pointsA(:, 1))'*skewsymm(pointsB(:, 1));
        M(2, :) = (R*pointsA(:, 2))'*skewsymm(pointsB(:, 2));
        t = null(M);
        if size(t, 2) > 1
            t = [];
        end
    end

    function [matskew] = skewsymm(mat)
        matskew = [0 -mat(3) mat(2);...
            mat(3) 0 -mat(1);...
            -mat(2) mat(1) 0];
    end


end

