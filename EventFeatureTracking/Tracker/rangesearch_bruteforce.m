function [idx, dist]=rangesearch_bruteforce(query, points, max_dist)
%RANGESEARCH_BRUTEFORCE is a vectorized rangesearch between a single query point and a set of points.
% This is faster than building a KDTree if you just want to query a single 
% point in a relatively small set of points.
%
% Syntax: [idx, dist]=RANGESEARCH_BRUTEFORCE(query, points, max_dist)
%
% Inputs: 
%    query - 2x1 query point.
%    points - 2xN point set.
%    max_dist - max dist for a point to count as an 'inlier'.
%
% Outputs: 
%    idx - indices in points within max_dist of query.
%    dist - distance of each point in idx from query.
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

dist = sqrt(sum(bsxfun(@minus, points, query).^2, 1));
fidx=dist < max_dist*max_dist;
idx=find(fidx);
dist=dist(fidx);