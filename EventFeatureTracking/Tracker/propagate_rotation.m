function [R_new, imu_iter] = propagate_rotation(...
    imu_iter, ...
    imu_times, ...
    curr_time, ...
    dt, ...
    ang_vel, ...
    R)
%PROPAGATE_ROTATION updates a 3D rotation by integrating a set of angular velocities.
%
% Syntax: [R_new, imu_iter] = PROPAGATE_ROTATION(...
%           imu_iter, ...
%           imu_times, ...
%           curr_time, ...
%           dt, ...
%           ang_vel, ...
%           R)
%
% Inputs:
%    imu_iter  - first index of ang_vel to be used.
%    imu_times - Nx1, timestamps for each ang_vel.
%    curr_time - timestamp to integrate up to.
%    dt        - equal to curr_time - last_time.
%    ang_vel   - 3xN, angular velocities.
%    R         - 3x3, initial rotation matrix.
%
% Outputs:
%    R_new    - 3x3, updated rotation matrix.
%    imu_iter - last index of ang_vel used.
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

% Find last index of ang_vel to use.
n_imu = find(imu_times(imu_iter:end) > curr_time, 1);

avg_angvel = mean(ang_vel(:, imu_iter:imu_iter + n_imu), 2);
imu_iter = imu_iter + n_imu;

% Integrate the ang_vel to update the rotation using quaternions.
q = rotm2quat(R)';
q_dot = 0.5 * omegaMat(avg_angvel * dt);
q_prop = q + q_dot * q;
q_prop = q_prop / norm(q_prop);
R_new = quat2rotm(q_prop');

    function [omega_hat] = omegaMat(omega)
        omega_hat = [0, omega'; -omega, -crossMat(omega)];
    end

    function [cross_mat] = crossMat(vec)
        cross_mat = [...
            0, -vec(3), vec(2); ...
            vec(3), 0, -vec(1); ...
            -vec(2), vec(1), 0];
    end
end

