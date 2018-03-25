#!/usr/bin/env python
import numpy as np
import cv2
import argparse
import yaml
from scipy.io import savemat
"""
generate_undistort_map_equi.py generates a set of undistort maps for a given camera with 
equidistant distortion. Works for both monocular and stereo cameras.

Author: Alex Zihao Zhu, University of Pennsylvania
Email: alexzhu(at)seas.upenn.edu
Copyright 2018 University of Pennsylvania 
Alex Zihao Zhu, Nikolay Atanasov, Kostas Daniilidis

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS, CONTRIBUTORS, AND THE 
TRUSTEES OF THE UNIVERSITY OF PENNSYLVANIA "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT OWNER, CONTRIBUTORS OR THE TRUSTEES OF 
THE UNIVERSITY OF PENNSYLVANIA BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

def main():
    parser = argparse.ArgumentParser(
        description="Generates an undistort map for a given set of camera parameters.")
    parser.add_argument("--camchain", dest="camchain",
                        help="Camchain yaml file path with camera info.",
                        required=True)
    path = parser.parse_args().camchain
    
    with open(path, 'r') as stream:
        try:
            data = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            
    # If stereo, calculate the rectification parameters
    is_stereo = 'cam0' in data and 'cam1' in data

    intrinsics0 = data['cam0']['intrinsics']
    K0 = np.array([[intrinsics0[0], 0, intrinsics0[2]],
                   [0, intrinsics0[1], intrinsics0[3]],
                   [0, 0, 1]])
    d0 = np.array(data['cam0']['distortion_coeffs'])
    resolution = data['cam0']['resolution']

    if is_stereo:
        intrinsics1 = data['cam1']['intrinsics']
        K1 = np.array([[intrinsics1[0], 0, intrinsics1[2]],
                             [0, intrinsics1[1], intrinsics1[3]],
                             [0, 0, 1]])
        d1 = np.array(data['cam1']['distortion_coeffs'])
        H_12 = np.array(data['cam1']['T_cn_cnm1'])
        R = H_12[:3, :3]
        t = H_12[:3, 3]

        R0, R1, P0, P1,Q = cv2.fisheye.stereoRectify(K0,
                                                     d0,
                                                     K1,
                                                     d1,
                                                     tuple(resolution),
                                                     R,
                                                     t,
                                                     cv2.CALIB_ZERO_DISPARITY)

    rectifyMap0 = np.zeros((resolution[0],
                            resolution[1],
                            2))
    n_points = resolution[0] * resolution[1]
    points_distorted = np.zeros((n_points, 2))
    point_iter = 0
    for r in range(resolution[1]):
        for c in range(resolution[0]):
            points_distorted[point_iter, :] = [c, r]
            point_iter += 1

    # Oddly, cv2 expects points_distorted to be a Nx2x1 dimensional matrix.
    if is_stereo:
        points_undistorted0 = cv2.fisheye.undistortPoints(points_distorted[:, None],
                                                          K0,
                                                          d0,
                                                          R=R0,
                                                          P=P0)
    else:
        points_undistorted0 = cv2.fisheye.undistortPoints(points_distorted[:, None],
                                                          K0,
                                                          d0)
        
    points_undistorted0 = np.reshape(points_undistorted0,
                                     (resolution[1],
                                      resolution[0],
                                      2))
    if is_stereo:
        points_undistorted1 = cv2.fisheye.undistortPoints(points_distorted[:, None],
                                                          K1,
                                                          d1,
                                                          R=R1,
                                                          P=P1)
        points_undistorted1 = np.reshape(points_undistorted1,
                                         (resolution[1],
                                          resolution[0],
                                          2))

    left_undistort_dict = { 'undistort_map_x' : points_undistorted0[..., 0],
                            'undistort_map_y' : points_undistorted0[..., 1] }
    savemat('left_undistort_map', left_undistort_dict)
    if is_stereo:
        right_undistort_dict = { 'undistort_map_x' : points_undistorted1[..., 0],
                                 'undistort_map_y' : points_undistorted1[..., 1] }
        savemat('right_undistort_map', right_undistort_dict)
        
if __name__ == "__main__":
    main()

