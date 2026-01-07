#!/usr/bin/env python3
# Modified by Raul Mur-Artal
# Automatically compute the optimal scale factor for monocular VO/SLAM.

# Software License Agreement (BSD License)
#
# Copyright (c) 2013, Juergen Sturm, TUM
# All rights reserved.

"""
This script computes the absolute trajectory error from the ground truth
trajectory and the estimated trajectory.
"""

import sys
import argparse
import numpy as np
import associate


def align(model: np.ndarray, data: np.ndarray):
    """Align two trajectories using the method of Horn (closed-form).

    Input:
    model -- first trajectory (3xn)
    data  -- second trajectory (3xn)

    Output:
    rot            -- rotation matrix (3x3)
    transGT        -- translation vector with scale correction (3x1)
    trans_errorGT  -- translational error per point with scale correction (n,)
    trans          -- translation vector without scale correction (3x1)
    trans_error    -- translational error per point without scale correction (n,)
    s              -- optimal scale factor (float)
    """
    np.set_printoptions(precision=3, suppress=True)

    model = np.asarray(model, dtype=float)
    data = np.asarray(data, dtype=float)

    model_mean = model.mean(axis=1, keepdims=True)
    data_mean = data.mean(axis=1, keepdims=True)

    model_zerocentered = model - model_mean
    data_zerocentered = data - data_mean

    W = np.zeros((3, 3), dtype=float)
    for i in range(model.shape[1]):
        W += np.outer(model_zerocentered[:, i], data_zerocentered[:, i])

    # SVD on W^T
    U, d, Vh = np.linalg.svd(W.T)

    S = np.eye(3, dtype=float)
    if np.linalg.det(U) * np.linalg.det(Vh) < 0:
        S[2, 2] = -1.0

    rot = U @ S @ Vh

    rotmodel = rot @ model_zerocentered
    dots = 0.0
    norms = 0.0
    for i in range(data_zerocentered.shape[1]):
        dots += float(data_zerocentered[:, i].T @ rotmodel[:, i])
        normi = np.linalg.norm(model_zerocentered[:, i])
        norms += float(normi * normi)

    s = float(dots / norms) if norms > 0 else 1.0

    transGT = data_mean - s * (rot @ model_mean)
    trans = data_mean - (rot @ model_mean)

    model_alignedGT = s * (rot @ model) + transGT
    model_aligned = (rot @ model) + trans

    alignment_errorGT = model_alignedGT - data
    alignment_error = model_aligned - data

    trans_errorGT = np.sqrt(np.sum(alignment_errorGT * alignment_errorGT, axis=0))
    trans_error = np.sqrt(np.sum(alignment_error * alignment_error, axis=0))

    return rot, transGT, trans_errorGT, trans, trans_error, s


def plot_traj(ax, stamps, traj, style, color, label):
    """
    Plot a trajectory using matplotlib.

    Input:
    ax     -- the plot
    stamps -- time stamps (list)
    traj   -- trajectory (nx3) or list of [x,y,z]
    style  -- line style
    color  -- line color
    label  -- plot legend
    """
    stamps = list(stamps)
    stamps.sort()

    interval = np.median([s - t for s, t in zip(stamps[1:], stamps[:-1])])
    x, y = [], []
    last = stamps[0]

    for i in range(len(stamps)):
        if stamps[i] - last < 2 * interval:
            x.append(traj[i][0])
            y.append(traj[i][1])
        elif len(x) > 0:
            ax.plot(x, y, style, color=color, label=label)
            label = ""
            x, y = [], []
        last = stamps[i]

    if len(x) > 0:
        ax.plot(x, y, style, color=color, label=label)


def rmse(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).reshape(-1)
    return float(np.sqrt(np.dot(x, x) / len(x)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This script computes the absolute trajectory error from the ground truth trajectory
    and the estimated trajectory.
    """)
    parser.add_argument("first_file", help="ground truth trajectory (format: timestamp tx ty tz qx qy qz qw)")
    parser.add_argument("second_file", help="estimated trajectory (format: timestamp tx ty tz qx qy qz qw)")
    parser.add_argument("--offset", default=0.0, type=float,
                        help="time offset added to the timestamps of the second file (default: 0.0)")
    parser.add_argument("--scale", default=1.0, type=float,
                        help="scaling factor for the second trajectory (default: 1.0)")
    parser.add_argument("--max_difference", default=20000000.0, type=float,
                        help="maximally allowed time difference for matching entries (default: 20000000 ns)")
    parser.add_argument("--save", help="save aligned second trajectory to disk (format: stamp2 x2 y2 z2)")
    parser.add_argument("--save_associations",
                        help="save associated first and aligned second trajectory to disk "
                             "(format: stamp1 x1 y1 z1 stamp2 x2 y2 z2)")
    parser.add_argument("--plot", help="plot the first and the aligned second trajectory to an image (format: png/pdf)")
    parser.add_argument("--verbose", action="store_true",
                        help="print all evaluation data (otherwise print CSV: rmse_noscale, scale, rmse_scale)")
    parser.add_argument("--verbose2", action="store_true",
                        help="print scale error and RMSE with and without scale correction")
    args = parser.parse_args()

    first_list = associate.read_file_list(args.first_file, False)
    second_list = associate.read_file_list(args.second_file, False)

    matches = associate.associate(first_list, second_list, float(args.offset), float(args.max_difference))
    if len(matches) < 2:
        sys.exit("Couldn't find matching timestamp pairs between groundtruth and estimated trajectory! "
                 "Did you choose the correct sequence?")

    first_xyz = np.array([[float(v) for v in first_list[a][0:3]] for a, b in matches], dtype=float).T
    second_xyz = np.array([[float(v) * float(args.scale) for v in second_list[b][0:3]] for a, b in matches],
                          dtype=float).T

    sorted_second_list = sorted(second_list.items())
    second_xyz_full = np.array(
        [[float(v) * float(args.scale) for v in sorted_second_list[i][1][0:3]] for i in range(len(sorted_second_list))],
        dtype=float
    ).T

    rot, transGT, trans_errorGT, trans, trans_error, scale = align(second_xyz, first_xyz)

    # scale-corrected alignment (recommended)
    second_xyz_aligned = scale * (rot @ second_xyz) + transGT
    second_xyz_full_aligned = scale * (rot @ second_xyz_full) + transGT

    # alignment without scale correction (for comparison)
    second_xyz_notscaled = (rot @ second_xyz) + trans
    second_xyz_notscaled_full = (rot @ second_xyz_full) + trans

    first_stamps = sorted(list(first_list.keys()))
    first_xyz_full = np.array([[float(v) for v in first_list[t][0:3]] for t in first_stamps], dtype=float).T

    second_stamps = sorted(list(second_list.keys()))
    second_xyz_full_raw = np.array([[float(v) * float(args.scale) for v in second_list[t][0:3]] for t in second_stamps],
                                   dtype=float).T
    second_xyz_full_raw_aligned = scale * (rot @ second_xyz_full_raw) + transGT

    if args.verbose:
        print(f"compared_pose_pairs {len(trans_error)} pairs")
        print(f"absolute_translational_error.rmse {rmse(trans_error):.6f} m")
        print(f"absolute_translational_error.mean {float(np.mean(trans_error)):.6f} m")
        print(f"absolute_translational_error.median {float(np.median(trans_error)):.6f} m")
        print(f"absolute_translational_error.std {float(np.std(trans_error)):.6f} m")
        print(f"absolute_translational_error.min {float(np.min(trans_error)):.6f} m")
        print(f"absolute_translational_error.max {float(np.max(trans_error)):.6f} m")
        print(f"max idx: {int(np.argmax(trans_error))}")
    else:
        # CSV: rmse_without_scale, optimal_scale, rmse_with_scale
        print(f"{rmse(trans_error):.6f},{scale:.6f},{rmse(trans_errorGT):.6f}")

    if args.verbose2:
        print(f"compared_pose_pairs {len(trans_error)} pairs")
        print(f"absolute_translational_error.rmse {rmse(trans_error):.6f} m")
        print(f"absolute_translational_errorGT.rmse {rmse(trans_errorGT):.6f} m")

    if args.save_associations:
        with open(args.save_associations, "w", encoding="utf-8") as f:
            lines = []
            for (a, b), (x1, y1, z1), (x2, y2, z2) in zip(
                matches, first_xyz.T, second_xyz_aligned.T
            ):
                lines.append(f"{float(a):.6f} {x1:.6f} {y1:.6f} {z1:.6f} {float(b):.6f} {x2:.6f} {y2:.6f} {z2:.6f}")
            f.write("\n".join(lines))

    if args.save:
        with open(args.save, "w", encoding="utf-8") as f:
            lines = []
            for stamp, xyz in zip(second_stamps, second_xyz_notscaled_full.T):
                lines.append(f"{float(stamp):.6f} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}")
            f.write("\n".join(lines))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)

        plot_traj(ax, first_stamps, first_xyz_full.T, "-", "black", "ground truth")
        plot_traj(ax, second_stamps, second_xyz_full_raw_aligned.T, "-", "blue", "estimated")

        label = "difference"
        for (a, b), (x1, y1, z1), (x2, y2, z2) in zip(matches, first_xyz.T, second_xyz_aligned.T):
            ax.plot([x1, x2], [y1, y2], "-", color="red", label=label)
            label = ""

        ax.legend()
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        plt.axis("equal")

        # 根据输出文件后缀选择格式（默认用后缀；没后缀就 png）
        out = args.plot
        fmt = "png"
        if "." in out:
            fmt = out.rsplit(".", 1)[-1].lower()
        plt.savefig(out, format=fmt)
