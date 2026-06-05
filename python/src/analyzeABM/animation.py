"""
TK 06.02.2026
"""
from analyzeABM.constants import *
import os
import numpy as np
import matplotlib.pyplot as plt
import shutil
import cv2
import pandas as pd
import multiprocessing
from functools import partial


def animate_tumor(csv_list, frames_dir="temp", fig_scale=1.5, delete_frames=True, vid_name="tumor_mov.avi"):
    """
    Automatically generate a video from a list of csv files that correspond to frames of spatial data. Colors according to the standard cell color cmap.
    :param csv_list: List of csv files containing spatial data
    :param frames_dir: Directory to store frames in (will be stored to "frames/{frames_dir}")
    :param fig_scale: How much larger than extent of metastasis at t = 0 to plot
    :param delete_frames: Whether to delete the directory storing the frames at the end of the function (recommended to save space)
    :param vid_name: Name of the video; must specify that the video is an AVI file
    :return: None
    """
    # Make frames directory
    path = f"frames/{frames_dir}"
    os.mkdir(path)

    # Get list of DataFrames
    frames = []
    for i, csv in enumerate(csv_list):
        frames.append([i, pd.read_csv(csv, usecols=[1, 2, 3, 4, 5], names=["cell_type", "x", "y", "radius", "state"])])

    # Determine xmin, xmax, ymin, and ymax
    fig_scale = 1 - fig_scale
    last_frame = frames[0][1]
    xmin = min(last_frame["x"])
    xmax = max(last_frame["x"])
    ymin = min(last_frame["y"])
    ymax = max(last_frame["y"])
    xspan = xmax - xmin
    yspan = ymax - ymin
    if xspan > yspan:
        diff = xspan - yspan
        xmax += fig_scale * xspan
        xmin -= fig_scale * xspan
        ymax += (diff / 2 + fig_scale * xspan)
        ymin -= (diff / 2 + fig_scale * xspan)
    else:
        diff = yspan - xspan
        ymax += fig_scale * yspan
        ymin -= fig_scale * yspan
        xmax += (diff / 2 + fig_scale * yspan)
        xmin -= (diff / 2 + fig_scale * yspan)

    # Generate frames
    frames.sort()
    frame_func = partial(plot_frame, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, frames_dir=frames_dir)
    pool = multiprocessing.Pool()
    pool.map(frame_func, enumerate(frames))

    # Generate video
    files = os.listdir(path)
    files.sort()
    imgnames = [os.path.join(path, file) for file in files]
    frame = cv2.imread(imgnames[0])
    h, w, l = frame.shape

    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    video = cv2.VideoWriter(vid_name, fourcc, 10, (w, h))  # 15 fps
    for fname in imgnames:
        video.write(cv2.imread(fname))
    video.release()

    if delete_frames:
        shutil.rmtree(path, ignore_errors=True)


def plot_frame(args, xmin=-1000, xmax=1000, ymin=-1000, ymax=1000, frames_dir="temp"):
    """
    Function that plots a frame used in making a video.
    :param args: Tuple of (frame_number, frame_df)
    :param xmin: xmin for plotting
    :param xmax: xmax for plotting
    :param ymin: ymin for plotting
    :param ymax: ymax for plotting
    :param frames_dir: Directory to store frames in
    :return: None
    """
    i, frame = args
    ts = frame[0]
    cells = frame[1]
    cells["state"] = cells["state"].replace({8: 7, 10: 8, 11: 9, 12: 10, 13: 11})

    scale = (6.4 / (max(xmax - xmin, ymax - ymin) * 2 + 1) * 300.) ** 2
    cells["adj_radius"] = scale * cells["radius"]

    fig, ax = plt.subplots(figsize=(10, 10 / 6.4 * 4.8), dpi=300)

    cells.plot.scatter("x", "y", c="state", s="adj_radius",
                       colormap=cell_cmap, vmin=-0.5, vmax=11.5,
                       xlim=[xmin, xmax], ylim=[ymin, ymax], ax=ax)

    f = plt.gcf()
    [pax, cax] = f.get_axes()
    pax.set(xlabel="x-coordinate", ylabel="y-coordinate", title=f"Metastasis at t = {ts:06.3f}")
    cax.set_ylabel('')
    cax.set_ylim(-0.5, num_cell_types - 0.5)

    cax.set_yticks(np.arange(0, num_cell_types), labels=cell_labels)
    f.savefig(f"frames/{frames_dir}/f_{i:05d}.png", dpi=300)
    plt.close()
