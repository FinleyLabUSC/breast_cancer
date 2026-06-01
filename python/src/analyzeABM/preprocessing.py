"""
TK 06.01.2026
"""
from analyzeABM.constants import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import alphashape
import shapely


def plot_alphashape(data, r=10):
    """
    Plots an alphashape given spatial data.
    :param data: DataFrame containing spatial and celltype data for a metastasis
    :param r: radius parameter for alphashape construction
    :return: alphashape object fit to data and approximate density of the metastasis
    """
    fig, ax = plt.subplots()
    data.plot.scatter("x", "y", ax=ax, c='k')
    alpha_shape = alphashape.alphashape(data[["x", "y"]], 1 / r)
    if type(alpha_shape) == shapely.geometry.multipolygon.MultiPolygon:
        polys = list(alpha_shape.geoms)
        for poly in polys:
            x, y = poly.exterior.xy
            ax.plot(x, y)
        ax.set_aspect('equal')
        print("Area = ", alpha_shape.area, " um^2")
        print("Est. density = ", len(data) / alpha_shape.area * 100, " cells/100 um^2")
    else:
        x, y = alpha_shape.exterior.xy
        print("Area = ", alpha_shape.area, " um^2")
        print("Est. density = ", len(data) / alpha_shape.area * 100, " cells/100 um^2")
        ax.plot(x, y)
        ax.set_aspect('equal')
    return alpha_shape, len(data) / alpha_shape.area * 100


def transform_spatial(data, r=10):
    """
    Transforms the spatial data so the density is approximately what would result from hexagonal packing of cells (i.e., no overlap).
    :param data: DataFrame containing spatial and celltype data for a metastasis
    :param r: radius parameter for alphashape construction
    :return: DataFrame containing transformed spatial and celltype data, the transformed alphashape object, and the scale factor
    """
    # Radii order is M0 / M1 / M2 / Cancer / CD4 / Treg / CD8 / NK / MDSC / Lymphoid / Myeloid / Stromal
    radii = {0:10, 1:10, 2:10, 3:10, 4:5, 5:5, 6:5, 8:7, 10:5, 11:5, 12:10, 13:10} # in um

    alpha_shape = plot_alphashape(data, r)  # Get the alpha shape & plot it
    density = len(data) / alpha_shape.area  # Get the approximate density of the alpha shape in cells / um^2
    data["area"] = data["state"].apply(lambda x: radii[x] ** 2 * np.pi)
    target_density = 1 / np.mean(data["area"]) * np.pi / np.sqrt(
        12)  # Get the avg. area of a cell & transform into density. Then multiply by factor to account for circle packing
    print("\nTarget density = ", target_density * 100, " cells/100 um^2")
    scale = np.sqrt(density / target_density)  # The change-in-coordinates needed is the sqrt of the ratio as we're in 2D
    print("Scale factor = ", scale, "\n")
    new_data = data
    new_data["x"] = data["x"] * scale
    new_data["y"] = data["y"] * scale
    new_alphashape = plot_alphashape(new_data, np.ceil(r * np.sqrt(scale)))  # Approximate what the new r should be; confirm using plots
    return new_data, new_alphashape, scale
