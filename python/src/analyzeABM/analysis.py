"""
TK 06.01.2026
"""
from analyzeABM.constants import *
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
from collections import defaultdict
import os


def extract_neighbors(input_var):
    """
    :param input_var: list containing path to data, met number, rep number, and hour
    :return: dataframe with celltype, cell UID, rep/met/hour info, and neighbors vector
    """
    # Extract data for this time point
    path, met, rep, hour = input_var
    data = pd.read_csv(path, usecols=[0, 2, 3, 5], names=["CellID", "x", "y", "state"])
    data["ctype"] = data["state"].map(state_to_name).astype("category")
    data["ctype"] = data["ctype"].cat.set_categories(list(state_to_name.values()))
    data["UID"] = [f"{met}_{hour}_{idx}" for idx in data["CellID"]]
    data["rep"] = rep
    data["met"] = met
    data["hour"] = hour

    # Calculate nearest neighbors, with a new instance each time
    # There are 15 cells in a window, but need to include self
    nbrs = NearestNeighbors(n_neighbors=14, algorithm="auto").fit(data[["x", "y"]])
    dist, indices = nbrs.kneighbors(data[["x", "y"]])

    nvects = []
    for i, idx in enumerate(indices):
        nvects.append(
            data.iloc[idx].ctype.value_counts().sort_index().to_numpy())  # This gives ordered list of neighbors
    data["neighbors"] = nvects

    return data


def retrieve_ts(folder_name, num_replicates=5, max_len=600):
    """
    Reads timeseries ABM output data into a DataFrame
    :param folder_name: base folder where ABM output is located
    :param num_replicates: number of replicates (must have run replicates using seeds from 0 to num_replicates - 1)
    :param max_len: max number of hours to extract
    :return: DataFrame with (nonspatial) timeseries data
    """
    files = []
    for i in range(num_replicates):
        file_name = folder_name + "/set_" + str(i) + "/populations_TS.csv"
        files.append(file_name)
    dataframe_dict = {}
    cols_to_get = ["m0", "m1", "m2", "c", "cd4_th", "cd4_treg", "cd8", "nk", "mdsc", "myeloid", "lymphoid", "stromal",
                   "radius"]

    # Create all the files
    for i, file_name in enumerate(files):
        df = pd.read_csv(file_name, header=0)

        # Append final state until max length
        last_time = max(df.index)
        max_index = last_time
        while max_index < max_len:
            df.loc[max_index + 1] = df.loc[last_time]
            max_index += 1

        if i == 0:
            for col in cols_to_get:
                dataframe_dict[col] = df[[col]].copy()
                dataframe_dict[col].rename({col: "0"}, axis=1, inplace=True)
        else:
            for col in cols_to_get:
                dataframe_dict[col][str(i)] = df[col].to_numpy()

    return dataframe_dict


def analyze_single_bdp(path_to_base, window=6):
    """
    Analyze a single replicate as a birth-death process.
    :param path_to_base: Path to ABM outputs for a specific replicate
    :param window: Number of hours to average over
    :return: DataFrame with cancer death data (instantaneous and averaged)
    """
    pop_ts = path_to_base + "/populations_TS.csv"
    deaths = path_to_base + "/cancer_death.csv"

    pop_df = pd.read_csv(pop_ts, header=0)
    cpop_list = pop_df["c"].to_list()
    dc = []
    for i in range(len(cpop_list) - 1):
        dc.append(cpop_list[i + 1] - cpop_list[i])

    death_df = pd.read_csv(deaths, header=0)
    death_df["tot_deaths"] = death_df["age"] + death_df["cd8"] + death_df["nk"]
    death_df["dc_dt"] = dc
    death_df["births"] = death_df["dc_dt"] + death_df["tot_deaths"]
    death_df["bdp"] = death_df["births"] - death_df["age"]

    # Rolling averages -- 6-hour window
    death_df["tot_deaths_a"] = death_df["tot_deaths"].rolling(window).mean()
    death_df["dc_dt_a"] = death_df["dc_dt"].rolling(window).mean()
    death_df["births_a"] = death_df["births"].rolling(window).mean()
    death_df["cd8_a"] = death_df["cd8"].rolling(window).mean()
    death_df["nk_a"] = death_df["nk"].rolling(window).mean()
    death_df["bdp_a"] = death_df["bdp"].rolling(window).mean()

    return death_df


def analyze_pop_bdp(path_to_base, num_replicates=5, max_len=600, window=6,
                    ax=None, err_type=("sd", 1), title=None):
    """
    Population-level equivalent of analyze_single_bdp(), with optional plotting.
    :param path_to_base: Path to ABM outputs for a specific metastasis
    :param num_replicates: Number of replicates (must have run replicates using seeds from 0 to num_replicates - 1)
    :param max_len: Max number of hours to analyze
    :param window: Number of hours to average over
    :param ax: Ax to plot on if you want to plot (if None, will not plot)
    :param err_type: What type of error to show if plotting
    :param title: Title of plot if plotting
    :return: Dict of DataFrames containing BDP data
    """
    base_paths = []
    for i in range(num_replicates):
        file_name = path_to_base + "/set_" + str(i)
        base_paths.append(file_name)
    dataframe_dict = {}

    # Create all the files
    for i, file_name in enumerate(base_paths):
        death_df = analyze_single_bdp(file_name, window)

        # Append final state until max length
        last_time = max(death_df.index)
        max_index = last_time
        while max_index < max_len:
            death_df.loc[max_index + 1] = death_df.loc[last_time]
            max_index += 1

        if i == 0:
            for col in death_df.columns.tolist():
                dataframe_dict[col] = death_df[[col]].copy()
                dataframe_dict[col].rename({col: "0"}, axis=1, inplace=True)
        else:
            for col in death_df.columns.tolist():
                dataframe_dict[col][str(i)] = death_df[col].to_numpy()

    # Return df if no axes to plot
    if not ax:
        return dataframe_dict

    # Else plot
    types_to_plot = ["tot_deaths_a", "dc_dt_a", "births_a", "cd8_a", "nk_a", "bdp_a"]
    melted = []

    for ctype in types_to_plot:
        melt = dataframe_dict[ctype].melt(var_name='Replicate', value_name=ctype, ignore_index=False)
        melt["Time"] = melt.index
        melted.append([ctype, melt])

    for entry in melted:
        sns.lineplot(data=entry[1].reset_index(), x="Time", y=entry[0], errorbar=err_type, ax=ax)

    if "radius" != types_to_plot:
        ax.set(xlabel="Time (hr)", ylabel=fr"$\Delta$ cells per time", title=title)

    return dataframe_dict


def retrieve_kill_ts(folder_name, num_replicates=5, max_len=600):
    """
    Retrieves the killing timeseries data from multiple replicates
    :param folder_name: Path to ABM outputs
    :param num_replicates: Number of replicates
    :param max_len: Max number of hours to analyze
    :return: DataFrame with killing data
    """
    files = []
    for i in range(num_replicates):
        file_name = folder_name + "/set_" + str(i) + "/cancer_death.csv"
        files.append(file_name)

    dataframe_dict = {}
    cols_to_get = ["age", "cd8", "nk"]

    # Create all the files
    for i, file_name in enumerate(files):
        df = pd.read_csv(file_name, header=0)

        # Append final state until max length
        last_time = max(df.index)
        max_index = last_time
        while max_index < max_len:
            df.loc[max_index + 1] = df.loc[last_time]
            max_index += 1

        if i == 0:
            for col in cols_to_get:
                dataframe_dict[col] = df[[col]].copy()
                dataframe_dict[col].rename({col: "0"}, axis=1, inplace=True)
        else:
            for col in cols_to_get:
                dataframe_dict[col][str(i)] = df[col].to_numpy()

    return dataframe_dict


def extract_lifespan_dist(path, ctype="CD8+"):
    """
    Extracts the distribution of lifespans of certain cell types.
    :param path: Path to ABM output (single simulation)
    :param ctype: What celltype to extract (CD8+ or NK only)
    :return: Returns dictionary of cell UIDs and corresponding lifespans in simulation
    """
    get_type = name_to_state[ctype]

    # Exclude cells at t = 0 and final timepoint
    duration = len(os.listdir(os.path.join(path, "cellLists")))
    data_at_zero = pd.read_csv(os.path.join(path, "cellLists", f"day_0", "cells.csv"), usecols=[0, 5],
                               names=["UID", "state"])
    data_at_end = pd.read_csv(os.path.join(path, "cellLists", f"day_{duration - 1}", "cells.csv"), usecols=[0, 5],
                              names=["UID", "state"])

    # Isolate cells of the correct type and get their UIDs to exclude
    keep_at_zero = data_at_zero[data_at_zero["state"] == get_type]
    uid_at_zero = list(keep_at_zero["UID"])
    keep_at_end = data_at_end[data_at_end["state"] == get_type]
    uid_at_end = list(keep_at_end["UID"])

    uid_to_lifespan = defaultdict(int)  # Dictionary of UID : lifespan pairs

    for i in range(len(os.listdir(os.path.join(path, "cellLists"))) - 1):
        fullpath = os.path.join(path, "cellLists", f"day_{i + 1}", "cells.csv")
        data = pd.read_csv(fullpath, usecols=[0, 5], names=["UID", "state"])
        keep = data[data["state"] == get_type]
        uids = list(keep["UID"])

        # If the UID is present, iterate by 1
        for uid in uids:
            uid_to_lifespan[uid] += 1

    for uid in uid_at_zero:
        uid_to_lifespan[uid] = 0
        del uid_to_lifespan[uid]

    for uid in uid_at_end:
        uid_to_lifespan[uid] = 0
        del uid_to_lifespan[uid]

    print(f"Deleted {len(uid_at_zero)} cells at t = 0 and {len(uid_at_end)} at t = {duration - 1}")
    return uid_to_lifespan


def gompertz(x, b):
    # Defines a Gompertzian distribution of observed lifespans
    # b is the exponential increase in hazard s.t. H(t) = ae^(bt)
    a = 1 / (14 * 24)  # The base death rate is fixed
    n = a / b
    return b * n * np.exp(n + b * x - n * np.exp(b * x))


def bi_gompertz(x, b1, b2, p):
    # Defines a weighted sum of 2 Gompertzian distributions of observed lifespans
    # b1 and b2 are the exponential increases in hazard
    # p is the proportion (from 0 to 1) of the distribution associated with b1
    a = 1 / (14 * 24)  # The base death rate is fixed
    n1 = a / b1
    n2 = a / b2
    g1 = b1 * n1 * np.exp(n1 + b1 * x - n1 * np.exp(b1 * x))
    g2 = b2 * n2 * np.exp(n2 + b2 * x - n2 * np.exp(b2 * x))
    return p * g1 + (1 - p) * g2


# Extract exhaustion trajectories for CD8+ cells
def extract_exhaustion_trajectories(path, ctype="CD8+", remove_end=False):
    """
    Extracts the exhaustion trajectories of certain cell types.

    path should be the base directory for the simulation.
    ctype can ONLY be "CD8+" or "NK"
    """
    get_type = name_to_state[ctype]

    # Exclude cells at t = 0 and final timepoint
    duration = len(os.listdir(os.path.join(path, "cellLists")))
    data_at_zero = pd.read_csv(os.path.join(path, "cellLists", f"day_0", "cells.csv"), usecols=[0, 5],
                               names=["UID", "state"])
    data_at_end = pd.read_csv(os.path.join(path, "cellLists", f"day_{duration - 1}", "cells.csv"), usecols=[0, 5],
                              names=["UID", "state"])

    # Isolate cells of the correct type and get their UIDs to exclude
    keep_at_zero = data_at_zero[data_at_zero["state"] == get_type]
    uid_at_zero = list(keep_at_zero["UID"])
    keep_at_end = data_at_end[data_at_end["state"] == get_type]
    uid_at_end = list(keep_at_end["UID"])

    uid_to_exh = defaultdict(list)  # Dictionary of UID : exhaustion-trajectory pairs

    for i in range(len(os.listdir(os.path.join(path, "cellLists"))) - 1):
        fullpath = os.path.join(path, "cellLists", f"day_{i + 1}", "cells.csv")
        data = pd.read_csv(fullpath, usecols=[0, 5, 7, 8, 9, 10, 11, 12],
                           names=["UID", "state", "PD1_tot", "PD1_avail", "speed", "p(kill)", "p(death)", "p(div)"])
        keep = data[data["state"] == get_type]
        uids = list(keep["UID"])
        keep = keep.set_index("UID")

        for uid in uids:
            store = keep.loc[uid].tolist()
            store.append(i)
            uid_to_exh[uid].append(store)  # Append cols 7-12 (can elim 0 & 5 later)

    for uid in uid_at_zero:
        uid_to_exh[uid] = []
        del uid_to_exh[uid]

    if remove_end:
        for uid in uid_at_end:
            uid_to_exh[uid] = []
            del uid_to_exh[uid]

    # Convert to DF
    for k, v in uid_to_exh.items():
        uid_to_exh[k] = pd.DataFrame(v,
                                     columns=["state", "PD1_tot", "PD1_avail", "speed", "p(kill)", "p(death)", "p(div)",
                                              "time"])

    return uid_to_exh


def align_on_exh(df_dict, feature="p(kill)", alignment="to_birth", max_time=600):
    """
    Aligns DF dict produced by extractor for a single exhaustion feature.

    feature must either be in ["PD1_tot", "PD1_avail", "speed", "p(kill)", "p(death)", "p(div)"]
    alignment is either "to_birth" or "to_time"
    """
    if feature not in ["PD1_tot", "PD1_avail", "speed", "p(kill)", "p(death)", "p(div)"]:
        raise ValueError(f"Feature {feature} cannot be used.")

    if alignment not in ["to_birth", "to_time"]:
        raise ValueError(f"Alignment {alignment} cannot be used.")

    data_list = []
    uid_list = []
    for k, v in df_dict.items():
        uid_list.append(k)
        if alignment == "to_time":
            start = min(v["time"])
            end = max(v["time"])
            data = np.concatenate((np.zeros(start),
                                   v[feature].to_numpy(),
                                   np.zeros(max_time - end)),
                                  axis=None)
        else:
            data = v[feature].to_numpy()

        data_list.append(data)

    if alignment == "to_birth":
        max_len = max([len(row) for row in data_list])
        data_list = np.array([np.pad(row, (0, max_len - len(row))) for row in data_list])
    else:
        data_list = np.array(data_list)

    return data_list, uid_list


def calc_neigh_traj(df, remove_start=True, remove_end=False):
    # Calculates MDSC, CD8+, and cancer neighbors by default

    trajectories = defaultdict(list)

    uid_at_zero = df[df.hour == 0].CellID.to_list()
    uid_at_end = df[df.hour == max(df.hour)].CellID.to_list()

    for hour in np.sort(df.hour.unique()):
        subset = df[df.hour == hour].set_index("CellID")
        for cell in subset.index.values:
            data = subset.loc[cell][["MDSC_neigh", "CD8+_neigh", "cancer_neigh"]].to_list()
            data.append(hour)
            trajectories[cell].append(data)

    if remove_start:
        for uid in uid_at_zero:
            trajectories[uid] = []
            del trajectories[uid]

    if remove_end:
        for uid in uid_at_end:
            trajectories[uid] = []
            del trajectories[uid]

    # Convert to DF
    for k, v in trajectories.items():
        trajectories[k] = pd.DataFrame(v, columns=["MDSC_neigh", "CD8+_neigh", "cancer_neigh", "time"])

    return trajectories


def align_on_neigh(df_dict, feature="MDSC_neigh", alignment="to_birth", max_time=600):
    """
    Aligns DF dict produced by extractor for a single exhaustion feature.

    feature must either be in ["MDSC_neigh", "CD8+_neigh", "cancer_neigh"]
    alignment is either "to_birth" or "to_time"
    """
    if feature not in ["MDSC_neigh", "CD8+_neigh", "cancer_neigh"]:
        raise ValueError(f"Feature {feature} cannot be used.")

    if alignment not in ["to_birth", "to_time"]:
        raise ValueError(f"Alignment method {alignment} cannot be used.")

    data_list = []
    uid_list = []
    for k, v in df_dict.items():
        uid_list.append(k)
        if alignment == "to_time":
            start = min(v["time"])
            end = max(v["time"])
            data = np.concatenate((np.zeros(start),
                                   v[feature].to_numpy(),
                                   np.zeros(max_time - end)),
                                  axis=None)
        else:
            data = v[feature].to_numpy()

        data_list.append(data)

    if alignment == "to_birth":
        max_len = max([len(row) for row in data_list])
        data_list = np.array([np.pad(row, (0, max_len - len(row))) for row in data_list])
    else:
        data_list = np.array(data_list)

    return data_list, uid_list


def extract_influences(input_var):
    # Extract data for the time point / met-rep combo specified in input variables
    path, met, rep, hour = input_var

    data_cells = pd.read_csv(path + "/cells.csv", usecols=[0, 2, 3, 5], names=["CellID", "x", "y", "state"])
    data_influences = pd.read_csv(path + "/influences.csv", usecols=[0, 1, 2, 3, 4, 5, 6, 7, 9, 11],
                                  names=["CellID", "M0", "M1", "M2", "cancer", "CD4+", "Treg", "CD8+", "NK", "MDSC"])
    data_merged = pd.merge(data_cells, data_influences, on="CellID")
    data_merged["ctype"] = data_merged["state"].map(state_to_name).astype("category")
    data_merged["ctype"] = data_merged["ctype"].cat.set_categories(list(state_to_name.values()))
    data_merged["UID"] = [f"{met}_{rep}_{hour}_{idx}" for idx in data_merged["CellID"]]
    data_merged["rep"] = rep
    data_merged["met"] = met
    data_merged["hour"] = hour
    return data_merged
