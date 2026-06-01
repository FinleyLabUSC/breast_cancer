"""
TK 06.01.2026
"""
from analyzeABM.constants import *
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
import matplotlib.pyplot as plt


def extract_neighbors(input_var):
    """
    :param input_var: list containing path to data, met number, rep number, and hour
    :return: dataframe with celltype, cell UID, rep/met/hour info, and neighbors vector
    """
    # Extract data for this time point
    path, met, rep, hour = input_var
    data = pd.read_csv(path, usecols=[0, 2, 3, 5], names = ["CellID", "x", "y", "state"])
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
        nvects.append(data.iloc[idx].ctype.value_counts().sort_index().to_numpy()) # This gives ordered list of neighbors
    data["neighbors"] = nvects

    return data

def retrieve_ts(folder_name, num_replicates = 5, max_len = 600):
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
    cols_to_get = ["m0", "m1", "m2", "c", "cd4_th", "cd4_treg", "cd8", "nk", "mdsc", "myeloid","lymphoid","stromal","radius"]

    # Create all the files
    for i, file_name in enumerate(files):
        df = pd.read_csv(file_name, header = 0)

        # Append final state until max length
        last_time = max(df.index)
        max_index = last_time
        while max_index < max_len:
            df.loc[max_index + 1] = df.loc[last_time]
            max_index += 1

        if i == 0:
            for col in cols_to_get:
                dataframe_dict[col] = df[[col]].copy()
                dataframe_dict[col].rename({col : "0"}, axis=1, inplace=True)
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
        dc.append(cpop_list[i+1]-cpop_list[i])

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
                dataframe_dict[col].rename({col : "0"}, axis=1, inplace=True)
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

def retrieve_kill_ts(folder_name, num_replicates = 5, max_len = 600):
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

    dataframe_dict = {};
    cols_to_get = ["age", "cd8", "nk"]

    # Create all the files
    for i, file_name in enumerate(files):
        df = pd.read_csv(file_name, header = 0)

        # Append final state until max length
        last_time = max(df.index)
        max_index = last_time
        while max_index < max_len:
            df.loc[max_index + 1] = df.loc[last_time]
            max_index += 1

        if i == 0:
            for col in cols_to_get:
                dataframe_dict[col] = df[[col]].copy()
                dataframe_dict[col].rename({col : "0"}, axis=1, inplace=True)
        else:
            for col in cols_to_get:
                dataframe_dict[col][str(i)] = df[col].to_numpy()

    return dataframe_dict

