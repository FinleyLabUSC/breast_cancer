"""
TK 06.01.2026
"""
import matplotlib.colors as mcolors

# Cell states, names, numerical labels, and colors
cell_colors = [
    'xkcd:lilac', # M0
    'xkcd:bright purple', # M1
    'tab:pink', # M2
    'tab:blue', # cancer
    'tab:green', # CD4 T helper
    'lime', # CD4 TReg
    'tab:orange', # CD8
    'tab:brown', # NK cell
    'xkcd:goldenrod', # MDSC
    'xkcd:light pink', # Myeloid
    'xkcd:beige', # Lymphoid
    'darkgray' # Stromal
]
dict_cell_colors = {
    "M0" : 'xkcd:lilac', # M0
    "M1" : 'xkcd:bright purple', # M1
    "M2" : 'tab:pink', # M2
    "cancer" : 'tab:blue', # cancer
    "CD4+" : 'tab:green', # CD4 T helper
    "Treg" : 'lime', # CD4 Treg
    "CD8+" : 'tab:orange', # CD8
    "NK" : 'tab:brown', # NK cell
    "MDSC" : 'xkcd:goldenrod', # MDSC
    "myeloid" : 'xkcd:light pink', # Myeloid
    "lymphoid" : 'xkcd:beige', # Lymphoid
    "stromal" : 'darkgray' # Stromal
}
cell_cmap = mcolors.ListedColormap(cell_colors)
state_to_name = {
    0:"M0",
    1:"M1",
    2:"M2",
    3:"cancer",
    4:"CD4+ T",
    5:"Treg",
    6:"CD8+ T",
    8:"NK",
    10:"MDSC",
    11:"myeloid",
    12:"lymphoid",
    13:"stromal"
}
name_to_state = {value: key for key, value in state_to_name.items()}
cellorder = ["M0", "M1", "M2", "cancer", "CD4+", "Treg", "CD8+", "NK", "MDSC", "myeloid", "lymphoid", "stromal"]
cellorder_nogen = ["M0", "M1", "M2", "cancer", "CD4+", "Treg", "CD8+", "NK", "MDSC"]
mets = [100, 106, 107, 108, 112, 113, 114, 117, 122, 123, 125, 128, 130, 132, 134, 140, 142, 143, 145]

# For display only
cell_labels = ["M0", "M1", "M2", "Cancer", "Th", "TReg", "CD8", "NK", "MDSC", "Myeloid", "Lymphoid", "Stromal"]
num_cell_types = 12
