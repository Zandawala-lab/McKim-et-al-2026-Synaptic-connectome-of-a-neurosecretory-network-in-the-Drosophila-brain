# %%
# Setup and import packages needed
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import yaml
import copy
import figurefirst as fifi

# pip install git+https://github.com/vanbreugel-lab/braid_tools.git
from braid_analysis import braid_filemanager
from braid_analysis import braid_slicing
from braid_analysis import braid_analysis_plots
from braid_analysis import flymath

# other dependencies for the above that may need to be installed
# pynumdiff, cvxpy, shapely, tables

# %%
# Project folder and path structure ────────────────────────────────────────────
from pathlib import Path

script_dir = Path(__file__).parent          # .../NSC Connectome/scripts/Figure_8
project_root = script_dir.parent.parent     # .../NSC Connectome

input_dir = project_root / "input"
output_dir = project_root / "output" / "Figure_8"

foldername = "Figure8_S1_female_flight_wind-tunnel"
input_file = input_dir / foldername

subfolder_eg4 = "empty-gal4xU9"
subfolder_crz = "CrzxU9"
subfolder_esg4 = "empty-split-gal4xU9"
subfolder_dng27 = "Dng27xU9"

output_dir.mkdir(exist_ok=True)
fig_path = output_dir / "Figure_8_S1DEF"

# %%
# Initialization for figure drawing
# =================================
plt.rcParams['font.sans-serif'] = ['Arial'] + plt.rcParams['font.sans-serif']
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 8
plt.rcParams["ps.usedistiller"] = 'xpdf'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'normal'
plt.rcParams["mathtext.fontset"] = 'cm'

tk_label_sz = 6
ax_stim_sz = 8
ax_label_sz = 6

# Figure/text variables
# =====================
# Figure dimensions
FIG_W, FIG_H = 6.5, 7.0
NROWS, NCOLS = 6, 4

# 0 = display figure only, 1 = display and save figure
save_fig = 1

# (Sub)Figure axis variables & labels
# ===================================
dt = 0.01 # exp frame rate

# Stimulation column labels
stim_txt = ["Flash","Sham"]
# stim_lab_x moved further left (-0.35 -> -0.38) so the 'Flash'/'Sham'
# labels sit clear to the left of the y-axis labels (no overlap), matching
# the original figure. This is an axes-fraction offset, so it is paired with
# the wider axes produced by the tighter wspace below.
stim_lab_x = -0.38
stim_lab_y = 0.35
# General x-axis label & tick positions
t_txt = 'Time relative to flash (s)'
t_ticks = [0, 1, 2, 3, 4, 5]
t_tick_lab = ['0','1','2','3','4','5']

# Fixed genotype txt
fix_gts = " > uas-csChrimson"
# Variable genotype txt
var_gts = ["Crz-gal4", "Empty-gal4", "Dng27-split-gal4", "Empty-split-gal4"]

# Course dir
course_ylab = "Course direction"
course_yticks= [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
# Labels must line up with course_yticks in ascending order; the previous
# ["-π","-π/2","0","π","π/2"] swapped the top two ticks (put "π" at π/2 and
# "π/2" at π). Corrected to match the original figure (top-to-bottom: π, π/2,
# 0, -π/2, -π).
course_ytick_lab = ["-π","-π/2","0","π/2","π"]

# XY-vel
xyspd_ylab = "Horizontal velocity (m/s)"
xyspd_yticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
xyspd_ytick_lab = ['0', '', '', '0.3', '', '', '0.6']

# Altitude
zalt_ylab = "Altitude (m)"
zalt_yticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
zalt_ytick_lab = ['0', '', '', '0.3', '', '', '0.6']

# Data plot options
# =================
ovrlay_type = 'median'

# Figure color options
# ====================
ovrlay_col = 'blue'

# Figure output variables
# =======================
save_type = ".png"

# Dimensions of wind tunnel -- can read these from preprocessing script too
# =========================================================================
xmin = -0.75
xmax = 0.75
ymin = -0.3
ymax = 0.3
zmin = 0
zmax = 0.6


# Adjust function for building filepath from initial setup
def get_filenames_fixed(path, contains, does_not_contain=['~', '.pyc']):
    all_filelist = os.listdir(path)
    filelist = []
    for filename in all_filelist:
        if contains in filename:
            fileok = True
            for nc in does_not_contain:
                if nc in filename:
                    fileok = False
            if fileok:
                filelist.append(os.path.join(path, filename))
    return filelist

braid_filemanager.get_filenames = get_filenames_fixed


# Load datasets
# =============
preprocessed_data_fname_suffix = 'preprocessed_optotrigger_trimmed.hdf'

preprocessed_directory_eg4 = os.path.join(input_file, subfolder_eg4)
preprocessed_directory_crz = os.path.join(input_file, subfolder_crz)
preprocessed_directory_esg4 = os.path.join(input_file, subfolder_esg4)
preprocessed_directory_dng27 = os.path.join(input_file, subfolder_dng27)

preprocessed_data_filename_eg4 = braid_filemanager.get_filename(preprocessed_directory_eg4,
                                                            preprocessed_data_fname_suffix)
preprocessed_data_filename_crz = braid_filemanager.get_filename(preprocessed_directory_crz,
                                                            preprocessed_data_fname_suffix)
preprocessed_data_filename_esg4 = braid_filemanager.get_filename(preprocessed_directory_esg4,
                                                            preprocessed_data_fname_suffix)
preprocessed_data_filename_dng27 = braid_filemanager.get_filename(preprocessed_directory_dng27,
                                                            preprocessed_data_fname_suffix)
print('Loading: ')
print(preprocessed_data_filename_eg4)
df_3d_eg4 = pd.read_hdf(preprocessed_data_filename_eg4)
print('Loading: ')
print(preprocessed_data_filename_crz)
df_3d_crz = pd.read_hdf(preprocessed_data_filename_crz)
print('Loading: ')
print(preprocessed_data_filename_esg4)
df_3d_esg4 = pd.read_hdf(preprocessed_data_filename_esg4)
print('Loading: ')
print(preprocessed_data_filename_dng27)
df_3d_dng27 = pd.read_hdf(preprocessed_data_filename_dng27)


# Filter data
# ===========
# ASSUMPTION --> flash, sham triggers are shared across datasets
trigger_exps = df_3d_eg4.trigger_exp.unique()
obj_id_key = 'obj_id_unique'
flash_key = 'trigger_exp'
flash_val = '2000_0_2000_100'
sham_val = '0_0_0_0'

# Create data variables to plot
# =============================
# empty-gal4
df_3d_eg4_flash = df_3d_eg4[df_3d_eg4[flash_key]==flash_val]
df_3d_eg4_sham = df_3d_eg4[df_3d_eg4[flash_key]==sham_val]

if (ovrlay_type == 'median'):
    eg4_xyspd_ovrl_flash = df_3d_eg4_flash.groupby('time_relative_to_flash').speed_xy.median()
    eg4_xyspd_ovrl_sham = df_3d_eg4_sham.groupby('time_relative_to_flash').speed_xy.median()
    eg4_zalt_ovrl_flash = df_3d_eg4_flash.groupby('time_relative_to_flash').z.median()
    eg4_zalt_ovrl_sham = df_3d_eg4_sham.groupby('time_relative_to_flash').z.median()
else:
    eg4_xyspd_ovrl_sham = df_3d_eg4_sham.groupby('time_relative_to_flash').speed_xy.mean()
    eg4_zalt_ovrl_flash = df_3d_eg4_flash.groupby('time_relative_to_flash').z.mean()
    eg4_zalt_ovrl_sham = df_3d_eg4_sham.groupby('time_relative_to_flash').z.mean()
    
# Crz
df_3d_crz_flash = df_3d_crz[df_3d_crz[flash_key]==flash_val]
df_3d_crz_sham = df_3d_crz[df_3d_crz[flash_key]==sham_val]

if (ovrlay_type == 'median'):
    crz_xyspd_ovrl_flash = df_3d_crz_flash.groupby('time_relative_to_flash').speed_xy.median()
    crz_xyspd_ovrl_sham = df_3d_crz_sham.groupby('time_relative_to_flash').speed_xy.median()
    crz_zalt_ovrl_flash = df_3d_crz_flash.groupby('time_relative_to_flash').z.median()
    crz_zalt_ovrl_sham = df_3d_crz_sham.groupby('time_relative_to_flash').z.median()
else:
    crz_xyspd_ovrl_sham = df_3d_crz_sham.groupby('time_relative_to_flash').speed_xy.mean()
    crz_zalt_ovrl_flash = df_3d_crz_flash.groupby('time_relative_to_flash').z.mean()
    crz_zalt_ovrl_sham = df_3d_crz_sham.groupby('time_relative_to_flash').z.mean()

# empty-split-gal4
df_3d_esg4_flash = df_3d_esg4[df_3d_esg4[flash_key]==flash_val]
df_3d_esg4_sham = df_3d_esg4[df_3d_esg4[flash_key]==sham_val]

if (ovrlay_type == 'median'):
    esg4_xyspd_ovrl_flash = df_3d_esg4_flash.groupby('time_relative_to_flash').speed_xy.median()
    esg4_xyspd_ovrl_sham = df_3d_esg4_sham.groupby('time_relative_to_flash').speed_xy.median()
    esg4_zalt_ovrl_flash = df_3d_esg4_flash.groupby('time_relative_to_flash').z.median()
    esg4_zalt_ovrl_sham = df_3d_esg4_sham.groupby('time_relative_to_flash').z.median()
else:
    esg4_xyspd_ovrl_sham = df_3d_esg4_sham.groupby('time_relative_to_flash').speed_xy.mean()
    esg4_zalt_ovrl_flash = df_3d_esg4_flash.groupby('time_relative_to_flash').z.mean()
    esg4_zalt_ovrl_sham = df_3d_esg4_sham.groupby('time_relative_to_flash').z.mean()

# Dng27
df_3d_dng27_flash = df_3d_dng27[df_3d_dng27[flash_key]==flash_val]
df_3d_dng27_sham = df_3d_dng27[df_3d_dng27[flash_key]==sham_val]

if (ovrlay_type == 'median'):
    dng27_xyspd_ovrl_flash = df_3d_dng27_flash.groupby('time_relative_to_flash').speed_xy.median()
    dng27_xyspd_ovrl_sham = df_3d_dng27_sham.groupby('time_relative_to_flash').speed_xy.median()
    dng27_zalt_ovrl_flash = df_3d_dng27_flash.groupby('time_relative_to_flash').z.median()
    dng27_zalt_ovrl_sham = df_3d_dng27_sham.groupby('time_relative_to_flash').z.median()
else:
    dng27_xyspd_ovrl_sham = df_3d_dng27_sham.groupby('time_relative_to_flash').speed_xy.mean()
    dng27_zalt_ovrl_flash = df_3d_dng27_flash.groupby('time_relative_to_flash').z.mean()
    dng27_zalt_ovrl_sham = df_3d_dng27_sham.groupby('time_relative_to_flash').z.mean()

# %%
# Trajectory counts for figure caption 
# =============
# # Each df_3d_<gt>_<flash/sham> dataframe contains the rows for that genotype x trigger-condition
# 'obj_id_unique_event' are all 3D positions belonging to a single fly
# trajectory through one trigger event, so the number of unique values of
# that column is the number of trajectories

traj_id_col = 'obj_id_unique_event'

trajectory_counts = {
    ('CRZ-Gal4',           'flash'): df_3d_crz_flash[traj_id_col].nunique(),
    ('CRZ-Gal4',           'sham'):  df_3d_crz_sham[traj_id_col].nunique(),
    ('empty-Gal4',         'flash'): df_3d_eg4_flash[traj_id_col].nunique(),
    ('empty-Gal4',         'sham'):  df_3d_eg4_sham[traj_id_col].nunique(),
    ('DNg27-Gal4',         'flash'): df_3d_dng27_flash[traj_id_col].nunique(),
    ('DNg27-Gal4',         'sham'):  df_3d_dng27_sham[traj_id_col].nunique(),
    ('empty-split-Gal4',   'flash'): df_3d_esg4_flash[traj_id_col].nunique(),
    ('empty-split-Gal4',   'sham'):  df_3d_esg4_sham[traj_id_col].nunique(),
}

print("Trajectory counts (panel D, course direction):")
for (gt, cond), n in trajectory_counts.items():
    print(f"  {gt:20s} {cond:6s} n = {n} trajectories")

# %%
# Set up figure
# =============
fig = plt.figure(figsize=(FIG_W, FIG_H))

# All cells will have equal width and height at this level
outer = gridspec.GridSpec(
    NROWS, NCOLS,
    figure=fig,
    left=0.08, right=0.98,
    top=0.97, bottom=0.06,
    # Tighter inter-panel spacing to reproduce the original figure. The grid
    # occupies the same overall span (0.90 wide x 0.91 tall of the figure),
    # so shrinking the gaps enlarges each panel to 94.65 x 63.26 pt (vs the
    # looser 83.41 x 55.59 pt before). Panels are more tightly packed and the
    # 6/8 pt fonts read correctly in proportion to the larger panels.
    hspace=0.25, wspace=0.15,
)
axes = {}  # store axes for later use

# Outer GridSpec: 6 rows x 4 columns
#-----------------------------------
# Row 0, Col 0 - crz course - Flash
# Row 0, Col 1 - egl4 course - Flash
# Row 0, Col 2 - dng27 course - Flash
# Row 0, Col 3 - esgl4 course - Flash
# Row 1, Col 0 - crz course - Sham
# Row 1, Col 1 - egl4 course - Sham
# Row 1, Col 2 - dng27 course - Sham
# Row 1, Col 3 - esgl4 course - Sham
# Row 2, Col 0 - crz xy-vel - Flash
# Row 2, Col 1 - egl4 xy-vel - Flash
# Row 2, Col 2 - dng27 xy-vel - Flash
# Row 2, Col 3 - esgl4 xy-vel - Flash
# Row 3, Col 0 - crz xy-vel - Sham
# Row 3, Col 1 - egl4 xy-vel - Sham
# Row 3, Col 2 - dng27 xy-vel - Sham
# Row 3, Col 3 - esgl4 xy-vel - Sham
# Row 4, Col 0 - crz alt - Flash
# Row 4, Col 1 - egl4 alt - Flash
# Row 4, Col 2 - dng27 alt - Flash
# Row 4, Col 3 - esgl4 alt - Flash
# Row 5, Col 0 - crz alt - Sham
# Row 5, Col 1 - egl4 alt - Sham
# Row 5, Col 2 - dng27 alt - Sham
# Row 5, Col 3 - esgl4 alt - Sham

# Plot figure
# ===========
for row in range(NROWS):
    for col in range(NCOLS):
        # print('Row: ' + str(row) + ', Col: ' + str(col))
        cell = outer[row, col]
        # Specify plot text/variables for current row
        curr_gt = var_gts[col]
        curr_stim_lab = stim_txt[row%2]

        # Specify data vectors to plot
        match row:  # Which variable to plot?
            case 0: # Course-direction, Flash stim
                match col:  # What genotype do we plot?
                    case 0:  # Crz
                        data = df_3d_crz_flash
                    case 1:  # empty-gal4
                        data = df_3d_eg4_flash
                    case 2:  # Dng27
                        data = df_3d_dng27_flash
                    case 3:  # empty split-gal4
                        data = df_3d_esg4_flash
            case 1: # Course-direction, Sham stim
                match col:  # What genotype do we plot?
                    case 0:  # Crz
                        data = df_3d_crz_sham
                    case 1:  # empty-gal4
                        data = df_3d_eg4_sham
                    case 2:  # Dng27
                        data = df_3d_dng27_sham
                    case 3:  # empty split-gal4
                        data = df_3d_esg4_sham
            case 2:  # XY-speed, Flash stim
                match col:  # What genotype do we plot?
                    case 0:  # Crz
                        data = df_3d_crz_flash
                        data_ovrl = crz_xyspd_ovrl_flash
                    case 1:  # empty-gal4
                        data = df_3d_eg4_flash
                        data_ovrl = eg4_xyspd_ovrl_flash
                    case 2:  # Dng27
                        data = df_3d_dng27_flash
                        data_ovrl = dng27_xyspd_ovrl_flash
                    case 3:  # empty split-gal4
                        data = df_3d_esg4_flash
                        data_ovrl = esg4_xyspd_ovrl_flash
            case 3:  # XY-speed, Sham stim
                match col:  # What genotype do we plot?
                    case 0:  # Crz
                        data = df_3d_crz_sham
                        data_ovrl = crz_xyspd_ovrl_sham
                    case 1:  # empty-gal4
                        data = df_3d_eg4_sham
                        data_ovrl = eg4_xyspd_ovrl_sham
                    case 2:  # Dng27
                        data = df_3d_dng27_sham
                        data_ovrl = dng27_xyspd_ovrl_sham
                    case 3:  # empty split-gal4
                        data = df_3d_esg4_sham
                        data_ovrl = esg4_xyspd_ovrl_sham
            case 4:  # Z-pos/Altitude, Flash stim
                match col:  # What genotype do we plot?
                    case 0:  # Crz
                        data = df_3d_crz_flash
                        data_ovrl = crz_zalt_ovrl_flash
                    case 1:  # empty-gal4
                        data = df_3d_eg4_flash
                        data_ovrl = eg4_zalt_ovrl_flash
                    case 2:  # Dng27
                        data = df_3d_dng27_flash
                        data_ovrl = dng27_zalt_ovrl_flash
                    case 3:  # empty split-gal4
                        data = df_3d_esg4_flash
                        data_ovrl = esg4_zalt_ovrl_flash
            case 5:  # Z-pos/Altitude, Sham stim
                match col:  # What genotype do we plot?
                    case 0:  # Crz
                        data = df_3d_crz_sham
                        data_ovrl = crz_zalt_ovrl_sham
                    case 1:  # empty-gal4
                        data = df_3d_eg4_sham
                        data_ovrl = eg4_zalt_ovrl_sham
                    case 2:  # Dng27
                        data = df_3d_dng27_sham
                        data_ovrl = dng27_zalt_ovrl_sham
                    case 3:  # empty split-gal4
                        data = df_3d_esg4_sham
                        data_ovrl = esg4_zalt_ovrl_sham
        # Create new axes
        ax = fig.add_subplot(cell)
        ax.tick_params(labelsize=8)

        match row:
            case 0|1: #Plot Course direction
                # Graph trajec course direction
                braid_analysis_plots.plot_column_vs_time(data,
                                         column='course_smoothish',
                                         time_key='time_relative_to_flash',
                                         norm_columns_to_min_max=True,
                                         norm_columns_to_min_max_smoothing=80,
                                         cmap='bone_r',
                                         vmin=0, vmax=1,
                                         bin_y=np.arange(-np.pi, np.pi + 0.05, 0.05),
                                         bin_x=np.arange(-0.5, 5, 0.03),
                                         ax=ax)
                # If flash stimulus, overlay red opto stimulus on
                if (row == 0):
                    braid_analysis_plots.add_stimulus_shading_to_ax(ax, data, 'lights_on',
                                                    np.pi, -np.pi,
                                                    time_key='time_relative_to_flash',
                                                    obj_id_key='obj_id_unique_event', cmap='bwr',
                                                    vmin=-1 * data.lights_on.max(),
                                                    vmax=data.lights_on.max(), alpha=0.3, zorder=1
                                                    )
                # Course dir graph settings
                ax.set_xticks(t_ticks)
                ax.set_xlabel('')
                ax.tick_params(axis='y', labelsize=tk_label_sz)
                if (col == 0):  # Only label y-axes and stimulus types of 1st column graphs!
                    ax.set_ylabel(course_ylab, size=ax_label_sz, rotation=90)
                    ax.text(stim_lab_x, stim_lab_y, curr_stim_lab, rotation=90, size=ax_stim_sz, transform=ax.transAxes)
                    ax.set_yticks(course_yticks)
                    ax.set_yticklabels(course_ytick_lab)
                    ax.yaxis.set_label_position('left')
                else:
                    ax.set_ylabel('')
                    ax.tick_params(axis='y', which='both', labelleft=False)
                ax.tick_params(axis='x', which='both', labelbottom=False)
                # If on flash row, label graph with GT as title
                if (row == 0):
                    ax.set_title(curr_gt)

            case 2|3: #Plot XY-speed + overlay
                braid_analysis_plots.plot_column_vs_time(data,
                                                        column='speed_xy',
                                                        time_key='time_relative_to_flash',
                                                        norm_columns_to_min_max=True,
                                                        norm_columns_to_min_max_smoothing=20,
                                                        cmap='bone_r',
                                                        vmin=0, vmax=1,
                                                        bin_y=np.arange(0, 0.6 + 0.01, 0.01),
                                                        bin_x=np.arange(-0.5 + dt / 2, 5 + dt + dt / 2, dt),
                                                        ax=ax)
                ax.plot(data_ovrl, color="blue", linewidth=2)

                # Add stimulus shading if Flash stim
                if (row == 2):
                    braid_analysis_plots.add_stimulus_shading_to_ax(ax, data, 'lights_on',
                                                                0.6, 0,
                                                                time_key='time_relative_to_flash',
                                                                obj_id_key='obj_id_unique_event', cmap='bwr',
                                                                vmin=-1 * data.lights_on.max(),
                                                                vmax=data.lights_on.max(), alpha=0.3, zorder=1,
                                                                )
                # Add plot-specific text/labels
                ax.set_xticks(t_ticks)
                ax.set_xlabel('')
                ax.tick_params(axis='y', labelsize=tk_label_sz)
                if (col == 0):
                    ax.set_yticks(xyspd_yticks)
                    ax.set_yticklabels(xyspd_ytick_lab)
                    ax.yaxis.set_label_position('left')
                    ax.set_ylabel(xyspd_ylab, size=ax_label_sz, rotation=90)
                    ax.text(stim_lab_x, stim_lab_y, curr_stim_lab, rotation=90, size=ax_stim_sz,
                                transform=ax.transAxes)
                else:
                    ax.set_ylabel('')
                    ax.tick_params(axis='y', which='both', labelleft=False)
                ax.tick_params(axis='x', which='both', labelbottom=False)
                ax.yaxis.set_label_position('left')

            case 4|5: # Plot Z-pos/Altitude + overlay
                braid_analysis_plots.plot_column_vs_time(data,
                                                 column='z',
                                                 time_key='time_relative_to_flash',
                                                 norm_columns_to_min_max=True,
                                                 norm_columns_to_min_max_smoothing=20,
                                                 cmap='bone_r',
                                                 vmin=0, vmax=1,
                                                 bin_y=np.arange(zmin, zmax + 0.01, 0.01),
                                                 bin_x=np.arange(-0.5 + dt / 2, 5 + dt + dt / 2, dt),
                                                 ax=ax)
                ax.plot(data_ovrl, color="blue", linewidth=2)
                # Add stimulus shading if Flash stim
                if (row == 4):
                    braid_analysis_plots.add_stimulus_shading_to_ax(ax, data, 'lights_on',
                                                        zmax, zmin,
                                                        time_key='time_relative_to_flash',
                                                        obj_id_key='obj_id_unique_event', cmap='bwr',
                                                        vmin=-1 * data.lights_on.max(),
                                                        vmax=data.lights_on.max(), alpha=0.3, zorder=1,
                                                        )
                # Add plot-specific text/labels
                ax.set_xticks(t_ticks)
                ax.set_xticklabels(t_tick_lab, size=tk_label_sz)

                if (row % 2 != 0):  # Label x-axes on last row
                    ax.set_xlabel(t_txt, size=ax_label_sz)
                else:
                    ax.set_xlabel('')
                    ax.tick_params(axis='x', which='both', labelbottom=False)

                ax.tick_params(axis='y', labelsize=tk_label_sz)
                if (col == 0):
                    ax.set_yticks(zalt_yticks)
                    ax.set_yticklabels(zalt_ytick_lab,size=tk_label_sz)
                    ax.set_ylabel(zalt_ylab, size=ax_label_sz, rotation=90)
                    ax.text(stim_lab_x, stim_lab_y, curr_stim_lab, rotation=90, size=ax_stim_sz,
                                transform=ax.transAxes)
                else:
                    ax.set_ylabel('')
                    ax.tick_params(axis='y', which='both', labelleft=False)


        # Save axes
        axes[(row, col)] = ax

# To make the spines look nicer, you can run this on each axis (adjust xticks, yticks, etc. as needed)
# fifi.mpl_functions.adjust_spines(ax, ['left', bottom'], xticks=[0, 1, 2, 3, 4, 5], spine_locations={'left': 5, bottom': 5}, tick_length=3, linewidth=1)

# Save figure
# ===========
if (save_fig):
    # output format is switched simply by editing save_type (".png" or ".svg").
    # bbox_inches='tight' is required: it is what the original figure used, and
    # it crops to the tight content bounds that give the matched viewBox size.
    savefile = fig_path.with_suffix(save_type)
    if save_type.lower() == ".png":
        fig.savefig(savefile, dpi=150, bbox_inches='tight')
    else:
        fig.savefig(savefile, bbox_inches='tight')
    print(f"Saved figure to: {savefile}")
plt.show()
