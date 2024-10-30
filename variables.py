##########################################################
### variables that can be changed for every simulation ###
##########################################################
# all units are given in µm!

import numpy as np

# name of this simulation
global_name = 'Example_Name3D'
# path to the folder the results should be saved in
path = r'F:\Example_Name3D'
# path to this variables.py file
variables_path = rf"C:\PathToProject\variables.py"


# Field of view in x direction in [µm]
FOV_x = 5
# Field of view in y direction in [µm]
FOV_y = 5


# protein density in [#/µm^2]
prot_dens = 10
# lipid density in [#/µm^2]
lip_dens = 10


# mu, alpha and scale factor of scipy statistics skew distribution
# average glycan height in [µm] (careful! This is not exactly the maximum in the scipy skew distribution!)
# for mucins, only lower and upper glycan height is important
mu_glycan_height = 0.008
# alpha value of glycan height
alpha_glycan_height = 2.5
# scale factor of glycan height
scale_glycan_height = 0.006
# lower bound glycan height in [µm]
lower_glycan_height = 0.002
# upper bound glycan height in [µm]
upper_glycan_height = 0.025


# probabilities for 0, 1, 2, ... glycans per protein
prob_x_glycoprots = np.array([0.3, 0.7])

# probabilities for 0, 1, 2, 3, 4, ... sias per glycan on glycoprotein (-> bi-/... antennary glycans)
prob_x_sias_prots = np.array([0.3, 0.7])
#
# probability for 0, 1 glycan per lipid
prob_x_glycolips = np.array([0.6, 0.4])
# probabilities for 0 or 1 sia per glycan on glycolipid
prob_x_sias_lips = np.array([0.65, 0.35])
#
#
# average distance between sias and their "glycoprotein"
mu_sias_glycoprots = 0.012
# standard deviation of the distance between glycans and their sias
sigma_sias_glycoprots = 0.002
#
#
# forbidden region around one glycan
forbidden_angle_glycans = np.deg2rad(10)
#
#
# average distance between two sias
mu_sias = 0.0026
# standard deviation of the distance between two sias
sigma_sias = 0.0008
#
cutoff_gaussians = 2
#
#

