from variables import *


def check_variables():
    """ checks whether the variables from variables.py are reasonable."""
    err = 0
    if (FOV_x < 0) or (FOV_y < 0):
        print("Your FOV has to be positive!")
        err += 1

    if (prot_dens < 0) or (lip_dens < 0):
        print("Your protein and lipid density have to be positive!")
        err += 1

    if (round(np.sum(prob_x_glycoprots), 6)) != 1 or (round(np.sum(prob_x_glycolips), 6) != 1):
        print("The probabilities for 0, 1, 2, ... glycans per protein and lipid have to sum up to 1!")
        err += 1

    if(round(np.sum(prob_x_sias_prots), 6)) != 1 or (round(np.sum(prob_x_sias_lips), 6) != 1):
        print("The probabilities for 0, 1, 2, ... sias per protein and lipid have to sum up to 1!")
        err += 1

    if (round(np.sum(prob_x_sias_prots), 6) != 1) or (round(np.sum(prob_x_sias_lips), 6) != 1):
        print("The probabilities for 0, 1, 2, ... sias per protein and lipid have to sum up to 1!")
        err += 1

    if(mu_sias or mu_sias_glycoprots) < 0:
        print("Your distances have to be positive!")
        err += 1

    if abs(0.5*(mu_sias + cutoff_gaussians*sigma_sias)/(mu_sias_glycoprots - cutoff_gaussians*sigma_sias_glycoprots)) > 1:
        print("Your distance between sias can get too large compared to the distance between sias and proteins."
              "Please consider either changing the gaussian cutoff, mean or standard deviation of one of the "
              "distributions!")
        err += 1

    if (mu_glycan_height or lower_glycan_height or upper_glycan_height) < 0:
        print("Your protein height as well as its distribution's lower and upper limit should be positive values!")
        err += 1

    if alpha_glycan_height < 0:
        print("Please note that your distribution is shifted to the right side. The code will run nonetheless. ")

    if err == 0:
        print("Your variables are valid :)")
    else:
        print("Please check your variables. The run was terminated.")
        quit()
