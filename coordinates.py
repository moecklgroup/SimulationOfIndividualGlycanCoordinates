#############################################################
### Define several functions to obtain random coordinates ###
#############################################################

# import packages
import numpy as np
import random
import scipy

import variables as v


def get_cut_skew_normal_random(loc, alpha, scale, lower_bound, upper_bound, point_nr=1):
    """
    Generate a random number from a cut skew normal distribution.

    Parameters:
    -----------------------
    loc:            float
                    approximate mean of the distribution
    alpha:          float
                    alpha parameter of skewnorm
    scale:          float
                    Scale parameter of skewnorm
    lower_bound:    float
                    lower limit of the distribution
    upper_bound:    float
                    upper limit of the distribution

    Returns
    ----------------------
    valid_points:   numpy array
                    random numbers from the skew normal distribution
    """
    skew_normal = scipy.stats.skewnorm(a=alpha, loc=loc, scale=scale)
    valid_points = []

    while len(valid_points) < point_nr:
        random_numbers = skew_normal.rvs(size=(point_nr - len(valid_points)))
        valid_points.extend([num for num in random_numbers if lower_bound <= num <= upper_bound])

    return valid_points


def random_coordinates(point_nr, x_size, y_size):
    """ creates a 2d numpy array with pairs of random coordinates that lie within a certain region.

        Parameters
        -----------
        point_nr:   int
                    number of points that should be added to the array
        x_size:     float
                    length in x-direction of the possible region for coordinates
        y_size:     float
                    length in y-direction of the possible region for coordinates

        Returns
        ----------
        coordinates_rand:    2d numpy array with random coordinates
    """
    # create the random x- and y-coordinates
    x_values = np.random.uniform(0, x_size, point_nr)
    y_values = np.random.uniform(0, y_size, point_nr)

    # stack them on top of each other
    coordinates_rand = np.column_stack([x_values, y_values])

    return coordinates_rand


def random_coordinates_3d(point_nr, x_size, y_size, mu_z, alpha_z, scale_z, lower_z, upper_z):
    """ creates a 2d numpy array with pairs of random coordinates that lie within a certain region.

        Parameters
        -----------
        point_nr:   int
                    number of points that should be added to the array
        x_size:     float
                    length in x-direction of the possible region for coordinates
        y_size:     float
                    length in y-direction of the possible region for coordinates
        mu_z:       float
                    mean z-value of coordinates
        alpha_z:    float
                    alpha value for coordinate height
        scale_z:    float
                    scale value for coordinate height
        lower_z:    float
                    lower limit for coordinate height
        upper_z:    float
                    upper limit for coordinate height

        Returns
        ----------
        coordinates_rand:    2d numpy array with random coordinates
    """
    # create the random x- and y-coordinates
    x_values = np.random.uniform(0, x_size, point_nr)
    y_values = np.random.uniform(0, y_size, point_nr)
    if mu_z == 0:
        z_values = np.zeros(point_nr)
    else:
        z_values = get_cut_skew_normal_random(mu_z, alpha_z, scale_z, lower_z, upper_z, point_nr)
    # stack them on top of each other
    coordinates_rand = np.column_stack([x_values, y_values, z_values])

    return coordinates_rand


def coordinates_grouped_around(ground_coords, prob_groups, prob_in_group, x_size, y_size, mu1, sig1, mu_sias,
                               sigma_sias, forbidden_angle=np.deg2rad(10)):
    """ generates point groups around coordinates of the ground_coords 2d numpy array.
       Parameters
       ----------
       ground_coords:       2d numpy array
                            array with the coordinates around which the new ones will be placed
       prob_groups:         np array
                            probability for 0, 1, 2, 3, ... new groups around one of the original coordinates
       prob_in_group:       np array
                            probability for 0, 1, 2, 3, ... coordinates within one new group
       x_size:              float
                            length in x-direction of the possible region for coordinates
       y_size:              float
                            length in y-direction of the possible region for coordinates
       mu1:                 float
                            mean of gaussian distribution for the distance between the original and one new group
                            Gaussian is cut at cutoff_gaussians (from variables.py) st. dev. left and right of the mean
       sig1:                float
                            sigma of gaussian distribution for the distance between the original and one new group
       mu_sias:             float
                            mean of gaussian distribution (if != 0) for the average angle between several new points
                            if mu2 == 0 -> take a uniform distribution -> random distribution of angles
                            Gaussian is cut at cutoff_gaussians (from variables.py) st. dev. left and right of the mean
       sigma_sias:          float
                            sigma of gaussian distribution for the distance between several new points
       forbidden_angle:     float
                            angle that should be left free left and right of a group

       Returns
       -------
       coordinates         2d numpy array
                           with the new coordinates

       """
    coordinates = []  # create an empty array

    if prob_groups[0] == 1 or prob_in_group[0] == 1:
        return None
    else:
        random_values = np.random.rand(len(ground_coords))

        for i in range(len(ground_coords)):  # for every original coordinate
            rand = random_values[i]  # produce a random number between 0 and 1
            prob = 0  # and set the probability value to 0

            # produce random radii in bulk
            nr_new_radii = 0.5 * len(prob_groups) * (len(prob_groups) + 1)
            radii = scipy.stats.truncnorm.rvs(-v.cutoff_gaussians, v.cutoff_gaussians, loc=mu1, scale=sig1,
                                              size=int(nr_new_radii))
            r_index = 0

            for j in range(len(prob_groups)):
                if rand < prob + prob_groups[j]:
                    # If the random number is within the probability for 0, 1, 2, ... groups generate j groups
                    angles_allowed = np.array([[0, 2 * np.pi]])

                    for _ in np.arange(j):
                        # for each group (glycan), determine a radius and an angle
                        r = radii[r_index]
                        r_index += 1

                        angle, angles_allowed = get_allowed_angle(angles_allowed, forbidden_angle)

                        # generate added angles in bulk for each group (glycan)
                        # -> radius of sias to protein should be the same for all of one protein
                        nr_new_angles = len(prob_in_group)
                        if mu_sias == 0:
                            angle_add = np.random.uniform(0, 2 * np.pi, nr_new_angles)
                        else:
                            sias_dis = scipy.stats.truncnorm.rvs(-v.cutoff_gaussians, v.cutoff_gaussians, loc=mu_sias,
                                                                 scale=sigma_sias, size=int(nr_new_angles))
                            angle_add = 2 * np.arcsin(0.5 * sias_dis / r)
                        angle_index = 0

                        # for each group (glycan), generate a random number to determine the number of sialin acids on it
                        rand2 = np.random.rand()
                        prob2 = 0

                        for k in range(len(prob_in_group)):
                            if rand2 < prob2 + prob_in_group[k]:
                                # generate k coordinates (sias) for that group (glycan)
                                for pt in range(k):
                                    x = r * np.cos(angle)
                                    y = r * np.sin(angle)
                                    new_point = [ground_coords[i, 0] + x, ground_coords[i, 1] + y]
                                    # if the random point lies within the FOV, take it as your new point
                                    if 0 <= (ground_coords[i, 0] + x) <= x_size and 0 <= (
                                            ground_coords[i, 1] + y) <= y_size:
                                        coordinates.append(new_point)

                                    # change the angle for the next sia and edit the allowed angles for the next group
                                    angle += angle_add[angle_index]
                                    angles_allowed = get_allowed_angles_left(angle, forbidden_angle, angles_allowed)
                                    angle_index += 1
                                break

                            else:
                                prob2 += prob_in_group[k]
                    break
                else:
                    prob += prob_groups[j]

    return np.array(coordinates)


def get_allowed_angle(former_allowed_angles, forbidden_angle):
    """ generates point groups around coordinates of the ground_coords 2d numpy array.

           Parameters
           ----------
           former_allowed_angles:   2d numpy array
                                    array with the ranges of allowed angles
           forbidden_angle:         float
                                    angle which has to be kept empty left and right of one group of coordinates

           Returns
           -------
           angle:                   float
                                    new angle
           angles_allowed_new:      2d numpy array
                                    array with the ranges of still allowed angles
    """
    degrees_left = 0
    # calculate the angle range that could still be hit
    for a in np.arange(len(former_allowed_angles)):
        degrees_left += (former_allowed_angles[a, 1] - former_allowed_angles[a, 0])
    if degrees_left == 0:
        print('Your blocked range around of ' + str(forbidden_angle) + ' could not be fulfilled.')
        degrees_left = 2 * np.pi
        former_allowed_angles = np.array([[0, 2 * np.pi]])
    # draw a random angle
    random_angle = np.random.uniform(0, degrees_left)

    # find the actual angle
    for i in np.arange(len(former_allowed_angles)):
        if former_allowed_angles[i, 0] + random_angle <= former_allowed_angles[i, 1]:
            # if the "random angle" is in the current range, determine the real angle and stop
            angle = former_allowed_angles[i, 0] + random_angle
            break
        else:
            # otherwise go to the next interval
            random_angle -= (former_allowed_angles[i, 1] - former_allowed_angles[i, 0])

    # determine the new array for allowed angles
    angles_allowed_new = get_allowed_angles_left(angle, forbidden_angle, former_allowed_angles)

    return angle, angles_allowed_new


def get_allowed_angles_left(angle, forbidden_angle, former_allowed_angles):
    """ Calculates the new array with sequences of allowed angles.

           Parameters
           ----------
           angle:                   float
                                    between 0 and 2*pi; gives the new angle
           forbidden_angle:         float
                                    angle which has to be kept empty left and right of one group of coordinates
           former_allowed_angles:   2d numpy array
                                    array with the ranges of allowed angles

           Returns
           -------
           angles_allowed_new:      2d numpy array
                                    array with the ranges of still allowed angles
    """
    # determine the ranges of theoretically allowed angles left if there would be only the new angle
    if ((angle - forbidden_angle) % (2 * np.pi)) < ((angle + forbidden_angle) % (2 * np.pi)):
        hypothetically_allowed = np.array([[0, angle - forbidden_angle], [angle + forbidden_angle, 2 * np.pi]])
    else:
        hypothetically_allowed = np.array(
            [[(angle + forbidden_angle) % (2 * np.pi), (angle - forbidden_angle) % (2 * np.pi)]])

    # Determine the intersections between the 2 arrays
    # Expand dimensions, so that one can do broadcasting
    former_allowed_angles = former_allowed_angles[:, None, :]
    hypothetically_allowed = hypothetically_allowed[None, :, :]

    # Find the lower and upper bounds of the intervals
    lower_bound = np.maximum(former_allowed_angles[:, :, 0], hypothetically_allowed[:, :, 0])
    upper_bound = np.minimum(former_allowed_angles[:, :, 1], hypothetically_allowed[:, :, 1])

    # Put valid intervals together
    mask = lower_bound <= upper_bound
    angles_allowed_new = np.column_stack([lower_bound[mask], upper_bound[mask]])

    return angles_allowed_new


def coordinates_reduce(coords, point_number_reduced):
    """ reduces a 2d array to a shorter one with point_number_reduced points.

    Parameters
    ----------
    coords                  2d numpy array
                            array with the original coordinates
    point_number_reduced    int
                            number of points left after reducing the array

    Returns
    -------
    coordinates_reduced     2d numpy array
                            with the coordinates left after reducing
    """
    # points are taken out randomly (necessary if coordinates have preferred locations)

    dimension = len(coords[0])

    if point_number_reduced > len(coords):
        print(
            "Your point number taking labeling efficiency into account is greater than the length of your input array. Please check for mistakes!")
        return coords

    coordinates_reduced = np.empty(
        shape=(point_number_reduced, dimension))  # empty array; length takes labeling efficiency into account
    # indices of the coordinates that should be kept
    indices = random.sample(range(len(coords)), point_number_reduced)
    reduced_index = 0
    for index in indices:
        coordinates_reduced[reduced_index] = coords[index]
        reduced_index += 1
    return coordinates_reduced


def spherical_to_cartesian(r, theta, phi):
    """
    converts spherical to cartesian coordinates.

    Parameters
    ----------
    r:      float
            radius
    theta:  float between 0 and pi
            polar angle
    phi:    float between 0 and 2*pi

    Returns
    --------
    x:      float
            x-coordinate
    y:      float
            y-coordinate
    z:      float
            z-coordinate
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def cartesian_to_spherical(x, y, z):
    """
    converts cartesian to spherical coordinates.

    Parameters
    ----------
    x:      float
            x-coordinate
    y:      float
            y-coordinate
    z:      float
            z-coordinate

    Returns
    --------
    r:      float
            radius
    theta:  float between 0 and pi
            polar angle
    phi:    float between 0 and 2*pi

    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    phi = phi if phi >= 0 else phi + 2 * np.pi  # to account for phi in range (0, 2pi)
    return r, theta, phi


def angle_between_vectors(v1, v2):
    """
    Calculates the angle (in radians) between two vectors.

    Parameter
    -----------------
    v1:     numpy array
            first vector
    v2:     numpy array
            second vector

    Returns
    -----------------
    angle:  float
            radian between the vectors
    """
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    cos_angle = dot_product / (norm_v1 * norm_v2)
    angle = np.arccos(np.clip(cos_angle, -1, 1))
    return angle


def generate_new_spherical_mask(Theta, Phi, mask, picked_point, cone_angle):
    """
    Update the mask to filter out points and their neighborhoods within a cone-shaped region around the picked point.

    Parameters
    -------------------
    Theta:      numpy meshgrid
                covers the theta part of the polar coordinates
    Phi:        numpy meshgrid
                covers the phi part of the polar coordinates
    mask:       numpy array
                mask for allowed polar coordinates

    Returns
    ------------------
    mask:       numpy array
                new mask
    """
    # Get the shape of the meshgrid
    num_rows, num_cols = Theta.shape

    # Convert picked point to spherical coordinates
    x_picked, y_picked, z_picked = picked_point

    # Update the mask to filter out points and their neighborhoods within the cone-shaped region around the picked point
    for i in range(num_rows):
        for j in range(num_cols):
            # Convert current point to spherical coordinates
            x_curr, y_curr, z_curr = spherical_to_cartesian(1, Theta[i, j], Phi[i, j])
            # Calculate the angle between the picked point's direction and the current point
            angle = angle_between_vectors(np.array([x_picked, y_picked, z_picked]), np.array([x_curr, y_curr, z_curr]))
            # If the angle is smaller than the cone angle, exclude the current point
            if angle < cone_angle:
                mask[i, j] = False

    if np.all(mask == False):
        print('The condition of forbidden angles could not be fulfilled. All the angles were set to allowed again.')
        mask = np.ones_like(Theta, dtype=bool)

    return mask


def random_spherical_point_from_meshgrid(Theta, Phi, mask):
    """
    Pick a random point from a meshgrid representing points on the surface of a sphere,
    excluding the points that have already been picked and their neighborhoods within a cone-shaped region.

    Parameters
    --------------
    Theta:          numpy meshgrid
                    covers the theta part of the polar coordinates
    Phi:            numpy meshgrid
                    covers the phi part of the polar coordinates
    mask:           numpy array
                    mask for allowed polar coordinates

    Returns
    ------------------
    random_theta:   float
                    randomly picked theta from allowed values
    random_phi:     float
                    randomly picked phi from allowed values
    """
    # Get the available choices (indices) for picking a random point
    available_choices = np.argwhere(mask)

    # Choose a random index from the available choices
    random_index = np.random.choice(len(available_choices))
    random_row, random_col = available_choices[random_index]

    # Get the corresponding coordinates
    random_theta = Theta[random_row, random_col]
    random_phi = Phi[random_row, random_col]

    return random_theta, random_phi


def rotate_point_around_z_axis(point, phi):
    """
    Rotate a point around the z-axis by an angle phi.

    Parameters
    --------------
    point:      numpy array
                Coordinates of the point (x, y, z)
    phi:        float
                Angle of rotation in radians

    Returns
    --------------
    rotated coordinates of the point (x', y', z')
    """
    rotation_matrix = np.array([[np.cos(phi), -np.sin(phi), 0],
                                [np.sin(phi), np.cos(phi), 0],
                                [0, 0, 1]])
    return np.dot(rotation_matrix, point)


def rotate_point_around_y_axis(point, theta):
    """
    Rotate a point around the y-axis by an angle theta.

    Parameters
    --------------
    point:      numpy array
                Coordinates of the point (x, y, z)
    theta:      float
                Angle of rotation in radians

    Returns
    --------------
    rotated coordinates of the point (x', y', z')
    """
    rotation_matrix = np.array([[np.cos(theta), 0, -np.sin(theta)],
                                [0, 1, 0],
                                [np.sin(theta), 0, np.cos(theta)]])
    return np.dot(rotation_matrix, point)


def rotate_point_theta_phi(point, theta, phi):
    """
    Rotate a point around the y-axis by an angle theta and around the z-axis by an angle phi

    Parameters
    --------------
    point:      numpy array
                Coordinates of the point (x, y, z)
    theta:      float
                Angle of rotation around y-axis in radians
    phi:        float
                Angle of rotation around z-axis in radians

    Returns
    --------------
    rotated coordinates of the point (x', y', z')
    """
    new_point = rotate_point_around_y_axis(point, -theta)
    final_point = rotate_point_around_z_axis(new_point, phi)

    return final_point


def coordinates_grouped_around_3d(ground_coords, prob_groups, prob_in_group, x_size, y_size, mu_radius_group,
                                  sig_radius_group, mu_sias,
                                  sigma_sias, forbidden_angle=np.deg2rad(10)):
    """ generates point groups around coordinates of the ground_coords 2d numpy array.
       Parameters
       ----------
       ground_coords:       2d numpy array
                            array with the coordinates around which the new ones will be placed
       prob_groups:         np array
                            probability for 0, 1, 2, 3, ... new groups around one of the original coordinates
       prob_in_group:       np array
                            probability for 0, 1, 2, 3, ... coordinates within one new group
       x_size:              float
                            length in x-direction of the possible region for coordinates
       y_size:              float
                            length in y-direction of the possible region for coordinates
       mu_radius_group:     float
                            mean of gaussian distribution for the distance between the original point and one new group
                            Gaussian is cut at cutoff_gaussians (from variables.py) st. dev. left and right of the mean
       sig_radius_group:    float
                            sigma of gaussian distribution for the distance between the original point and one new group
       mu_sias:             float > 0
                            mean of gaussian distribution for the average angle between several new points
                            Gaussian is cut at cutoff_gaussians (from variables.py) st. dev. left and right of the mean
       sigma_sias:          float
                            sigma of gaussian distribution for the distance between several new points
       forbidden_angle:     float
                            angle that should be left free left and right of a group

       Returns
       -------
       coordinates         2d numpy array
                           with the new coordinates

       """
    coordinates = []  # create an empty array

    if prob_groups[0] == 1 or prob_in_group[0] == 1:  # if there are no sias, return an empty array
        return None
    else:
        random_values = np.random.rand(len(ground_coords))

        for i in range(len(ground_coords)):  # for every original coordinate
            rand = random_values[i]  # use one of the random numbers between 0 and 1
            prob = 0  # and set the probability value to 0

            # produce random radii in bulk
            nr_new_radii = 0.5 * len(prob_groups) * (len(prob_groups) + 1)
            radii = scipy.stats.truncnorm.rvs(-v.cutoff_gaussians, v.cutoff_gaussians, loc=mu_radius_group,
                                              scale=sig_radius_group,
                                              size=int(nr_new_radii))
            r_index = 0

            for j in range(len(prob_groups)):
                if rand < prob + prob_groups[j]:
                    # If the random number is within the probability for 0, 1, 2, ... groups generate j groups
                    # Define the spherical grid
                    theta = np.linspace(0, 0.5 * np.pi, 100)  # Polar angle (latitude), only upper half sphere
                    phi = np.linspace(0, 2 * np.pi, 400)  # Azimuthal angle (longitude)
                    Theta, Phi = np.meshgrid(theta, phi)

                    mask = np.ones_like(Theta,
                                        dtype=bool)  # create a mask that saves which points on the sphere are still available

                    for _ in np.arange(j):
                        # for each group (glycan), determine a radius
                        r = radii[r_index]
                        r_index += 1

                        # generate added angles/points in bulk for each group (glycan)
                        # -> radius of sias to protein should be the same for all of one protein
                        nr_new_angles = len(prob_in_group) - 2
                        sias_dis = scipy.stats.truncnorm.rvs(-v.cutoff_gaussians, v.cutoff_gaussians, loc=mu_sias,
                                                             scale=sigma_sias, size=int(nr_new_angles))

                        theta_unrotated = np.arccos((2 * r ** 2 - sias_dis ** 2) / (2 * r ** 2))
                        phi_unrotated = np.random.uniform(0, 2 * np.pi, int(nr_new_angles))
                        points_unrotated = np.transpose(spherical_to_cartesian(r, theta_unrotated, phi_unrotated))
                        point_index = 0

                        # for each group (glycan), generate a random number to determine the number of sialin acids on it
                        rand2 = np.random.rand()
                        prob2 = 0

                        for k in range(len(prob_in_group)):
                            if rand2 < prob2 + prob_in_group[k]:
                                # generate k coordinates (sias) for that group (glycan)
                                current_rotating_angle = 0
                                for pt in range(k):
                                    if pt == 0:  # create random (possible) angles for the first point and update the mask of possible angles
                                        random_theta, random_phi = random_spherical_point_from_meshgrid(Theta, Phi,
                                                                                                        mask)  # will be in radians
                                        mask = generate_new_spherical_mask(Theta, Phi, mask,
                                                                           spherical_to_cartesian(1, random_theta,
                                                                                                  random_phi),
                                                                           forbidden_angle)
                                        # convert the point to cartesian coordinates and save the theta as rotation angle
                                        new_x, new_y, new_z = spherical_to_cartesian(r, random_theta, random_phi)
                                        current_rotating_angle_theta = random_theta  # in rad
                                        current_rotating_angle_phi = random_phi  # in rad
                                    else:
                                        # use one of the points that were used around the center and rotate them accordingly
                                        new_x, new_y, new_z = rotate_point_theta_phi(points_unrotated[point_index],
                                                                                     current_rotating_angle_theta,
                                                                                     current_rotating_angle_phi)
                                        point_index += 1
                                        _, current_rotating_angle_theta, current_rotating_angle_phi = cartesian_to_spherical(
                                            new_x, new_y, new_z)
                                        mask = generate_new_spherical_mask(Theta, Phi, mask, spherical_to_cartesian(1,
                                                                                                                    current_rotating_angle_theta,
                                                                                                                    current_rotating_angle_phi),
                                                                           forbidden_angle)

                                    # translate the coordinates from ground coordinate center to original coordinate system
                                    new_point = [ground_coords[i, 0] + new_x, ground_coords[i, 1] + new_y,
                                                 ground_coords[i, 2] + new_z]
                                    # if the random point lies within the FOV, take it as your new point
                                    if 0 <= (ground_coords[i, 0] + new_x) <= x_size and 0 <= (
                                            ground_coords[i, 1] + new_y) <= y_size and ground_coords[i, 2] + new_z > 0:
                                        coordinates.append(new_point)
                                break

                            else:
                                prob2 += prob_in_group[k]
                    break
                else:
                    prob += prob_groups[j]

    return np.array(coordinates)


def coordinates_grouped_around_3d_mucinextension(ground_coords, prob_groups, prob_in_group, x_size, y_size,
                                                 mu_radius_group, sig_radius_group, mu_sias,
                                                 sigma_sias, lower_glycan_height=v.lower_glycan_height,
                                                 upper_glycan_height=v.upper_glycan_height,
                                                 forbidden_angle=np.deg2rad(10)):
    """ generates point groups around coordinates of the ground_coords 2d numpy array.
       Parameters
       ----------
       ground_coords:       2d numpy array
                            array with the coordinates around which the new ones will be placed
       prob_groups:         np array
                            probability for 0, 1, 2, 3, ... new groups around one of the original coordinates
       prob_in_group:       np array
                            probability for 0, 1, 2, 3, ... coordinates within one new group
       x_size:              float
                            length in x-direction of the possible region for coordinates
       y_size:              float
                            length in y-direction of the possible region for coordinates
       mu_radius_group:     float
                            mean of gaussian distribution for the distance between the original point and one new group
                            Gaussian is cut at cutoff_gaussians (from variables.py) st. dev. left and right of the mean
       sig_radius_group:    float
                            sigma of gaussian distribution for the distance between the original point and one new group
       mu_sias:             float > 0
                            mean of gaussian distribution for the average angle between several new points
                            Gaussian is cut at cutoff_gaussians (from variables.py) st. dev. left and right of the mean
       sigma_sias:          float
                            sigma of gaussian distribution for the distance between several new points
       forbidden_angle:     float
                            angle that should be left free left and right of a group

       Returns
       -------
       coordinates         2d numpy array
                           with the new coordinates
       """
    coordinates = []  # create an empty array

    if prob_groups[0] == 1 or prob_in_group[0] == 1:  # if there are no sias, return an empty array
        return None
    else:
        random_values = np.random.rand(len(ground_coords))

        for i in range(len(ground_coords)):  # for every original coordinate
            rand = random_values[i]  # use one of the random numbers between 0 and 1
            prob = 0  # and set the probability value to 0

            # produce random radii in bulk
            nr_new_radii = 0.5 * len(prob_groups) * (len(prob_groups) + 1)
            radii = scipy.stats.truncnorm.rvs(-v.cutoff_gaussians, v.cutoff_gaussians, loc=mu_radius_group,
                                              scale=sig_radius_group,
                                              size=int(nr_new_radii))
            heights = np.random.uniform(lower_glycan_height, upper_glycan_height, size=int(nr_new_radii))
            r_index = 0
            h_index = 0

            for j in range(len(prob_groups)):
                if rand < prob + prob_groups[j]:
                    # If the random number is within the probability for 0, 1, 2, ... groups generate j groups
                    # Define the spherical grid
                    theta = np.linspace(0, 0.5 * np.pi, 100)  # Polar angle (latitude), only upper half sphere
                    phi = np.linspace(0, 2 * np.pi, 400)  # Azimuthal angle (longitude)
                    Theta, Phi = np.meshgrid(theta, phi)

                    mask = np.ones_like(Theta,
                                        dtype=bool)  # create a mask that saves which points on the sphere are still available

                    for _ in np.arange(j):
                        # for each group (glycan), determine a radius
                        r = radii[r_index]
                        r_index += 1
                        h_index += 1

                        # generate added angles/points in bulk for each group (glycan)
                        # -> radius of sias to protein should be the same for all of one protein
                        nr_new_angles = len(prob_in_group)
                        sias_dis = scipy.stats.truncnorm.rvs(-v.cutoff_gaussians, v.cutoff_gaussians, loc=mu_sias,
                                                             scale=sigma_sias, size=int(nr_new_angles))

                        theta_unrotated = np.arccos((2 * r ** 2 - sias_dis ** 2) / (2 * r ** 2))
                        phi_unrotated = np.random.uniform(0, 2 * np.pi, int(nr_new_angles))
                        points_unrotated = np.transpose(spherical_to_cartesian(r, theta_unrotated, phi_unrotated))
                        point_index = 0

                        # for each group (glycan), generate a random number to determine the number of sialin acids on it
                        rand2 = np.random.rand()
                        prob2 = 0

                        for k in range(len(prob_in_group)):
                            if rand2 < prob2 + prob_in_group[k]:
                                # generate k coordinates (sias) for that group (glycan)
                                current_rotating_angle = 0
                                for pt in range(k):
                                    if pt == 0:  # create random (possible) angles for the first point and update the mask of possible angles
                                        random_theta, random_phi = random_spherical_point_from_meshgrid(Theta, Phi,
                                                                                                        mask)
                                        mask = generate_new_spherical_mask(Theta, Phi, mask,
                                                                           spherical_to_cartesian(1, random_theta,
                                                                                                  random_phi),
                                                                           forbidden_angle)
                                        # convert the point to cartesian coordinates and save the theta as rotation angle
                                        new_x, new_y, new_z = spherical_to_cartesian(r, random_theta, random_phi)
                                        current_rotating_angle_theta = random_theta
                                        current_rotating_angle_phi = random_phi
                                    else:
                                        # use one of the points that were used around the center and rotate them accordingly
                                        new_x, new_y, new_z = rotate_point_theta_phi(points_unrotated[point_index],
                                                                                     current_rotating_angle_theta,
                                                                                     current_rotating_angle_phi)
                                        point_index += 1
                                        _, current_rotating_angle_theta, current_rotating_angle_phi = cartesian_to_spherical(
                                            new_x, new_y, new_z)
                                        mask = generate_new_spherical_mask(Theta, Phi, mask, spherical_to_cartesian(1,
                                                                                                                    current_rotating_angle_theta,
                                                                                                                    current_rotating_angle_phi),
                                                                           forbidden_angle)

                                    # translate the coordinates from ground coordinate center to original coordinate system
                                    new_point = [ground_coords[i, 0] + new_x, ground_coords[i, 1] + new_y,
                                                 heights[h_index] + new_z]
                                    # if the random point lies within the FOV, take it as your new point
                                    if 0 <= (ground_coords[i, 0] + new_x) <= x_size and 0 <= (
                                            ground_coords[i, 1] + new_y) <= y_size and heights[h_index] + new_z > 0:
                                        coordinates.append(new_point)
                                break

                            else:
                                prob2 += prob_in_group[k]
                    break
                else:
                    prob += prob_groups[j]

    return np.array(coordinates)
