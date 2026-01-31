# file of important functions when using disperse; note many require a specific run output of .RaDecZ.a.NDskl
import numpy as np
from scipy.spatial import KDTree

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D

import plotly.graph_objects as go

from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

# ------------------------------------------------------------------------------------------------------------------------------- #

def filaments_extract(file_name):
    """Extracts the ra, dec, and z values outlining the filaments from a skeleton file output

    file_name: string of a file which is in .RaDecZ.a.NDskl format

    Returns: array of the separate filaments, which contain each an array of the ra, dec, and z values outlining their shape
    """
# open resulting skeleton file
    filaments = ([])

    with open(file_name) as f:
        for line in f:
            # go until reach the filaments header
            if line.strip().startswith("[FILAMENTS]"):
                break

        # the line right underneath the filaments header is the number of data points
        n_fil = int(next(f)) 

        for i in range(n_fil):
            # each section has a header describing the filament ID, critical point, and number of points
            header = next(f).strip()

            # moves the read to the next line after that
            if not header:
                continue

            parts = header.split()
            # it should be of a length 3 (ra, dec, z)
            if len(parts) != 3:
                continue
            filament_id, cp_id, n_points = map(int, parts)

            coords = ([])
            for i in range(n_points):
                values = list(map(float, next(f).split()))
                coords.append(values)

            filaments.append(np.array(coords))

    return filaments

# ------------------------------------------------------------------------------------------------------------------------------- #

def ra_dec_z_convert(ra_values, dec_values, z_values, reference = None):
    """Converts ra, dec, and z values into cartesian coords in units of Mpc

    ra_values: array of right ascension values in degrees
    dec_values: array of declination values in degrees
    z_values: array of redshift values
    reference: list [ra, dec, z] of the center of the field, automatically set to zero
    
    Returns: distances in x, y, and z from reference point in Mpc"""

    if reference is None:
        ra_ref, dec_ref, z_ref = [0, 0, 0]
    else:
        ra_ref, dec_ref, z_ref = reference

    # convert from degrees to radians
    ra = np.radians(ra_values)
    dec = np.radians(dec_values)
    ra0 = np.radians(ra_ref)
    dec0 = np.radians(dec_ref)

    # using comoving distance for converting redshift values
    cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)

    r = cosmo.comoving_distance(z_values).to(u.Mpc).value
    r0 = cosmo.comoving_distance(z_ref).to(u.Mpc).value

    # this is just trig; using r to convert
    x = r * np.cos(dec) * np.cos(ra)
    y = r * np.cos(dec) * np.sin(ra)
    z = r * np.sin(dec)

    x0 = r0 * np.cos(dec0) * np.cos(ra0)
    y0 = r0 * np.cos(dec0) * np.sin(ra0)
    z0 = r0 * np.sin(dec0)

    # make it in terms of reference point, does nothing if set to zero
    dist_x = x - x0
    dist_y = y - y0
    dist_z = z - z0
    
    return dist_x, dist_y, dist_z

# ------------------------------------------------------------------------------------------------------------------------------- #

def plot_filaments2D(filename, center_reference = False, ax = None):
    """Takes a file and condition if a non-zero reference is being used, and outputs plot of filaments

    filename: string of a file which is in .RaDecZ.a.NDskl format
    center_reference: bool of whether the reference is zero or the center of the field (typical use in astronomy)
    ax: allows ability to plot as subplots, automatically set to None

    Returns: 2D plot of filaments
    """
    filaments = filaments_extract(filename)

    if center_reference:
        all_ra = np.concatenate([f[:, 0] for f in filaments])
        all_dec = np.concatenate([f[:, 1] for f in filaments])
        all_z  = np.concatenate([f[:, 2] for f in filaments])
        reference = (np.median(all_ra), np.median(all_dec), np.median(all_z))
    else:
        reference = None

    # if axes not specified, just outputs the regular one plot
    if ax is None:
        fig, ax = plt.subplots(figsize = (10, 8))

    filament_distances = []

    for filament in filaments:
        ra = filament[:, 0]
        dec = filament[:, 1]
        zr = filament[:, 2]

        x, y, z = ra_dec_z_convert(ra, dec, zr, reference = reference)
        filament_distances.append(np.median(z))

    normalization = colors.Normalize(vmin = min(filament_distances), vmax = max(filament_distances))
    colormap = cm.viridis

    # plotting each filament
    for filament, distance in zip(filaments, filament_distances):
        ra = filament[:, 0]
        dec = filament[:, 1]
        zr = filament[:, 2]

        x, y, z = ra_dec_z_convert(ra, dec, zr, reference)

        color = colormap(normalization(distance))   # color based on the distance
        ax.plot(x, y, linewidth = 0.4, color = color)

    mappable = cm.ScalarMappable(norm = normalization, cmap = colormap)
    plt.colorbar(mappable, ax = ax, label = "median filament z (Mpc from reference)")

    ax.grid(visible = True)
    ax.set_aspect('equal')
    ax.set_xlabel("x (Mpc from reference)")
    ax.set_ylabel("y (Mpc from reference)")    

    return ax

# ------------------------------------------------------------------------------------------------------------------------------- #

def plot_filaments3D(filename, center_reference = False):
    """Takes file name and whether a non-zero reference is being used, and outputs a 3D plot of filaments

    filename: string of a file which is in .RaDecZ.a.NDskl format
    center_reference: bool of whether the reference is zero or the center of the field

    Returns: 3D plot of the filaments
    """
    filaments = filaments_extract(filename)

    if center_reference:
        all_ra = np.concatenate([f[:, 0] for f in filaments])
        all_dec = np.concatenate([f[:, 1] for f in filaments])
        all_z  = np.concatenate([f[:, 2] for f in filaments])
        reference = (np.median(all_ra), np.median(all_dec), np.median(all_z))
    else:
        reference = None

    filament_distances = []

    for filament in filaments:
        ra = filament[:, 0]
        dec = filament[:, 1]
        zr = filament[:, 2]

        x, y, z = ra_dec_z_convert(ra, dec, zr, reference = reference)
        filament_distances.append(np.median(z))

    normalization = colors.Normalize(vmin = min(filament_distances), vmax = max(filament_distances))
    colormap = cm.viridis

    fig = plt.figure(figsize = (12, 8))
    ax = fig.add_subplot(projection = "3d")

    for filament, distance in zip(filaments, filament_distances):
        ra = filament[:, 0]
        dec = filament[:, 1]
        zr = filament[:, 2]

        x, y, z = ra_dec_z_convert(ra, dec, zr, reference)

        color = colormap(normalization(distance))

        ax.plot(x, y, z, linewidth = 0.4, color = color)

    ax.grid(visible = True)
    ax.set_aspect("equal")

    ax.set_xlabel("x (Mpc from reference)")
    ax.set_ylabel("y (Mpc from reference)")
    ax.set_zlabel("z (Mpc from reference)")

    return ax

# ------------------------------------------------------------------------------------------------------------------------------- #

def plot_filaments_interactive(filename, center_reference = False):
    """Plot filaments in 3D with interactive environment"""
    filaments = filaments_extract(filename)

    if center_reference:
        all_ra = np.concatenate([f[:, 0] for f in filaments])
        all_dec = np.concatenate([f[:, 1] for f in filaments])
        all_z  = np.concatenate([f[:, 2] for f in filaments])
        reference = (np.median(all_ra), np.median(all_dec), np.median(all_z))
    else:
        reference = None

    fig = go.Figure()

    for filament in filaments:
        ra = filament[:, 0]
        dec = filament[:, 1]
        zr = filament[:, 2]

        x, y, z = ra_dec_z_convert(ra, dec, zr, reference)

        fig.add_trace(go.Scatter3d(x = x, y = y, z = z, 
                                   mode = "lines", 
                                   line = dict(color = "slateblue", width = 2), 
                                   hoverinfo = "none"))

    fig.update_layout(
        scene = dict(
            xaxis_title = "x (Mpc from reference)", 
            yaxis_title = "y (Mpc from reference)", 
            zaxis_title = "z (Mpc from reference)", 
            aspectmode = "data"
        ),
        margin = dict(l = 0, r = 0, b = 0, t = 0)
    )

    return fig

# ------------------------------------------------------------------------------------------------------------------------------- #

def filaments_convert(filename, center_reference = False):
    """Converts the filaments from ra, dec, and z into cartesian coords in Mpc while keeping each filament separate

    filename: the file
    center_reference: bool, optional to set the reference as the center or just [0, 0, 0]

    Returns: the filaments but in x, y, z coords
    """
    filaments = filaments_extract(filename)

    if center_reference:
        all_ra = np.concatenate([f[:, 0] for f in filaments])
        all_dec = np.concatenate([f[:, 1] for f in filaments])
        all_z  = np.concatenate([f[:, 2] for f in filaments])
        reference = (np.median(all_ra), np.median(all_dec), np.median(all_z))
    else:
        reference = None

    filaments_coords = []

    for filament in filaments:
        ra = filament[:, 0]
        dec = filament[:, 1]
        zr = filament[:, 2]

        x, y, z = ra_dec_z_convert(ra, dec, zr, reference = reference)

        filament = np.column_stack((x, y, z))
        filaments_coords.append(filament)

    return filaments_coords

# ------------------------------------------------------------------------------------------------------------------------------- #

def filament_tree_build(filaments):
    """Builds a KDtree for the filaments using each galaxy as a point

    filaments: array containing each filament

    Returns: KDtree, the points that were used, and the information of the indexes
    """
    points = []
    index_map = []

    for fil_id, fil in enumerate(filaments):
        for v_id, point in enumerate(fil):
            points.append(point)
            index_map.append((fil_id, v_id))

    points = np.array(points)

    # build tree
    tree = KDTree(points)

    return tree, points, index_map
    
# ------------------------------------------------------------------------------------------------------------------------------- #

def point_to_segment(point, a, b):
    """Minimizes a vector projection

    point: an integer or float of the point we want to find the distance of to the segment
    a: first bound of the segment
    b: second bound of the segment

    returns: minimum distance from the point to the segment
    """
    ab = b - a 

    t = np.dot(point - a, ab) / np.dot(ab, ab)

    # need to make sure that it stays on the segment
    t = np.clip(t, 0.0, 1.0)

    closest_t = a + t * ab
    distance = np.linalg.norm(point - closest_t)

    return distance

# ------------------------------------------------------------------------------------------------------------------------------- #

def distances_to_filaments(galaxy_points, filaments):
    """Uses the KDtree to find the nearest distance of the set of galaxies to a filament

    filaments: array containing each filament

    Returns: array of all the minimized distances    
    """
    # build the tree
    tree, fil_points, index_map = filament_tree_build(filaments)

    distances = np.zeros(len(galaxy_points))

    # loop over all points
    for i, p in enumerate(galaxy_points):
        # find the nearest filament point
        distance, idx = tree.query(p)
        fil_id, v_id = index_map[idx]

        filament = filaments[fil_id]

        # examine closeby filaments
        candidates = []

        if v_id > 0:
            candidates.append((filament[v_id - 1], filament[v_id]))
        if v_id < len(filament) - 1:
            candidates.append((filament[v_id], filament[v_id + 1]))

        # then compute the distances 
        best_distance = min(point_to_segment(p, a, b) for a, b in candidates)
        distances[i] = best_distance

    return distances

# ------------------------------------------------------------------------------------------------------------------------------- #

def cumul_dist(distances):
    """Finds the cumulative distribution of the inputted distances and the point where it reaches 50%
    
    distance: array of distances from running the distances_to_filaments function
    
    Returns: distances (sorted by size), cumulative distribution, and 50% mark"""
    sorted_dist = np.sort(distances)
    cdf = np.arange(1, len(sorted_dist) + 1) / len(sorted_dist)

    halfway_point = max(sorted_dist[cdf <= 0.5])

    return sorted_dist, cdf, halfway_point

