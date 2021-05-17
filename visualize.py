# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 12:31:52 2021

@author: hboog
"""

import dist_calc as dc
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import progressbar
from textwrap import wrap

# -----------------------------------------------------------------------------


def read_out_cells(receptor_directory, nucleus_directory, resolution_x,
                   resolution_y, resolution_z, output_filepath,
                   figtitle=0):
    receptor_filelist = dc.COM_read_directory_files(receptor_directory,
                                                    ".comred")
    nucleus_filelist = dc.COM_read_directory_files(nucleus_directory,
                                                   ".comred")
    with PdfPages(output_filepath) as pdf:
        print("Writing plots...")
        with progressbar.ProgressBar(max_value=len(receptor_filelist)) as bar:
            count = 1
            for rec, nuc in zip(receptor_filelist, nucleus_filelist):
                # Do I need the real resolutions here? Yes bcs of different
                # possible resolutions in x, y, and z
                recvals = dc.read_comred(rec)
                nucvals = dc.read_comred(nuc)
                x, y, z, virtualmass = [], [], [], []
                for line in recvals:
                    x.append((float(line[5]) - float(nucvals[0][5])) *
                             float(resolution_x))
                    y.append((float(line[6]) - float(nucvals[0][6])) *
                             float(resolution_y))
                    z.append((float(line[7]) - float(nucvals[0][7])) *
                             float(resolution_z))
                    # division by 100 here to make less cluttered plots
                    virtualmass.append(float(line[1])/100)
                rng = np.random.RandomState(0)
                colors = rng.rand(len(x))
                fig = plt.figure(figsize=(8, 20))
                if not figtitle:
                    figtitle = ("\n".join(
                        wrap("Visualization of file {}".format(rec), 60)))
                plt.suptitle(figtitle)
                ax1 = plt.subplot(311)
                ax1.set_title("X-Y Plane")
                ax1.set_xlabel("x coordinate [µm]")
                ax1.set_ylabel("y coordinate [µm]")
                plt.scatter(x, y, c=colors, s=virtualmass, alpha=0.3,
                            cmap='viridis')
                ax2 = plt.subplot(312)
                ax2.set_title("X-Z Plane")
                ax2.set_xlabel("x coordinate [µm]")
                ax2.set_ylabel("z coordinate [µm]")
                plt.scatter(x, z, c=colors, s=virtualmass, alpha=0.3,
                            cmap='viridis')
                ax3 = plt.subplot(313)
                ax3.set_title("Y-Z Plane")
                ax3.set_xlabel("y coordinate [µm]")
                ax3.set_ylabel("z coordinate [µm]")
                plt.scatter(y, z, c=colors, s=virtualmass, alpha=0.3,
                            cmap='viridis')
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig()
                plt.close(fig)
                bar.update(count)
                count += 1
    print("Finished writing plots")


# -----------------------------------------------------------------------------

def vis_args():
    description = ("ComRed 0.3 - Last update: 2021-05-14 - Script for "
                   "visualizing comreds calculations. Plots 3 plots per "
                   "analyzed cell in x-y, x-z, and y-z planes. The size of "
                   "the scatter plots corresponds to their relative virtual "
                   "mass (the product of size and mean intensity). "
                   "As such larger circles in the plot correspond to a higher"
                   " weight in the final COM calculation. The circles color "
                   "is randomly pulled from a lookup table and is used to "
                   "differ between overlapping circles.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-o", "--output_directory", metavar=("Output "
                                                             "directory"),
                        help=("Set the output directory, where results"
                              " will be saved, required parameter"),
                        required=True)
    parser.add_argument("-f", "--output_filename", metavar=("Output "
                                                            "filename"),
                        help=("Set a filename scheme for output files"))
    parser.add_argument("-r", "--receptors", metavar="Receptor directory",
                        help=("Directory containing the receptor files, "
                              "all .comred files in the directory will be "
                              "read in, so use a different directory for"
                              " the reference files, required parameter"),
                        required=True)
    parser.add_argument("-n", "--reference", metavar=("Nucleus/Reference"
                                                      " directory"),
                        help=("Directory containing the reference (usually"
                              " nucleus) files, all .comred files in the "
                              "directory will be read in, so use a "
                              "different directory for the receptor"
                              " files, required parameter"), required=True)
    parser.add_argument("-x", "--x_resolution", metavar=("X-resolution"),
                        help=("x Resolution in micrometer per pixel"
                              ", required parameter"),
                        required=True, type=float)
    parser.add_argument("-y", "--y_resolution", metavar=("Y-resolution"),
                        help=("y Resolution in micrometer per pixel"
                              ", required parameter"),
                        required=True, type=float)
    parser.add_argument("-z", "--z_resolution", metavar=("Z-resolution"),
                        help=("Z-Resolution in micrometer per pixel"
                              ", required parameter"),
                        required=True, type=float)
# =============================================================================
#     parser.add_argument('-t', "--title", metavar="Title",
#                         help=("Titles for the plots. Wrap titles in ' ' to "
#                               "allow for spaces in titles"), type=str,
#                         required=False)
# =============================================================================
    args = parser.parse_args()
    output_filename = args.output_filename + ".pdf"
    output_filepath = os.path.join(args.output_directory, output_filename)
    read_out_cells(args.receptors, args.reference, args.x_resolution,
                   args.y_resolution, args.z_resolution, output_filepath)


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    vis_args()
