#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import progressbar


# -----------------------------------------------------------------------------

class dist_calc:

    """COMRed Version 0.1.1 - Dist_Calc"""

    def __init__(self):
        '''Constructor'''
        # self.filename = filename
        1

# -----------------------------------------------------------------------------

    def read_comred(filename):
        with open(filename, "r") as comred_file:
            comvals = []
            iter_comred_file = iter(comred_file)
            next(iter_comred_file)
            for line in iter_comred_file:
                line = line.strip("\n")
                line = line.split("\t")
                virtualmass = float(line[4]) * float(line[5])
                masspos_x = float(line[1]) * virtualmass
                masspos_y = float(line[2]) * virtualmass
                masspos_z = float(line[3]) * virtualmass
                values = [line[0], virtualmass, masspos_x, masspos_y,
                          masspos_z, line[1], line[2], line[3], line[4],
                          line[5]]
                comvals.append(values)
            return comvals

# -----------------------------------------------------------------------------

    def calculate_com_from_comvals(comvals):
        sum_virtualmass, sum_masspos_x, sum_masspos_y = 0, 0, 0
        sum_masspos_z = 0
        for element in comvals:
            sum_virtualmass = sum_virtualmass + element[1]
            sum_masspos_x = sum_masspos_x + element[2]
            sum_masspos_y = sum_masspos_y + element[3]
            sum_masspos_z = sum_masspos_z + element[4]
        center_of_mass_pic_x = sum_masspos_x / sum_virtualmass
        center_of_mass_pic_y = sum_masspos_y / sum_virtualmass
        center_of_mass_pic_z = sum_masspos_z / sum_virtualmass
        center_of_mass_pic = [center_of_mass_pic_x, center_of_mass_pic_y,
                              center_of_mass_pic_z]
        return center_of_mass_pic

# -----------------------------------------------------------------------------

    def calculate_spread_from_comvals(comvals, resolution_x, resolution_y,
                                      resolution_z):
        # Takes all accumulations and compares them with the center of mass of
        # the same distribution (weighted by virtual mass)
        # giving a mean distance of receptors
        # return center_of_mass_pic -> put distance std here???
        com_pic = dist_calc.calculate_com_from_comvals(comvals)
        distance_list = []
        for element in comvals:
            x = dist_calc.Conversion_Pixel_micrometer(
                (float(element[5]) - com_pic[0]), resolution_x)
            y = dist_calc.Conversion_Pixel_micrometer(
                (float(element[6]) - com_pic[1]), resolution_y)
            z = dist_calc.Conversion_Pixel_micrometer(
                (float(element[7]) - com_pic[2]), resolution_z)
            distance = (x**2 + y**2 + z**2)**0.5
            distance_list.append(distance)
        average_distance = numpy.mean(distance_list)
        return average_distance

# -----------------------------------------------------------------------------

    def COM_comparison_with_nucleus(receptor_file, nucleus_file, resolution_x,
                                    resolution_y, resolution_z):
        # ====================================================================
        # mppxy = 0.0763298 #Resolution in xy direction in micrometer
        # mppz = 0.6337175 #Resolution in z direction in micrometer
        # ====================================================================
        receptor_COM = dist_calc.calculate_com_from_comvals(
                dist_calc.read_comred(receptor_file))
        nucleus_COM = dist_calc.calculate_com_from_comvals(
                dist_calc.read_comred(nucleus_file))
        difference_x = dist_calc.Conversion_Pixel_micrometer(
                        (receptor_COM[0] - nucleus_COM[0]), resolution_x)
        difference_y = dist_calc.Conversion_Pixel_micrometer(
                        (receptor_COM[1] - nucleus_COM[1]), resolution_y)
        difference_z = dist_calc.Conversion_Pixel_micrometer(
                        (receptor_COM[2] - nucleus_COM[2]), resolution_z)
        Arith_mean = (difference_x**2 + difference_y**2 + difference_z**2)**0.5
        difference = [difference_x, difference_y, difference_z, Arith_mean]
        return difference

# -----------------------------------------------------------------------------

# Reads all files from a directory with the .csv expansion, returns list of
# files
    def COM_read_directory_files(directory_name, ending):
        filelist = []
        for filename in os.listdir(directory_name):
            if filename.endswith(ending):
                filelist.append(os.path.join(directory_name, filename))
        return filelist

# -----------------------------------------------------------------------------

# Calculates all COM difference values for files in 2 directories, make sure
# the number of files are equal!!!
    def Calculate_COMs_from_directories(receptor_directory, nucleus_directory,
                                        resolution_x, resolution_y,
                                        resolution_z, verbose=0):
        receptor_filelist = dist_calc.COM_read_directory_files(
                receptor_directory, ".comred")
        nuclei_filelist = dist_calc.COM_read_directory_files(
                nucleus_directory, ".comred")
        COM_value_list = []
        distance_value_list = []
        rn_filelist = []
        filenumber = 0
        for filenumber in range(0, len(receptor_filelist)):
            new_entry = [receptor_filelist[filenumber],
                         nuclei_filelist[filenumber]]
            rn_filelist.append(new_entry)
        with progressbar.ProgressBar(max_value=len(rn_filelist)) as bar:
            count = 0
            for filepair in rn_filelist:
                if verbose:
                    print("Calculating filepair: {} and {}".format(
                        filepair[0], filepair[1]))
                COM_value_list.append(dist_calc.COM_comparison_with_nucleus(
                        filepair[0], filepair[1], resolution_x, resolution_y,
                        resolution_z))
                distance_value_list.append(
                    dist_calc.calculate_spread_from_comvals(
                        dist_calc.read_comred(filepair[0]), resolution_x,
                        resolution_y, resolution_z))
                count += 1
                bar.update(count)
        return COM_value_list, distance_value_list

# -----------------------------------------------------------------------------

# Calculates mean and std of files returned from the
# Calculate_COMs_from_directories funcion
    def COMs_statistics_csv(value_list):
        # for the list of values coming out of Calculate_COM_from_directories
        values_x, values_y, values_z, Arith_mean = [], [], [], []
        for value in value_list:
            values_x.append(value[0])
            values_y.append(value[1])
            values_z.append(value[2])
            Arith_mean.append(value[3])
        data_x = [numpy.mean(values_x), numpy.std(values_x)]
        data_y = [numpy.mean(values_y), numpy.std(values_y)]
        data_z = [numpy.mean(values_z), numpy.std(values_z)]
        data_Arith_mean = [numpy.mean(Arith_mean), numpy.std(Arith_mean)]
        end_data = [data_x, data_y, data_z, data_Arith_mean]
        return end_data

# -----------------------------------------------------------------------------

# Calculates mean and std of files returned from the
# Calculate_COMs_from_directories funcion
    def COMs_statistics_distance(value_list):
        # for the list of values coming out of Calculate_COM_from_directories
        values_distance = []
        for value in value_list:
            values_distance.append(value)
        end_data = [numpy.mean(values_distance), numpy.std(values_distance)]
        return end_data

# -----------------------------------------------------------------------------

# Convertes Pixel values in micrometer
    def Conversion_Pixel_micrometer(number_in_pixels, micrometer_per_pixel):
        Micrometer = float(number_in_pixels) * float(micrometer_per_pixel)
        return Micrometer

# -----------------------------------------------------------------------------

# Makes a file in output_directory with output_filename.txt, uses
# Calculate_COMs_from_directories so two file directories need to be given
    def make_file_with_COMs(output_filename, output_directory,
                            receptor_directory, nucleus_directory,
                            resolution_x, resolution_y, resolution_z,
                            verbose=0):
        print("Constructing whole cell center of mass file...")
        output_filepath = "".join([output_directory, "\\", output_filename,
                                   ".txt"])
        values_com, value_distance = dist_calc.Calculate_COMs_from_directories(
                receptor_directory, nucleus_directory, resolution_x,
                resolution_y, resolution_z, verbose)
        end_values = dist_calc.COMs_statistics_csv(values_com)
        end_vals_d = dist_calc.COMs_statistics_distance(value_distance)
        with open(output_filepath, "w") as output_file:
            output_file.write("#All data in micrometer.\n\n"
                              "#Calculated means:\n")
            alph = ["x", "y", "z"]
            for count, element in enumerate(end_values):
                if count <= 2:
                    mean_r = ("#Mean_Difference_{}:\t{}\n#"
                              "Standarddeviation_{}:\t{}\n").format(
                                  alph[count], element[0], alph[count],
                                  element[1])
                    output_file.write(mean_r)
                else:
                    mean_r = ("#Mean_Arithmetic_difference:\t{}\n"
                              "#Standarddeviation_Arithmetic_difference:"
                              "\t{}\n").format(element[0], element[1])
                    output_file.write(mean_r)
            write_distance = ("#Mean spread distance:\t{}\n"
                              "#Standarddeviation:\t{}\n\n").format(
                                  end_vals_d[0], end_vals_d[1])
            output_file.write(write_distance)
            heading = ("#Difference_x\tDifference_y"
                       "\tDifference_z\tArithmetic_difference\t"
                       "Mean_spread_distance\n")
            output_file.write(heading)
            count = 0
            for line in zip(values_com, value_distance):
                write_val = "{}\t{}\t{}\t{}\t{}\n".format(line[0][0],
                                                          line[0][1],
                                                          line[0][2],
                                                          line[0][3],
                                                          line[1])
                output_file.write(write_val)
        print("File successfully written")

# -----------------------------------------------------------------------------

    def make_totals_from_columnvalues(columnvalues_list):
        total_difference_x, total_difference_y, total_difference_z = [], [], []
        total_object_volume, total_Arith_mean = [], []
        for columnvalues in columnvalues_list:
            for element in columnvalues[0]:
                total_difference_x.append(element)
            for element in columnvalues[1]:
                total_difference_y.append(element)
            for element in columnvalues[2]:
                total_difference_z.append(element)
            for element in columnvalues[3]:
                total_object_volume.append(element)
            for element in columnvalues[4]:
                total_Arith_mean.append(element)
        total_columnvals = [total_difference_x, total_difference_y,
                            total_difference_z, total_object_volume,
                            total_Arith_mean]
        return total_columnvals

# -----------------------------------------------------------------------------

    def make_histogram(data_list, xlabel, title, color):
        bins = numpy.linspace(min(data_list), max(data_list), 20)
        plt.hist(data_list, density=0, bins=bins,
                 color=color, edgecolor="none")
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel('Count')
        plt.grid(axis="y", color="0.9", linestyle="-")

# -----------------------------------------------------------------------------

    def make_boxplot(data_list, ylabel, title, coloris):
        box = plt.boxplot(data_list, 0, patch_artist=True,
                          flierprops={"markersize": 4,
                                      "markeredgecolor": "None",
                                      "marker": ".",
                                      "markerfacecolor": coloris})
        for element in ["whiskers", 'fliers', 'caps']:
            plt.setp(box[element], color=coloris)
        for element in ["means", "medians"]:
            plt.setp(box[element], color="black")
        plt.setp(box["boxes"], color="none")
        for patch, color in zip(box['boxes'], coloris):
            patch.set_facecolor(coloris)
        plt.title(title)
        plt.ylabel(ylabel)

# -----------------------------------------------------------------------------

    def read_output_data_whole_cell_COM(filepath):
        with open(filepath, "r") as file:
            end_values = []
            for i in range(16):
                next(file)
            for line in file:
                values = []
                line = line.strip()
                if line:
                    line = line.split("\t")
                    for element in line:
                        values.append(float(element))
                end_values.append(values)
            return end_values

# -----------------------------------------------------------------------------

    def sort_whole_cell_COM_values(filename):
        values = dist_calc.read_output_data_whole_cell_COM(filename)
        difference_x, difference_y, difference_z, Arith_mean = [], [], [], []
        for elements in values:
            if elements:
                difference_x.append(elements[0])
                difference_y.append(elements[1])
                difference_z.append(elements[2])
                Arith_mean.append(elements[3])
        columnvalues = [difference_x, difference_y, difference_z, Arith_mean]
        return columnvalues

# -----------------------------------------------------------------------------

    def sort_distance_COM_values(filename):
        values = dist_calc.read_output_data_whole_cell_COM(filename)
        distances = []
        for elements in values:
            if elements:
                distances.append(elements[4])
        return distances

# -----------------------------------------------------------------------------

    def create_plots_whole_cell_COMs(output_filename, output_directory,
                                     input_filepath, color):
        print("Creating Whole Cell Center of Mass plots...")
        output_filepath = "".join([output_directory, "\\", output_filename,
                                  ".pdf"])
        with PdfPages(output_filepath) as pdf:
            with progressbar.ProgressBar(max_value=10) as bar:
                total_values = dist_calc.sort_whole_cell_COM_values(
                    input_filepath)
                fig = plt.figure(figsize=(8, 5))
                figtitle = ("Histograms of total values from "
                            + str(input_filepath))
                plt.suptitle(figtitle)
                plt.subplot(221)
                dist_calc.make_histogram(total_values[0], "Difference [µm]",
                                         "x-Direction", color)
                bar.update(1)
                plt.subplot(222)
                dist_calc.make_histogram(total_values[1], "Difference [µm]",
                                         "y-Direction", color)
                bar.update(2)
                plt.subplot(223)
                dist_calc.make_histogram(total_values[2], "Difference [µm]",
                                         "z-Direction", color)
                bar.update(3)
                plt.subplot(224)
                dist_calc.make_histogram(total_values[3], "Difference [µm]",
                                         "Arithmetic mean", color)
                bar.update(4)
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig()
                plt.close(fig)
                fig = plt.figure(figsize=(8, 5))
                figtitle = ("Boxplots of total values from "
                            + str(input_filepath))
                plt.suptitle(figtitle)
                plt.subplot(221)
                dist_calc.make_boxplot(total_values[0], "Difference [µm]",
                                       "x-Direction", color)
                bar.update(5)
                plt.subplot(222)
                dist_calc.make_boxplot(total_values[1], "Difference [µm]",
                                       "y-Direction", color)
                bar.update(6)
                plt.subplot(223)
                dist_calc.make_boxplot(total_values[2], "Difference [µm]",
                                       "z-Direction", color)
                bar.update(7)
                plt.subplot(224)
                dist_calc.make_boxplot(total_values[3], "Difference [µm]",
                                       "Arithmetic mean", color)
                bar.update(8)
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig()
                plt.close(fig)
                distances = dist_calc.sort_distance_COM_values(input_filepath)
                fig = plt.figure(figsize=(5, 8))
                figtitle = ("Mean distance of accumulations in "
                            + str(input_filepath))
                plt.suptitle(figtitle)
                plt.subplot(211)
                dist_calc.make_histogram(distances, "Difference [µm]",
                                         ("Mean distance of receptors to"
                                          " their COM"), color)
                bar.update(9)
                plt.subplot(212)
                dist_calc.make_boxplot(distances, "Difference [µm]",
                                       ("Mean distance of receptors to"
                                        " their COM"), color)
                bar.update(10)
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig()
                plt.close(fig)
        print("Finished generating plots")

# -----------------------------------------------------------------------------

    def generate_COM_data_from_directories():
        description = ("ComRed 0.2 - Last update: 2021-02-03 - Script for "
                       "calculating the Center of Mass distances between COM "
                       "of the receptors and a reference, as well as the "
                       "mean distance of receptors to their Center of Mass")
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
        parser.add_argument("-c", "--color", metavar=("Plot color"),
                            help=("Color of histograms and boxplots"))
        parser.add_argument("-v", "--verbose", action=("store_true"),
                            help=("More details while running the program"))
        args = parser.parse_args()
        dist_calc.make_file_with_COMs(args.output_filename,
                                      args.output_directory,
                                      args.receptors, args.reference,
                                      args.x_resolution, args.y_resolution,
                                      args.z_resolution, args.verbose)
        input_filepath_whole_cell = (args.output_directory + "\\"
                                     + args.output_filename + ".txt")
        dist_calc.create_plots_whole_cell_COMs(args.output_filename,
                                               args.output_directory,
                                               input_filepath_whole_cell,
                                               args.color)


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    dist_calc.generate_COM_data_from_directories()
