#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scipy.stats as sps
import math
from dist_calc import dist_calc as dc
import argparse
import progressbar


# -----------------------------------------------------------------------------

def stars(p):
    if p < 0.0001:
        return "****"
    elif (p < 0.001):
        return "***"
    elif (p < 0.01):
        return "**"
    elif (p < 0.05):
        return "*"
    else:
        return "-"


# -----------------------------------------------------------------------------

def mwu_effectsize(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_std = np.std(x)
    y_std = np.std(y)
    n_x = len(x)
    n_y = len(y)
    u, pval = sps.mannwhitneyu(x, y, use_continuity=True,
                               alternative="two-sided")
    d = math.fabs((x_mean - y_mean)/(math.sqrt(((n_x - 1) * x_std ** 2
                                                + (n_y + 1) * y_std ** 2)
                                               / (n_x + n_y - 2))))
    z = (u - n_x * n_y / 2)/(math.sqrt(n_x * n_y * (n_x + n_y + 1)/12))
    r = math.fabs(z/math.sqrt(n_x + n_y))
    return pval, d, r, x_mean, y_mean, x_std, y_std, n_x, n_y


# -----------------------------------------------------------------------------

def mwu_text(x, y, name_x, name_y):
    pval, d, r, x_mean, y_mean, x_std, y_std, n_x, n_y = mwu_effectsize(x, y)
    textstr = ("p = {}\nCohens d = {}\nPearsons r = {}\n{} mean: {}\n{} "
               "mean: {}\n{} std: {}\n{} std: {}\n{} number data points: {}\n"
               "{} number data points: {}").format(pval, d, r, name_x, x_mean,
                                                   name_y, y_mean, name_x,
                                                   x_std, name_y, y_std,
                                                   name_x, n_x, name_y, n_y)
    return textstr


# -----------------------------------------------------------------------------

def make_space_above(axes, topmargin=1):
    """ increase figure size to make topmargin (in inches) space for
        titles, without changing the axes sizes"""
    fig = axes.flatten()[0].figure
    s = fig.subplotpars
    w, h = fig.get_size_inches()

    figh = h - (1-s.top)*h + topmargin
    fig.subplots_adjust(bottom=s.bottom*h/figh, top=1-topmargin/figh)
    fig.set_figheight(figh)


# -----------------------------------------------------------------------------

def create_comparison_COM_plot(output_filename, output_directory,
                               input_filepath_whole_1, input_filepath_whole_2,
                               datanames):
    print("Create comparison plots")
    output_filepath = "".join([output_directory, "/", output_filename, ".pdf"])
    with PdfPages(output_filepath) as pdf:
        with progressbar.ProgressBar(max_value=5) as bar:
            whole_cell_values = [
                    dc.sort_whole_cell_COM_values(input_filepath_whole_1),
                    dc.sort_whole_cell_COM_values(input_filepath_whole_2)]
            fig = plt.figure(figsize=(12, 20))
            figtitle = ("Comparison of whole cell totals from "
                        "\n{}\n and \n{}").format(input_filepath_whole_1,
                                                  input_filepath_whole_2)
            plt.suptitle(figtitle, fontsize=9)
            plt.subplot(321)
            make_boxplot_comp_bw(whole_cell_values[0][0],
                                 whole_cell_values[1][0], "Difference [µm]",
                                 "x-Direction", datanames)
            bar.update(1)
            plt.subplot(322)
            make_boxplot_comp_bw(whole_cell_values[0][1],
                                 whole_cell_values[1][1], "Difference [µm]",
                                 "y-Direction", datanames)
            bar.update(2)
            plt.subplot(323)
            make_boxplot_comp_bw(whole_cell_values[0][2],
                                 whole_cell_values[1][2], "Difference [µm]",
                                 "z-Direction", datanames)
            bar.update(3)
            plt.subplot(324)
            make_boxplot_comp_bw(whole_cell_values[0][3],
                                 whole_cell_values[1][3], "Difference [µm]",
                                 "Arithmetic mean", datanames)
            bar.update(4)
            textplot2 = plt.subplot(325)
            mwu_textlist = []
            for r in range(4):
                mwu_textlist.append(mwu_text(whole_cell_values[0][r],
                                             whole_cell_values[1][r],
                                             datanames[0], datanames[1]))
            figtxt = ("x-Direction: \n{}\n\ny-Direction: \n{}\n\nz-Direction:"
                      " \n{}\n\nArithmetic mean: \n{}").format(*mwu_textlist)
            textplot2.text(0, -0.1, figtxt, fontsize=7)
            textplot2.axis("off")
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig()
            plt.close(fig)
            # -----------
            distance_cell_values = [
                    dc.sort_distance_COM_values(input_filepath_whole_1),
                    dc.sort_distance_COM_values(input_filepath_whole_2)]
            fig = plt.figure(figsize=(12, 8))
            figtitle = ("\n".join([
                "Comparison of average accumulation distances in ", "{}",
                " and ", "{}"])).format(input_filepath_whole_1,
                                        input_filepath_whole_2)
            plt.suptitle(figtitle, fontsize=9)
            plt.subplot(121)
            make_boxplot_comp_bw(distance_cell_values[0],
                                 distance_cell_values[1], "Difference [µm]",
                                 "Mean distance of receptors to COM",
                                 datanames)
            bar.update(5)
            textplot2 = plt.subplot(122)
            mwu_textlist = []
            mwu_textlist.append(mwu_text(distance_cell_values[0],
                                         distance_cell_values[1],
                                         datanames[0], datanames[1]))
            figtxt = ("Average distance between receptors and their COM:\n"
                      "{}").format(*mwu_textlist)
            textplot2.text(0.1, -0.1, figtxt, fontsize=7)
            textplot2.axis("off")
            plt.tight_layout(rect=[0, 0.03, 1, 0.88])
            pdf.savefig()
            plt.close(fig)
    print("Done")


# -----------------------------------------------------------------------------

def make_boxplot_comp_bw(data_list_1, data_list_2, ylabel, title, datanames):
    p = mwu_effectsize(data_list_1, data_list_2)[0]
    sp = stars(p)
    box = plt.boxplot([data_list_1, data_list_2], patch_artist=True,
                      flierprops={"markersize": 4, "markeredgecolor": "None",
                                  "marker": "."})
    for elementa, color in zip(box["whiskers"], ["black", "black", "black",
                                                 "black"]):
        elementa.set_color(color)
    for elementi, color in zip(box['fliers'], ["black", "black", "black",
                                               "black"]):
        elementi.set_markerfacecolor(color)
    for elementu, color in zip(box["caps"], ["black", "black", "black",
                                             "black"]):
        elementu.set_color(color)
    for elemento, colori in zip(["means", "medians"], ["grey", "grey"]):
        plt.setp(box[elemento], color=colori)
    plt.setp(box["boxes"], color="black")
    for patch, color in zip(box['boxes'], ["white", "black"]):
        patch.set_facecolor(color)
    y_max = np.max(np.concatenate((data_list_1, data_list_2)))
    y_min = np.min(np.concatenate((data_list_1, data_list_2)))
    plt.annotate("", xy=(1, y_max), xycoords='data', xytext=(2, y_max),
                 textcoords='data', arrowprops=dict(arrowstyle="-",
                 ec='#aaaaaa', connectionstyle="bar,fraction=0.1"))
    plt.text(1.5, y_max + (math.fabs(y_min) + math.fabs(y_max)) * 0.1, sp,
             horizontalalignment='center', verticalalignment='center')
    axes = plt.gca()
    ymax = y_max + (math.fabs(y_min) + math.fabs(y_max)) * 0.2
    axes.set_ylim([y_min, ymax])
    plt.title(title)
    plt.ylabel(ylabel)
    x = [1, 2]
    labels = datanames
    plt.xticks(x, labels, rotation='vertical')


# -----------------------------------------------------------------------------

def invisplot():
    frame1 = plt.gca()
    frame1.axes().get_xaxis().set_visible(False)
    frame1.axes().get_yaxis().set_visible(False)


# -----------------------------------------------------------------------------

def make_comparison_boxplots():
    description = ("ComRed 0.2 - Last update: 2021-03-24 - Script for "
                   "comparing two sets of receptor datasets and produce "
                   "boxplot figures")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-o", "--output_directory", metavar=("Output "
                                                             "directory"),
                        help=("Set the output directory, where results"
                              " will be saved, required parameter"),
                        required=True)
    parser.add_argument("-f", "--output_filename", metavar=("Output filename"),
                        help=("Set a filename scheme for output files"))
    parser.add_argument("-i", "--input", nargs=2,
                        metavar=("Input file 1", "Input file 2"),
                        required=True, default=("File 1", "File 2"),
                        help=("The filepath to the output files from dist_calc"
                              " for the comparison, required parameter,"
                              " divide the two file names with a space and "
                              "put '' around files with space in the names"))
    parser.add_argument("-c", "--caption", nargs=2,
                        metavar=("Left boxplot caption", "Right boxplot "
                                 "caption"), default=("File 1", "File 2"),
                        help=("Sets the boxplot captions under the single "
                              "boxplots, first for the left, then for the "
                              "right boxplot, divide the captions with a "
                              "space, if you want to incorporate spaces "
                              "put the captions in '', if left out file 1 "
                              "and file 2 will be chosen as captions"))
    args = parser.parse_args()
    create_comparison_COM_plot(output_filename=args.output_filename,
                               output_directory=args.output_directory,
                               input_filepath_whole_1=args.input[0],
                               input_filepath_whole_2=args.input[1],
                               datanames=[args.caption[0], args.caption[1]])


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    make_comparison_boxplots()
