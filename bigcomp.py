# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scipy.stats as sps
import math
import argparse
from dist_calc import dist_calc as dc
import progressbar


# -----------------------------------------------------------------------------

def stars(p):
    # used to produce the stars on the boxplots
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
    # calculating the mann-whitney-u-test and effect sizes
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

def invisplot():
    frame1 = plt.gca()
    frame1.axes().get_xaxis().set_visible(False)
    frame1.axes().get_yaxis().set_visible(False)


# -----------------------------------------------------------------------------

def boxplot_all_bw(ylabel, title, datanames, datalists, legend):
    # datalists shall have the form [[data1, data2], [data1, data2], ...],
    # that you want to have compared
    x = 1
    tickpos = [0]
    yy = []
    ymin = []
    for element in datalists:
        y_max = np.max(np.concatenate((element[0], element[1])))
        y_min = np.min(np.concatenate((element[0], element[1])))
        yy.append(y_max)
        ymin.append(y_min)
    y2 = np.max(yy)
    ymax = y2 + (math.fabs(np.min(ymin)) + math.fabs(y2)) * 0.33
    # * 0.25 or 0.35
    for element in datalists:
        p = mwu_effectsize(element[0], element[1])[0]
        sp = stars(p)
        box = plt.boxplot([element[0], element[1]], positions=[x, (x+1)],
                          patch_artist=True, flierprops={"markersize": 4,
                          "markeredgecolor": "None", "marker": "."},
                          widths=0.9)
        # coloris1 = [colors[0], colors[0], colors[1], colors[1]]
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
        y_max = np.max(np.concatenate((element[0], element[1])))
        y_min = np.min(np.concatenate((element[0], element[1])))
        plt.annotate("", xy=(x, y_max + ymax * 0.01), xycoords='data',
                     xytext=((x+1), y_max + ymax * 0.01), textcoords='data',
                     arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                     connectionstyle="bar,fraction=0.3"))
        plt.text(((2*x+1)/2), y_max + ymax * 0.08, sp, fontsize=14,
                 horizontalalignment='center', verticalalignment='center')
        tickpos.append(((2*x+1)/2))
        x += 3
    axes = plt.gca()
    axes.set_ylim([y_min, ymax])
    axes.set_xlim([0, x-1])
    plt.title(title, fontsize=20)
    plt.ylabel(ylabel, fontsize=14)
    # plt.xlabel("Conditions")
    labels = [""]
    for elementastimo in datanames:
        labels.append(elementastimo)
    plt.xticks(tickpos, labels, rotation=45, fontsize=14)
    plt.yticks(fontsize=12)
    plt.text(0.98, 0.98, legend[0], size=14, transform=axes.transAxes,
             ha="right", va="top", bbox=dict(boxstyle="round", ec="black",
                                             facecolor="white"))
    plt.text(0.98, 0.90, legend[1], size=14, transform=axes.transAxes,
             ha="right", va="top", color="white",
             bbox=dict(boxstyle="round", ec="black", facecolor="black"))


# -----------------------------------------------------------------------------

def create_big_plot_bw(output_filename, output_directory, filelist, rows, cols,
                       figtitle, titles, legends):
    # filelist is a list containing: [[bp_caption, filepath1, filepath2], ...]
    print("Creating plots...")
    output_filepath = output_directory + "/" + output_filename + ".pdf"
    with PdfPages(output_filepath) as pdf:
        with progressbar.ProgressBar(max_value=2*rows) as bar:
            # Creates shift plots
            fig = plt.figure(figsize=(2 * cols + 2.5, 8 * rows))
            plt.suptitle(figtitle, fontsize=18)
            x = 0
            counter = 0
            while x < rows:
                plt.subplot(rows, 1, (x + 1))
                labels = []
                flist = []
                y = 0
                while y < cols:
                    labels.append(filelist[counter][0])
                    flist.append([dc.sort_whole_cell_COM_values(
                            filelist[counter][1])[3],
                            dc.sort_whole_cell_COM_values(
                            filelist[counter][2])[3]])
                    y += 1
                    counter += 1
                subtitles = titles[x]
                boxplot_all_bw(r"$ \Delta $-COM (receptor - nucleus)",
                               subtitles, labels, flist, legends[x])
                bar.update(x)
                x += 1
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig()
            plt.close(fig)
            # Creates spread plots
            fig = plt.figure(figsize=(2 * cols + 2.5, 8 * rows))
            plt.suptitle(figtitle, fontsize=18)
            x = 0
            counter = 0
            while x < rows:
                plt.subplot(rows, 1, (x + 1))
                labels = []
                flist = []
                y = 0
                while y < cols:
                    labels.append(filelist[counter][0])
                    flist.append([dc.sort_distance_COM_values(
                            filelist[counter][1]),
                            dc.sort_distance_COM_values(
                            filelist[counter][2])])
                    y += 1
                    counter += 1
                subtitles = titles[x]
                boxplot_all_bw(r"Mean distance of objects to their COM",
                               subtitles, labels, flist, legends[x])
                bar.update(rows+x)
                x += 1
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig()
            plt.close(fig)
    print("Finished creating plots")


# -----------------------------------------------------------------------------

def make_big_comparison_plot():
    description = ("ComRed 0.3 - Last update: 2021-03-25 - Script for "
                   "comparing n pairs of receptor distributions. The plot "
                   "can be divided into subplots (for example for comparing "
                   "different distributions) with n pairs that are compared.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-o", "--output_directory", metavar=("Output "
                                                             "directory"),
                        help=("Set the output directory, where results"
                              " will be saved, required parameter"),
                        required=True)
    parser.add_argument("-f", "--output_filename", metavar=("Output filename"),
                        help=("Set a filename scheme for output files"),
                        required=True)
    parser.add_argument('-i', "--input", nargs="+", metavar="Input filepath",
                        help=("Add the input filepath for every input file."
                              " Divide different files with a space. Put a"
                              " a file in '' to allow spaces in the path."),
                        required=True,  type=str)
    parser.add_argument("-s", "--subplots", metavar=("Number of subplots"),
                        help=("Number of subplots in the final plot"),
                        required=True, type=int)
    parser.add_argument("-p", "--pairs", metavar=("Number of pairs"),
                        help=("Number of compared pairs per subplot"),
                        required=True, type=int)
    parser.add_argument("-t", "--title", metavar=("Figure title"),
                        help=("Head title of the final figure"),
                        required=True)
    parser.add_argument("-c", "--caption", metavar=("Subplot caption"),
                        nargs="+", required=True,
                        help=("Caption for every subplot (for every row),"
                              " number of captions is equal to the number of"
                              " rows"))
    parser.add_argument("-a", "--subtitle", metavar=("Subtitle"),
                        nargs="+", required=True,
                        help=("Number of subtitles needs to be equal to the"
                              " number of subplots"))
    parser.add_argument("-l", "--legends", nargs=2, action="append",
                        metavar=("Condition1", "Condition2"), required=True,
                        help=("Legends for comparison pairs for every subplot,"
                              " e.g. BMPRIb BMPRII. Number of legends needs "
                              "to be equal to twice the number of subplots."
                              " For every subplot you need to add more"
                              " legends by invoking a new -l with two "
                              "conditions"))
    args = parser.parse_args()
    filelist = []
    # sorts the arguments to be usable by create_big_plot_bw
    legends = [*args.legends]
    for subplot in range(args.subplots):
        for pair in range(0, args.pairs*2, 2):
            subpltlist = [args.caption[int(subplot*args.pairs+pair/2)]]
            subpltlist.append(args.input[subplot*args.pairs*2+pair])
            subpltlist.append(args.input[subplot*args.pairs*2+pair+1])
            filelist.append(subpltlist)
    create_big_plot_bw(output_filename=args.output_filename,
                       output_directory=args.output_directory,
                       filelist=filelist, rows=args.subplots, cols=args.pairs,
                       figtitle=args.title, titles=args.subtitle,
                       legends=legends)


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    make_big_comparison_plot()
