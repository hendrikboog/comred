# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 13:57:50 2018

@author: Hendrik
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import scipy as scp
import scipy.interpolate as interpolate
from matplotlib import gridspec
import string
import argparse
from pathlib import Path
from datetime import date
import progressbar


# ----------------------------------------------------------------------------

def model_real(data, n_bins, n_samples):
    # inverse transition sampling for any own data provided, would need testing
    # for hyperparameter optimization bins >= 1000, n_samples
    hist, bin_edges = np.histogram(data, bins=n_bins, density="True")
    cumulative_values = np.zeros(bin_edges.shape)
    cumulative_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    inv_cdf = interpolate.interp1d(cumulative_values, bin_edges)
    r = np.random.rand(n_samples)
    return inv_cdf(r)


# ----------------------------------------------------------------------------

def model_real_rs(data, n_bins, n_samples, std_samples):
    # inverse transition sampling for any own data provided, would need
    # testing for hyperparameter optimization bins >= 1000, n_samples
    n = np.random.normal(n_samples, std_samples)
    hist, bin_edges = np.histogram(data, bins=n_bins, density="True")
    cumulative_values = np.zeros(bin_edges.shape)
    cumulative_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    inv_cdf = interpolate.interp1d(cumulative_values, bin_edges)
    r = np.random.rand(n)
    return inv_cdf(r)


# ----------------------------------------------------------------------------

def make_histogram_o(data_list, xlabel, title, color, bins, linestyle, label):
    n, outbins, patches = plt.hist(data_list, density=0, bins=bins,
                                   color=color, histtype="step", lw=4,
                                   linestyle=linestyle, label=label)
    plt.ylim(ymin=0)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Count', fontsize=20)
    # plt.grid(axis = "y", color = "0.9", linestyle="-")
    return outbins, n


# ----------------------------------------------------------------------------

def read_comred(filename):
    # reads a comred file and transfers the data into a list
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


# ----------------------------------------------------------------------------

# Reads all files from a directory with the .csv expansion,
# returns list of files
def COM_read_directory_files(directory_name, ending):
    filelist = []
    for filename in os.listdir(directory_name):
        if filename.endswith(ending):
            filelist.append(os.path.join(directory_name, filename))
    return filelist


# ---------------------------------------------------------------------------

def create_data(directory_name):
    filelist = COM_read_directory_files(directory_name, ".comred")
    dataallfiles = []
    datapoint_stat = []
    for element in filelist:
        data_in_file = read_comred(element)
        datapoint_stat.append(len(data_in_file))
        dataallfiles.append(data_in_file)
    x = []
    y = []
    z = []
    vol = []
    meanint = []
    mean_datpoint = np.mean(datapoint_stat)
    std_datpoint = np.std(datapoint_stat)
    for element in dataallfiles:
        for lement in element:
            x.append(float(lement[5]))
            y.append(float(lement[6]))
            z.append(float(lement[7]))
            vol.append(float(lement[8]))
            meanint.append(float(lement[9]))
    return x, y, z, vol, meanint, mean_datpoint, std_datpoint


# ----------------------------------------------------------------------------

def hyperparam_search(data):
    # searches for hyperparameters for available data
    val = {}
    for element in range(1, 10000, 1):
        model = model_real(data, element, len(data))
        histo = histogram_intersection(data, model)
        val[element] = histo
    # bestbins = max(val.items(), key=operator.itemgetter(1))[0]
    bestvalue = max(val.items())
    return bestvalue


# ----------------------------------------------------------------------------

def histogram_intersection(hist1, hist2):  # not working correctly??
    minima = np.minimum(hist1, hist2)
    intersection = np.sum(minima) / np.sum(hist2)
    return intersection


# ----------------------------------------------------------------------------

def hellinger(f1, f2):
    # f1 and f2 need to be np.histogram hists with the same bin edges
    hell = 0
    for element in range(0, len(f1) - 1, 1):
        hell += (np.sqrt(f1[element]) - np.sqrt(f2[element]))**2
    hell = 1 / (np.sqrt(2)*len(f1)) * np.sqrt(hell)
    return hell


# ----------------------------------------------------------------------------

def double_hist(f1, f2, xmin, xmax):
    y1, bn1 = np.histogram(f1, bins=50)
    y2, bn2 = np.histogram(f2, bins=50)
    y = []
    bins = []
    for el in y1:
        y.append(el)
    for el in y2:
        y.append(el)
    for el in bn1:
        bins.append(el)
    for el in bn2:
        bins.append(el)
    bins.sort()
    binmin = min(bins)
    binmax = max(bins)
    outbins = np.linspace(binmin, binmax)
    o1, n1 = make_histogram_o(f1, "", "", "black", outbins, "-",
                              "Original data")
    o2, n2 = make_histogram_o(f2, "", "", "orange", outbins, "-",
                              "Model data")
    n = []
    for el in n1:
        n.append(el)
    for el in n2:
        n.append(el)
    # if len(f1) == len(f2):
    cosine = scp.spatial.distance.pdist([f1, f2], metric="cosine")[0]
    # else:
    #    cosine = "Not available: Unequal vector length"
    hist1, bins1 = np.histogram(f1, bins=1000)
    hist2, bins2 = np.histogram(f2, bins=bins1)
    histol = histogram_intersection(hist1, hist2)
    hell = hellinger(hist1, hist2)
    textges = ("Cosine distance: " + str(round(cosine, 5)) +
               "\nHistogram overlap: " + str(round(histol, 5)) +
               "\nHellinger distance: " + str(round(hell, 5)))
    # xbound = binmax - 0.14*binmax
    handle1 = Line2D([], [], c='black')
    handle2 = Line2D([], [], c='orange')
    plt.legend(handles=[handle1, handle2], labels=["Experimental data",
                                                   "Simulated data"],
               loc="upper left", edgecolor="black", fontsize=15)
    axes = plt.gca()
    axes.set_ylim([0, 1.1*max(n)])
    axes.set_xlim([xmin, xmax])
    axes.tick_params(labelsize=15)
    # ybound = 1.1*max(n) - 0.05*1.1*max(n)
    plt.text(0.98, 0.98, textges, size=15, transform=axes.transAxes,
             ha="right", va="top", bbox=dict(boxstyle="round", ec="black",
                                             facecolor="white", alpha=0.5))


# ----------------------------------------------------------------------------

def compare_histograms(x, y, z, vol, meanint, output_filepath, title,
                       n, change="shift", num_bins=1000):
    if change == "shift":
        xchanged = change_data_shift(x, n)
        xchanged = model_real(xchanged, num_bins, len(xchanged))
        ychanged = change_data_shift(y, n)
        ychanged = model_real(ychanged, num_bins, len(ychanged))
        zchanged = change_data_shift(z, n)
        zchanged = model_real(zchanged, num_bins, len(zchanged))
        volchanged = change_data_shift(vol, n)
        volchanged = model_real(volchanged, num_bins, len(volchanged))
        meanintchanged = change_data_shift(meanint, n)
        meanintchanged = model_real(meanintchanged, num_bins,
                                    len(meanintchanged))
    elif change == "spread":
        xchanged = change_data_spread(x, n)
        xchanged = model_real(xchanged, num_bins, len(xchanged))
        ychanged = change_data_spread(y, n)
        ychanged = model_real(ychanged, num_bins, len(ychanged))
        zchanged = change_data_spread(z, n)
        zchanged = model_real(zchanged, num_bins, len(zchanged))
        volchanged = change_data_spread(vol, n)
        volchanged = model_real(volchanged, num_bins, len(volchanged))
        meanintchanged = change_data_spread(meanint, n)
        meanintchanged = model_real(meanintchanged, num_bins,
                                    len(meanintchanged))
    else:
        raise Exception("Type '{}' not recognized.".format(change),
                        " Use 'shift' or 'spread'")
    xmin = min(min(xchanged), min(x))
    xmax = max(max(xchanged), max(x))
    ymin = min(min(ychanged), min(y))
    ymax = max(max(ychanged), max(y))
    zmin = min(min(zchanged), min(z))
    zmax = max(max(zchanged), max(z))
    volmin = min(min(volchanged), min(vol))
    volmax = max(max(volchanged), max(vol))
    meanmin = min(min(meanintchanged), min(meanint))
    meanmax = max(max(meanintchanged), max(meanint))
    with PdfPages(output_filepath) as pdf:
        fig = plt.figure(figsize=(10, 40))
        plt.suptitle(title, fontsize=25, wrap=True)
        ax = plt.subplot(5, 1, 1)
        plt.text(0.5, 1.05, "X-Coordinate", fontsize=20,
                 horizontalalignment="center", verticalalignment="center",
                 transform=ax.transAxes)
        double_hist(x, xchanged, xmin, xmax)
        ax2 = plt.subplot(5, 1, 2)
        plt.text(0.5, 1.05, "Y-Coordinate", fontsize=20,
                 horizontalalignment="center", verticalalignment="center",
                 transform=ax2.transAxes)
        double_hist(y, ychanged, ymin, ymax)
        ax3 = plt.subplot(5, 1, 3)
        plt.text(0.5, 1.05, "Z-Coordinate", fontsize=20,
                 horizontalalignment="center", verticalalignment="center",
                 transform=ax3.transAxes)
        double_hist(z, zchanged, zmin, zmax)
        ax4 = plt.subplot(5, 1, 4)
        plt.text(0.5, 1.05, "Volume", fontsize=20,
                 horizontalalignment="center", verticalalignment="center",
                 transform=ax4.transAxes)
        double_hist(vol, volchanged, volmin, volmax)
        ax5 = plt.subplot(5, 1, 5)
        plt.text(0.5, 1.05, "Mean Intensity", fontsize=20,
                 horizontalalignment="center", verticalalignment="center",
                 transform=ax5.transAxes)
        double_hist(meanint, meanintchanged, meanmin, meanmax)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        pdf.savefig()
        plt.close(fig)


# ----------------------------------------------------------------------------

def write_files_receptor(input_filepath, output_filepath, output_filename,
                         num_bins=1000, n=0, change_type="shift", n2=0):
    filelist = COM_read_directory_files(input_filepath, ".comred")
    x, y, z, vol, meanint, mean_datpoint, std_datpoint = create_data(
        input_filepath)
    if change_type == "shift":
        nx = change_data_shift(x, n)
        ny = change_data_shift(y, n)
        nz = change_data_shift(z, n)
        nvol = change_data_shift(vol, n2)
        nmi = change_data_shift(meanint, n2)
    elif change_type == "spread":
        nx = change_data_spread(x, n)
        ny = change_data_spread(y, n)
        nz = change_data_spread(z, n)
        nvol = change_data_shift(vol, n2)
        nmi = change_data_shift(meanint, n2)
    print("Write receptor {} files...".format(change_type))
    # all vectors (x,y,z, etc) should have the same amount of data point
    with progressbar.ProgressBar(max_value=len(filelist)) as bar:
        for element in range(1, len(filelist)+1, 1):
            fullpath = (output_filepath + "\\" + output_filename +
                        str(element) + ".comred")
            with open(fullpath, "w") as file:
                file.write("#Count\tX\tY\tZ\tVolume\tMean Intensity\n")
                number_data_per_file = int(np.random.normal(mean_datpoint,
                                                            std_datpoint))
                modelx = model_real(nx, num_bins, number_data_per_file)
                modely = model_real(ny, num_bins, number_data_per_file)
                modelz = model_real(nz, num_bins, number_data_per_file)
                modelvol = model_real(nvol, num_bins, number_data_per_file)
                modelmeanint = model_real(nmi, num_bins, number_data_per_file)
                for el in range(0, number_data_per_file, 1):
                    file.write(str(el) + "\t" + str(modelx[el]) + "\t" +
                               str(modely[el]) + "\t" + str(modelz[el]) +
                               "\t" + str(modelvol[el]) + "\t" +
                               str(modelmeanint[el]) + "\n")
            bar.update(element)
    print("Done")


# ----------------------------------------------------------------------------

def write_files_nucleus(input_filepath, output_filepath, output_filename,
                        num_bins=1000, n=0):
    filelist = COM_read_directory_files(input_filepath, ".comred")
    x, y, z, vol, meanint, mean_datpoint, std_datpoint = create_data(
        input_filepath)
    nx = change_data_shift(x, n)
    ny = change_data_shift(y, n)
    nz = change_data_shift(z, n)
    nvol = change_data_shift(vol, n)
    nmi = change_data_shift(meanint, n)
    number_data_per_file = int(len(x) / len(filelist))
    print("Write reference files...")
    # all vectors (x,y,z, etc) should have the same amount of data points
    with progressbar.ProgressBar(max_value=len(filelist)) as bar:
        for element in range(1, len(filelist)+1, 1):
            fullpath = (output_filepath + "\\" + output_filename +
                        str(element) + ".comred")
            with open(fullpath, "w") as file:
                file.write("#Count\tX\tY\tZ\tVolume\tMean Intensity\n")
                modelx = model_real(nx, num_bins, number_data_per_file)
                modely = model_real(ny, num_bins, number_data_per_file)
                modelz = model_real(nz, num_bins, number_data_per_file)
                modelvol = model_real(nvol, num_bins, number_data_per_file)
                modelmeanint = model_real(nmi, num_bins, number_data_per_file)
                for el in range(0, number_data_per_file, 1):
                    file.write(str(el) + "\t" + str(modelx[el]) + "\t" +
                               str(modely[el]) + "\t" + str(modelz[el]) +
                               "\t" + str(modelvol[el]) + "\t" +
                               str(modelmeanint[el]) + "\n")
            bar.update(element)
    print("Done")


# ----------------------------------------------------------------------------

def create_output_folders(output_directory, directory_list):
    directories = directory_list
    headlist = (date.today().strftime("%Y%m%d"), "comred_data_simulation")
    headfolder = "_".join(headlist)
    mainoutdir = os.path.join(output_directory, headfolder)
    for item in directories:
        dirname = "n=" + item
        finalpath = os.path.join(mainoutdir, dirname)
        path_receptors = os.path.join(finalpath, "receptors")
        path_reference = os.path.join(finalpath, "reference")
        path_recshift = os.path.join(path_receptors, "shift")
        path_recspread = os.path.join(path_receptors, "spread")
        try:
            Path(finalpath).mkdir(parents=True)
            Path(path_recshift).mkdir(parents=True)
            Path(path_recspread).mkdir(parents=True)
            Path(path_reference).mkdir(parents=True)
        except FileExistsError:
            print("Folder already exists")
            break
    return mainoutdir


# ----------------------------------------------------------------------------

def change_data_shift(data, n):
    new_data = []
    mean = np.mean(data)
    for x in data:
        y = x + n*mean
        new_data.append(y)
    return new_data


# ----------------------------------------------------------------------------

def change_data_spread(data, n):
    # like change_data_gauss w. same distribution (input distribution),
    # parameterized between -x and x
    new_data = []
    distrib = model_real(data, 1000, len(data))
    meandat = np.mean(data)
    hist, edges = np.histogram(data, bins=50)
    ledge = meandat - edges[0]
    redge = edges[-1] - meandat
    if ledge > redge:
        bin_edge = ledge
    else:
        bin_edge = redge
    for counter in range(0, len(data), 1):
        rand = (distrib[counter] - meandat) / bin_edge
        y = data[counter] + n*rand*meandat
        new_data.append(y)
    return new_data


# ----------------------------------------------------------------------------

def receptor_dist_mod_big_plot(input_directory, output_filepath, title,
                               n_list, check=0):
    print("Beginning task...\nReading in data...")
    x, y, z, vol, meanint, mean_datpoint, std_datpoint = create_data(
        input_directory)
    num_bins = 1000
    print("Creating histograms")
    with PdfPages(output_filepath) as pdf:
        sizey = 8 * len(n_list)
        fig = plt.figure(figsize=(25, sizey))
        height_ratios = [2]
        for element in n_list:
            height_ratios.append(4)
        gs = gridspec.GridSpec(len(n_list) + 1, 3, width_ratios=[1, 4, 4],
                               height_ratios=height_ratios)
        plt.suptitle(title, fontsize=35)
        # plt.subplot(gs[0])
        form1 = plt.subplot(gs[1])
        form1txt = r"$ x_{new} = x_{real} + n \cdot \overline{x} $"
        form1.text(0.5, 0.5, form1txt, fontsize=30,
                   horizontalalignment="center", verticalalignment="center")
        form1.axis("off")
        form2 = plt.subplot(gs[2])
        form2txt = (r"$ x_{new} = x_{real} + $"
                    r"$\frac{n \cdot \overline{x} \cdot (x_{rnd} - "
                    r"\overline{x})}{b} $")
        form2.text(0.5, 0.5, form2txt, fontsize=30,
                   horizontalalignment="center", verticalalignment="center")
        form2.axis("off")
        counter = 3
        models_shift = []
        models_spread = []
        maxs = []
        mins = []
        for element in n_list:
            element = float(element)
            nc1 = change_data_shift(x, element)
            nc2 = change_data_spread(x, element)
            nm1 = model_real(nc1, num_bins, len(nc1))
            nm2 = model_real(nc2, num_bins, len(nc2))
            models_shift.append(nm1)
            models_spread.append(nm2)
            if len(nm1) and len(nm2):
                maxs.append(max(max(nm1), max(nm2)))
                mins.append(min(min(nm1), min(nm2)))
            else:
                raise Exception("Data simulation contains an empty list."
                                " Please check the input folders - "
                                "Do they contain .comred files?")
            if check:
                parent_path = os.path.dirname(output_filepath)
                check_filename_spread = "check_n={}_spread.pdf".format(element)
                check_filename_shift = "check_n={}_shift.pdf".format(element)
                check_output_shift = os.path.join(parent_path,
                                                  check_filename_shift)
                check_output_spread = os.path.join(parent_path,
                                                   check_filename_spread)
                title_spreadlist = [title, "Inverse Transform Sampling check",
                                    "spread n={}".format(element)]
                title_shiftlist = [title, "Inverse Transform Sampling check",
                                   "spread n={}".format(element)]
                title_check_spread = " ".join(title_spreadlist)
                title_check_shift = " ".join(title_shiftlist)
                compare_histograms(x, y, z, vol, meanint,
                                   output_filepath=check_output_shift,
                                   title=title_check_shift, n=element,
                                   change="shift")
                compare_histograms(x, y, z, vol, meanint,
                                   output_filepath=check_output_spread,
                                   title=title_check_spread, n=element,
                                   change="spread")
        xmax = max(maxs)
        xmin = min(mins)
        for shift, spread, element in zip(models_shift, models_spread, n_list):
            lefttext = plt.subplot(gs[counter])
            lefttxt = (string.ascii_lowercase[int(counter/3)-1] + ") " +
                       "n = " + str(element))
            lefttext.text(0.5, 0.5, lefttxt, fontsize=35,
                          horizontalalignment="center",
                          verticalalignment="center")
            lefttext.axis("off")
            counter += 1
            plt.subplot(gs[counter])
            double_hist(x, shift, xmin, xmax)
            counter += 1
            plt.subplot(gs[counter])
            double_hist(x, spread, xmin, xmax)
            counter += 1
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        pdf.savefig()
        plt.close(fig)
    print("Task completed")


# ----------------------------------------------------------------------------

def write_files(output_directory, input_directory, reference_directory,
                filehandle, n_list):
    print("Initializing writing of data...\nCreating output folders...")
    mainoutdir = create_output_folders(output_directory, n_list)
    for n in n_list:
        nf = float(n)
        dirname = "n=" + n
        recshiftdir = os.path.join(mainoutdir, dirname, "receptors", "shift")
        print("Calculating shifted receptor distributions with"
              " n={}...".format(n))
        write_files_receptor(input_directory, recshiftdir, filehandle, n=nf,
                             change_type="shift")
        recspreaddir = os.path.join(mainoutdir, dirname, "receptors", "spread")
        print("Calculating spread out receptor",
              " distributions with n={}...".format(n))
        write_files_receptor(input_directory, recspreaddir, filehandle, n=nf,
                             change_type="spread")
        referencedir = os.path.join(mainoutdir, dirname, "reference")
        print("Calculating reference points...")
        write_files_nucleus(reference_directory, referencedir, filehandle, n=0)


# ----------------------------------------------------------------------------

def data_simulation():
    description = ("ComRed 0.3 - Last update: 2021-05-14 - ComReds Data"
                   " simulation script - The script uses original data to "
                   "generate simulated distributions via inverse transform"
                   " sampling.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-o", "--output", metavar=("Output directory"),
                        help=("Set the output directory, where results"
                              " will be saved, required parameter"),
                        required=True)
    parser.add_argument("-f", "--filename", metavar=("Output pdf filename"),
                        help=("Set an output filename for the checkup "
                              "plots pdf"),
                        required=True)
    parser.add_argument('-i', "--input", metavar="Input directory",
                        help=("Input directory filepath. All .comred files"
                              " in the input directory will be read."),
                        required=True,  type=str)
    parser.add_argument('-t', "--title", metavar="Title",
                        help=("Title of the overview plot. Wrap the title "
                              "in '' to allow spaces in the title."), type=str,
                        required=True)
    parser.add_argument("-n", "--n_list", metavar=("List of n-values"),
                        help=("List of n-values for doing simulations."),
                        required=True, type=str, nargs="+")
    parser.add_argument("-w", "--write", action="store_true",
                        help=("Toggles writing files to single folder before"
                              " when using data simulation plots."))
    parser.add_argument("-a", "--filehandle", metavar=("Filehandle for"
                                                       " writing files"),
                        help=("Toggles writing files to single folder before"
                              " when using data simulation plots. Only"
                              " required, when using -w to write files."),
                        type=str)
    parser.add_argument('-r', "--reference", metavar="Reference directory",
                        help=("Directory filepath for your reference files."
                              " All .comred files in the input directory will"
                              " be read. Only required with -w to write"
                              " files"),
                        required=True,  type=str)
    parser.add_argument("-c", "--control", action="store_true",
                        help=("Toggles writing a control plot to check "
                              "the accuracy of the inverse transform "
                              "sampling."))
    args = parser.parse_args()
    pdf_handle = args.filename + ".pdf"
    pdf_filepath = os.path.join(args.output, pdf_handle)
    receptor_dist_mod_big_plot(input_directory=args.input,
                               output_filepath=pdf_filepath,
                               title=args.title, n_list=args.n_list,
                               check=args.control)
    if args.write:
        write_files(output_directory=args.output, input_directory=args.input,
                    reference_directory=args.reference,
                    filehandle=args.filehandle, n_list=args.n_list)


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    data_simulation()
