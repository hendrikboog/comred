# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 10:29:15 2021

@author: hboog
"""

import dist_calc as dc
import parser_for_comred as par
import bigcomp as bc
import data_simulation as ds
import two_comp as tc
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import visualize as vis

# -----------------------------------------------------------------------------


def comred_args():
    description = ("ComRed 0.3 - Last update: 2021-05-14\n\n"
                   "Welcome to ComRed - a script to analyze receptor "
                   "distributions using center of mass (COM) calculations\n"
                   "Written by: Hendrik Boog\n\n"
                   "See: https://github.com/hendrikboog/comred \n"
                   "for more info!\n\n"
                   "If you want to analyze your data you first need to\n"
                   "parse your files with the ComRed parser (working with\n"
                   ".tsv and .csv files). Then you can go on and analyze\n"
                   "your files using the calc script. Afterwards you can\n"
                   "either do a pairwise comparison with comp or a\n"
                   "multi-comparsion using plot. If you want to see what\n"
                   "differences in your data you can detect with this\n"
                   "ComRed you can play around with the data simulation\n"
                   "tool: sim.\n\n"
                   "Get help on the tools:\n"
                   "python comred.py tool -h\n"
                   "e.g.: python comred.py parse -h")
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    subparsers = parser.add_subparsers(help=("ComRed script type - options "
                                             "are parse, calc, comp, plot, "
                                             "sim, and visualize."),
                                       dest="script")
    pardes = ("ComRed 0.3 - Last update: 2021-05-14 - Parser for tsv "
              "or csv files")
    parseparse = subparsers.add_parser("parse", description=pardes)
    calcdes = ("ComRed 0.3 - Last update: 2021-05-14 - Script for "
               "calculating the Center of Mass distances between COM "
               "of the receptors and a reference, as well as the "
               "mean distance of receptors to their Center of Mass")
    calcparse = subparsers.add_parser("calc", description=calcdes)
    compdes = ("ComRed 0.3 - Last update: 2021-05-14 - Script for "
               "comparing two sets of receptor datasets and produce "
               "boxplot figures")
    compparse = subparsers.add_parser("comp", description=compdes)
    plotdes = ("ComRed 0.3 - Last update: 2021-05-14 - Script for "
               "comparing n pairs of receptor distributions. The plot "
               "can be divided into subplots (for example for comparing "
               "different distributions) with n pairs that are compared.")
    plotparse = subparsers.add_parser("plot", description=plotdes)
    simdes = ("ComRed 0.3 - Last update: 2021-05-14 - ComReds Data"
              " simulation script - The script uses original data to "
              "generate simulated distributions via inverse transform"
              " sampling.")
    simparse = subparsers.add_parser("sim", description=simdes)
    visdes = ("ComRed 0.3 - Last update: 2021-05-14 - Script for "
              "visualizing comreds calculations. Plots 3 plots per "
              "analyzed cell in x-y, x-z, and y-z planes. The size of "
              "the scatter plots corresponds to their relative virtual "
              "mass (the product of size and mean intensity). "
              "As such larger circles in the plot correspond to a higher"
              " weight in the final COM calculation. The circles color "
              "is randomly pulled from a lookup table and is used to "
              "differ between overlapping circles.")
    visparse = subparsers.add_parser("visualize", description=visdes)
    # parse arguments
    parseparse.add_argument("-d", "--directory", metavar="Input directory",
                            help=("Only use either -d or -f, files are read in"
                                  " automatically"), type=str)
    parseparse.add_argument("-f", "--file", metavar="Input file",
                            help="Only use either -d or -f", type=str)
    parseparse.add_argument("-t", "--type",
                            choices=[".csv", "csv", ".tsv", "tsv"],
                            default=".csv", metavar=("filetype"))
    parseparse.add_argument("-x", "--x_column", type=int, required=True,
                            help=("Column-number for x-coordinate data of"
                                  " receptors, Required argument"),
                            metavar="X Column")
    parseparse.add_argument("-y", "--y_column", type=int, required=True,
                            help=("Column-number for y-coordinate data of"
                                  " receptors, Required argument"),
                            metavar="Y Column")
    parseparse.add_argument("-z", "--z_column", type=int, required=True,
                            help=("Column-number for z-coordinate data of"
                                  " receptors, Required argument"),
                            metavar="Z Column")
    parseparse.add_argument("-s", "--size_column", type=int, required=True,
                            help=("Column-number for size of receptors"
                                  " accumulations in voxels, "
                                  "Required argument"),
                            metavar="Volume column")
    parseparse.add_argument("-m", "--meanint_column", type=int, required=True,
                            help=("Column-number for mean intensity data of"
                                  " receptors, Required argument"),
                            metavar="Mean intensity Column")
    parseparse.add_argument("-o", "--output", required=True,
                            metavar="Output folder",
                            help=("Output folder for resulting .comred files"),
                            type=str)
    parseparse.add_argument("-c", "--counting", action="store_true",
                            help=("Start counting the columns from 1 not 0"))
    parseparse.add_argument("-v", "--verbose", action=("store_true"),
                            help=("More details while running the program"))
    # calc arguments
    calcparse.add_argument("-o", "--output_directory", metavar=("Output "
                                                                "directory"),
                           help=("Set the output directory, where results"
                                 " will be saved, required parameter"),
                           required=True, type=str)
    calcparse.add_argument("-f", "--output_filename", metavar=("Output "
                                                               "filename"),
                           help=("Set a filename scheme for output files"),
                           required=True, type=str)
    calcparse.add_argument("-r", "--receptors", metavar="Receptor directory",
                           help=("Directory containing the receptor files, "
                                 "all .comred files in the directory will be "
                                 "read in, so use a different directory for"
                                 " the reference files, required parameter"),
                           required=True, type=str)
    calcparse.add_argument("-n", "--reference", metavar=("Nucleus/Reference"
                                                         " directory"),
                           help=("Directory containing the reference (usually"
                                 " nucleus) files, all .comred files in the "
                                 "directory will be read in, so use a "
                                 "different directory for the receptor"
                                 " files, required parameter"), required=True,
                           type=str)
    calcparse.add_argument("-x", "--x_resolution", metavar=("X-resolution"),
                           help=("x Resolution in micrometer per pixel"
                                 ", required parameter"),
                           required=True, type=float)
    calcparse.add_argument("-y", "--y_resolution", metavar=("Y-resolution"),
                           help=("y Resolution in micrometer per pixel"
                                 ", required parameter"),
                           required=True, type=float)
    calcparse.add_argument("-z", "--z_resolution", metavar=("Z-resolution"),
                           help=("Z-Resolution in micrometer per pixel"
                                 ", required parameter"),
                           required=True, type=float)
    calcparse.add_argument("-c", "--color", metavar=("Plot color"),
                           help=("Color of histograms and boxplots"),
                           required=True, type=str)
    calcparse.add_argument("-v", "--verbose", action="store_true",
                           help=("More details while running the program"))
    # comp arguments
    compparse.add_argument("-o", "--output_directory", metavar=("Output "
                                                                "directory"),
                           help=("Set the output directory, where results"
                                 " will be saved, required parameter"),
                           required=True)
    compparse.add_argument("-f", "--output_filename",
                           metavar=("Output filename"),
                           help=("Set a filename scheme for output files"))
    compparse.add_argument("-i", "--input", nargs=2,
                           metavar=("Input file 1", "Input file 2"),
                           required=True, default=("File 1", "File 2"),
                           help=("The filepath to the output files from "
                                 "dist_calc for the comparison, required "
                                 "parameter, divide the two file names with "
                                 "a space and put '' around files with "
                                 "space in the names"))
    compparse.add_argument("-c", "--caption", nargs=2,
                           metavar=("Left boxplot caption", "Right boxplot "
                                    "caption"), default=("File 1", "File 2"),
                           help=("Sets the boxplot captions under the single "
                                 "boxplots, first for the left, then for the "
                                 "right boxplot, divide the captions with a "
                                 "space, if you want to incorporate spaces "
                                 "put the captions in '', if left out file 1 "
                                 "and file 2 will be chosen as captions"))
    # plot arguments
    plotparse.add_argument("-o", "--output_directory", metavar=("Output "
                                                                "directory"),
                           help=("Set the output directory, where results"
                                 " will be saved, required parameter"),
                           required=True)
    plotparse.add_argument("-f", "--output_filename",
                           metavar=("Output filename"),
                           help=("Set a filename scheme for output files"),
                           required=True)
    plotparse.add_argument('-i', "--input", nargs="+",
                           metavar="Input filepath",
                           help=("Add the input filepath for every input file."
                                 " Divide different files with a space. Put a"
                                 " a file in '' to allow spaces in the path."),
                           required=True,  type=str)
    plotparse.add_argument("-s", "--subplots", metavar=("Number of subplots"),
                           help=("Number of subplots in the final plot"),
                           required=True, type=int)
    plotparse.add_argument("-p", "--pairs", metavar=("Number of pairs"),
                           help=("Number of compared pairs per subplot"),
                           required=True, type=int)
    plotparse.add_argument("-t", "--title", metavar=("Figure title"),
                           help=("Head title of the final figure"),
                           required=True)
    plotparse.add_argument("-c", "--caption", metavar=("Subplot caption"),
                           nargs="+", required=True,
                           help=("Caption for every subplot (for every row),"
                                 " number of captions is equal to the number "
                                 "of rows"))
    plotparse.add_argument("-a", "--subtitle", metavar=("Subtitle"),
                           nargs="+", required=True,
                           help=("Number of subtitles needs to be equal to the"
                                 " number of subplots"))
    plotparse.add_argument("-l", "--legends", nargs=2, action="append",
                           metavar=("Condition1", "Condition2"), required=True,
                           help=("Legends for comparison pairs for every "
                                 "subplot, e.g. BMPRIb BMPRII. Number of "
                                 "legends needs to be equal to twice the "
                                 "number of subplots. For every subplot you "
                                 "need to add more legends by invoking a "
                                 "new -l with two conditions"))
    # sim arguments
    simparse.add_argument("-o", "--output", metavar=("Output directory"),
                          help=("Set the output directory, where results"
                                " will be saved, required parameter"),
                          required=True)
    simparse.add_argument("-f", "--filename", metavar=("Output pdf filename"),
                          help=("Set an output filename for the checkup "
                                "plots pdf"),
                          required=True)
    simparse.add_argument('-i', "--input", metavar="Input directory",
                          help=("Input directory filepath. All .comred files"
                                " in the input directory will be read."),
                          required=True,  type=str)
    simparse.add_argument('-t', "--title", metavar="Title",
                          help=("Title of the overview plot. Wrap the title "
                                "in '' to allow spaces in the title."),
                          type=str, required=True)
    simparse.add_argument("-n", "--n_list", metavar=("List of n-values"),
                          help=("List of n-values for doing simulations."),
                          required=True, type=str, nargs="+")
    simparse.add_argument("-w", "--write", action="store_true",
                          help=("Toggles writing files to single folder before"
                                " when using data simulation plots."))
    simparse.add_argument("-a", "--filehandle", metavar=("Filehandle for"
                                                         " writing files"),
                          help=("Toggles writing files to single folder before"
                                " when using data simulation plots. Only"
                                " required, when using -w to write files."),
                          type=str)
    simparse.add_argument('-r', "--reference", metavar="Reference directory",
                          help=("Directory filepath for your reference files."
                                " All .comred files in the input directory "
                                "will be read. Only required with -w to write"
                                " files"), type=str)
    simparse.add_argument("-c", "--control", action="store_true",
                          help=("Toggles writing a control plot to check "
                                "the accuracy of the inverse transform "
                                "sampling."))
    # visualize arguments
    visparse.add_argument("-o", "--output_directory", metavar=("Output "
                                                               "directory"),
                          help=("Set the output directory, where results"
                                " will be saved, required parameter"),
                          required=True)
    visparse.add_argument("-f", "--output_filename", metavar=("Output "
                                                              "filename"),
                          help=("Set a filename scheme for output files"))
    visparse.add_argument("-r", "--receptors", metavar="Receptor directory",
                          help=("Directory containing the receptor files, "
                                "all .comred files in the directory will be "
                                "read in, so use a different directory for"
                                " the reference files, required parameter"),
                          required=True)
    visparse.add_argument("-n", "--reference", metavar=("Nucleus/Reference"
                                                        " directory"),
                          help=("Directory containing the reference (usually"
                                " nucleus) files, all .comred files in the "
                                "directory will be read in, so use a "
                                "different directory for the receptor"
                                " files, required parameter"), required=True)
    visparse.add_argument("-x", "--x_resolution", metavar=("X-resolution"),
                          help=("x Resolution in micrometer per pixel"
                                ", required parameter"),
                          required=True, type=float)
    visparse.add_argument("-y", "--y_resolution", metavar=("Y-resolution"),
                          help=("y Resolution in micrometer per pixel"
                                ", required parameter"),
                          required=True, type=float)
    visparse.add_argument("-z", "--z_resolution", metavar=("Z-resolution"),
                          help=("Z-Resolution in micrometer per pixel"
                                ", required parameter"),
                          required=True, type=float)
    args = parser.parse_args()
    if args.script == "parse":
        print("Parsing files...")
        if args.directory:
            par.parse_sv_directory(
                directory_name=args.directory, column_x=int(args.x_column),
                column_y=int(args.y_column), column_z=int(args.z_column),
                column_vol=int(args.size_column),
                column_meanint=int(args.meanint_column),
                output_directory_name=args.output, counting=args.counting,
                dattype=args.type, verbose=args.verbose)
        elif args.file:
            par.parse_sv(args.file, args.x_column, args.y_column,
                         args.z_column, args.size_column, args.meanint_column,
                         output_file_name=args.output, counting=args.counting,
                         dattype=args.type, verbose=args.verbose)
        else:
            print("Please provide a file (-f) or a directory (-d)\n")
            sys.exit()
        print("Finished parsing files")
    elif args.script == "calc":
        dc.make_file_with_COMs(args.output_filename, args.output_directory,
                               args.receptors, args.reference,
                               args.x_resolution, args.y_resolution,
                               args.z_resolution, args.verbose)
        input_filename = f"{args.output_filename}.txt"
        input_filepath_whole_cell = os.path.join(args.output_directory,
                                                 input_filename)
        dc.create_plots_whole_cell_COMs(args.output_filename,
                                        args.output_directory,
                                        input_filepath_whole_cell, args.color)
    elif args.script == "comp":
        tc.create_comparison_COM_plot(output_filename=args.output_filename,
                                      output_directory=args.output_directory,
                                      input_filepath_whole_1=args.input[0],
                                      input_filepath_whole_2=args.input[1],
                                      datanames=[args.caption[0],
                                                 args.caption[1]])
    elif args.script == "plot":
        filelist = []
        legends = [*args.legends]
        for subplot in range(args.subplots):
            for pair in range(0, args.pairs*2, 2):
                subpltlist = [args.caption[int(subplot*args.pairs+pair/2)]]
                subpltlist.append(args.input[subplot*args.pairs*2+pair])
                subpltlist.append(args.input[subplot*args.pairs*2+pair+1])
                filelist.append(subpltlist)
        bc.create_big_plot_bw(output_filename=args.output_filename,
                              output_directory=args.output_directory,
                              filelist=filelist, rows=args.subplots,
                              cols=args.pairs, figtitle=args.title,
                              titles=args.subtitle, legends=legends)
    elif args.script == "sim":
        pdf_handle = f"{args.filename}.pdf"
        pdf_filepath = os.path.join(args.output, pdf_handle)
        ds.receptor_dist_mod_big_plot(input_directory=args.input,
                                      output_filepath=pdf_filepath,
                                      title=args.title, n_list=args.n_list,
                                      check=args.control)
        if args.write:
            ds.write_files(output_directory=args.output,
                           input_directory=args.input,
                           reference_directory=args.reference,
                           filehandle=args.filehandle, n_list=args.n_list)
    elif args.script == "visualize":
        output_filename = f"{args.output_filename}.pdf"
        output_filepath = os.path.join(args.output_directory, output_filename)
        vis.read_out_cells(args.receptors, args.reference, args.x_resolution,
                           args.y_resolution, args.z_resolution,
                           output_filepath)


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    comred_args()
