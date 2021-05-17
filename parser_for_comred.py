#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import sys
import progressbar


# -----------------------------------------------------------------------------

# Module that parses one csv or tsv file
def parse_sv(input_filename, column_x, column_y, column_z, column_vol,
             column_meanint, counting=0, output_filename=False,
             dattype=".csv", verbose=0):
    with open(input_filename, "r") as csv_file:
        iter_csv_file = iter(csv_file)
        next(iter_csv_file)
        if counting:
            column_x -= counting
            column_y -= counting
            column_z -= counting
            column_meanint -= counting
            column_vol -= counting
        if output_filename:
            output_filename = output_filename
            with open(output_filename, "w") as output:
                output.write(
                        "#Count\tX\tY\tZ\tVolume\tMean Intensity\n")
                for count, line in enumerate(iter_csv_file):
                    if dattype in ("csv", ".csv"):
                        line = line.split(",")
                    elif dattype in ("tsv", "tsv"):
                        line = line.split("\t")
                    else:
                        print("""Format not recognized try csv or tsv
                              files""")
                        sys.exit()
                    values = [count, line[column_x], line[column_y],
                              line[column_z], line[column_vol],
                              line[column_meanint]]
                    values = ("\t".join(str(val) for val in values) + "\n")
                    output.write(values)
        else:
            for count, line in enumerate(iter_csv_file):
                output = []
                if dattype in ("csv", ".csv"):
                    line = line.split(",")
                elif dattype in ("tsv", "tsv"):
                    line = line.split("\t")
                values = [count, line[column_x], line[column_y],
                          line[column_z], line[column_vol],
                          line[column_meanint]]
                output.append(values)
                if verbose:
                    print(("File {}"
                           " successfully parsed").format(input_filename))
            return output
    if verbose:
        print("File {} successfully parsed".format(input_filename))

# -----------------------------------------------------------------------------

# Module which takes parse_sv and uses it to a whole directory
def parse_sv_directory(directory_name, column_x, column_y, column_z,
                       column_vol, column_meanint, counting=0,
                       output_directory_name=False,
                       output_filename=False, dattype=".csv", verbose=0):
    filelist = []
    column_x -= counting
    column_y -= counting
    column_z -= counting
    column_meanint -= counting
    column_vol -= counting
    for filename in os.listdir(directory_name):
        if filename.endswith(dattype):
            filelist.append(os.path.join(directory_name, filename))
    if not output_directory_name:
        datalist = []
    with progressbar.ProgressBar(max_value=len(filelist)) as bar:
        for counter, file in enumerate(filelist):
            if output_directory_name:
                if output_filename:
                    output_filename = os.path.join(output_directory_name,
                                                   output_filename + "_"
                                                   + str(counter)
                                                   + ".comred")
                    parse_sv(file, column_x, column_y, column_z, column_vol,
                             column_meanint, output_filename=output_filename,
                             dattype=dattype, verbose=verbose)
                    output_filename = False
                    bar.update(counter)
                else:
                    output_filename = os.path.join(
                        output_directory_name,
                        os.path.split(file)[1] + "_output" + ".comred")
                    parse_sv(file, column_x, column_y, column_z, column_vol,
                             column_meanint, output_filename=output_filename,
                             dattype=dattype, verbose=verbose)
                    output_filename = False
                    bar.update(counter)
            else:
                datalist.append(parse_sv(file, column_x, column_y, column_z,
                                         column_vol, column_meanint, dattype,
                                         verbose=verbose))
    if not output_directory_name:
        if verbose:
            print(("Directory {}"
                  " successfully parsed").format(directory_name))
        return datalist
    if verbose:
        print("Directory {} successfully parsed".format(directory_name))


# -----------------------------------------------------------------------------

# Main module - takes in sys arguments and starts up the other functions
def main():
    description = ("ComRed 0.3 - Last update: 2021-05-14 - "
                   "Parser for tsv or csv files")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-d", "--directory", metavar="Input directory",
                        help=("Only use either -d or -f, files are read in"
                              " automatically"))
    parser.add_argument("-f", "--file", metavar="Input file",
                        help="Only use either -d or -f")
    parser.add_argument("-t", "--type", choices=[".csv", "csv", ".tsv", "tsv"],
                        default=".csv", metavar=("filetype"))
    parser.add_argument("-x", "--x_column", type=int, required=True,
                        help=("Column-number for x-coordinate data of"
                              " receptors, Required argument"),
                        metavar="X Column")
    parser.add_argument("-y", "--y_column", type=int, required=True,
                        help=("Column-number for y-coordinate data of"
                              " receptors, Required argument"),
                        metavar="Y Column")
    parser.add_argument("-z", "--z_column", type=int, required=True,
                        help=("Column-number for z-coordinate data of"
                              " receptors, Required argument"),
                        metavar="Z Column")
    parser.add_argument("-s", "--size_column", type=int, required=True,
                        help=("Column-number for size of receptors"
                              " accumulations in voxels, Required argument"),
                        metavar="Volume column")
    parser.add_argument("-m", "--meanint_column", type=int, required=True,
                        help=("Column-number for mean intensity data of"
                              " receptors, Required argument"),
                        metavar="Mean intensity Column")
    parser.add_argument("-o", "--output", required=True,
                        metavar="Output folder",
                        help=("Output folder for resulting .comred files"))
    parser.add_argument("-c", "--counting", action=("store_true"),
                        help=("Start counting the columns from 1 not 0"))
    parser.add_argument("-v", "--verbose", action=("store_true"),
                        help=("More details while running the program"))
    args = parser.parse_args()
    print("Parsing files...")
    if args.directory:
        parse_sv_directory(
            directory_name=args.directory, column_x=int(args.x_column),
            column_y=int(args.y_column), column_z=int(args.z_column),
            column_vol=int(args.size_column),
            column_meanint=int(args.meanint_column),
            output_directory_name=args.output, counting=args.counting,
            dattype=args.type, verbose=args.verbose)
    elif args.file:
        parse_sv(args.file, args.x_column, args.y_column, args.z_column,
                 args.size_column, args.meanint_column,
                 output_file_name=args.output, counting=args.counting,
                 dattype=args.type, verbose=args.verbose)
    else:
        print("Please provide a file (-f) or a directory (-d)\n")
        sys.exit()
    print("Finished parsing files")


# -----------------------------------------------------------------------------

# Required for starting with the script with the command line
if __name__ == "__main__":
    main()
