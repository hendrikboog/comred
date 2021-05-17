# ComRed - Center Of Mass REceptor Distribution 

COMRED is a receptor distribution measurement tool, that uses center of mass calculations. It is intended to be used to get differences in distribution between pools of microscopy images of fixed cells. It is usually used with .csv files produced by ilastik, but might also be used with files from other programs.

## Installation

### Requirements
* python 3.8 or higher
* numpy
* scipy
* matplotlib
* argparse
* progressbar2

### Procedure
Install using `git clone https://github.com/hendrikboog/comred.git` or download and unpack this project.

## Usage
You need two different files for your positions for any comparisons - for every cell you need a file for the receptor positions and a file for your reference position. If you have several types of receptors you need to provide them in a new file each.
First parse your files from ilastik (or another program) using the parse command for ComRed. You need to have a .csv or .tsv file, that provides a x, y, and z coordinates a volume for every object (in pixel) and a mean intensity value for pixels in the object.
	
Get help on any command of ComRed: `python comred.py [command] -h`
For example: `python comred.py parse -h`
```
usage: comred.py parse [-h] [-d Input directory] [-f Input file] [-t filetype] -x X Column -y Y Column -z Z Column -s Volume column -m Mean intensity Column -o Output folder [-c] [-v]

ComRed 0.3 - Parser for tsv or csv files

optional arguments:
  -h, --help            show this help message and exit
  -d Input directory, --directory Input directory
                        Only use either -d or -f, files are read in automatically
  -f Input file, --file Input file
                        Only use either -d or -f
  -t filetype, --type filetype
  -x X Column, --x_column X Column
                        Column-number for x-coordinate data of receptors, Required argument
  -y Y Column, --y_column Y Column
                        Column-number for y-coordinate data of receptors, Required argument
  -z Z Column, --z_column Z Column
                        Column-number for z-coordinate data of receptors, Required argument
  -s Volume column, --size_column Volume column
                        Column-number for size of receptors accumulations in voxels, Required argument
  -m Mean intensity Column, --meanint_column Mean intensity Column
                        Column-number for mean intensity data of receptors, Required argument
  -o Output folder, --output Output folder
                        Output folder for resulting .comred files
  -c, --counting        Start counting the columns from 1 not 0
  -v, --verbose         More details while running the program
```

After you've parsed you should get two folders with .comred files, one for the receptors, one for the reference. You can analyze these with the calc command.

```
usage: comred.py calc [-h] -o Output directory [-f Output filename] -r Receptor directory -n Nucleus/Reference directory -x X-resolution -y Y-resolution -z Z-resolution
                      [-c Plot color] [-v]

ComRed 0.3 - Script for calculating the Center of Mass distances between COM of the receptors and a reference, as well as the mean distance of
receptors to their Center of Mass

optional arguments:
  -h, --help            show this help message and exit
  -o Output directory, --output_directory Output directory
                        Set the output directory, where results will be saved, required parameter
  -f Output filename, --output_filename Output filename
                        Set a filename scheme for output files
  -r Receptor directory, --receptors Receptor directory
                        Directory containing the receptor files, all .comred files in the directory will be read in, so use a different directory for the reference files,
                        required parameter
  -n Nucleus/Reference directory, --reference Nucleus/Reference directory
                        Directory containing the reference (usually nucleus) files, all .comred files in the directory will be read in, so use a different directory for
                        the receptor files, required parameter
  -x X-resolution, --x_resolution X-resolution
                        x Resolution in micrometer per pixel, required parameter
  -y Y-resolution, --y_resolution Y-resolution
                        y Resolution in micrometer per pixel, required parameter
  -z Z-resolution, --z_resolution Z-resolution
                        Z-Resolution in micrometer per pixel, required parameter
  -c Plot color, --color Plot color
                        Color of histograms and boxplots
  -v, --verbose         More details while running the program
```

You will get an output file for the receptor distributions. If you want to compare two different receptor distributions (and such their files you got from calc) you can use the comp command.

```
usage: comred.py comp [-h] -o Output directory [-f Output filename] -i Input file 1 Input file 2 [-c Left boxplot caption Right boxplot caption]

ComRed 0.3 - Script for comparing two sets of receptor datasets and produce boxplot figures

optional arguments:
  -h, --help            show this help message and exit
  -o Output directory, --output_directory Output directory
                        Set the output directory, where results will be saved, required parameter
  -f Output filename, --output_filename Output filename
                        Set a filename scheme for output files
  -i Input file 1 Input file 2, --input Input file 1 Input file 2
                        The filepath to the output files from dist_calc for the comparison, required parameter, divide the two file names with a space and put '' around
                        files with space in the names
  -c Left boxplot caption Right boxplot caption, --caption Left boxplot caption Right boxplot caption
                        Sets the boxplot captions under the single boxplots, first for the left, then for the right boxplot, divide the captions with a space, if you want
                        to incorporate spaces put the captions in '', if left out file 1 and file 2 will be chosen as captions
```

Another way to get a plot is to use the plot command - it uses pair comparisons to generate a lot of plots. However, you will not get the additional information you would get from the comp command. The program is set to generate n subplots with m pairs in it. Usually one subplot would compare two conditions with each other, for example two distributions of receptors in the same cell type would be one pair - but you add another pair for the same comparison in another cell type. Then in the next subplot you would compare the two cell types, but only the distribution of one type of receptor per pair.

```
usage: comred.py plot [-h] -o Output directory -f Output filename -i Input filepath [Input filepath ...] -s Number of subplots -p Number of pairs -t Figure title -c Subplot caption [Subplot caption ...] -a Subtitle [Subtitle ...] -l Condition1 Condition2

ComRed 0.3 - Script for comparing n pairs of receptor distributions. The plot can be divided into subplots (for example for comparing different
distributions) with n pairs that are compared.

optional arguments:
  -h, --help            show this help message and exit
  -o Output directory, --output_directory Output directory
                        Set the output directory, where results will be saved, required parameter
  -f Output filename, --output_filename Output filename
                        Set a filename scheme for output files
  -i Input filepath [Input filepath ...], --input Input filepath [Input filepath ...]
                        Add the input filepath for every input file. Divide different files with a space. Put a a file in '' to allow spaces in the path.
  -s Number of subplots, --subplots Number of subplots
                        Number of subplots in the final plot
  -p Number of pairs, --pairs Number of pairs
                        Number of compared pairs per subplot
  -t Figure title, --title Figure title
                        Head title of the final figure
  -c Subplot caption [Subplot caption ...], --caption Subplot caption [Subplot caption ...]
                        Caption for every subplot (for every row), number of captions is equal to the number of rows
  -a Subtitle [Subtitle ...], --subtitle Subtitle [Subtitle ...]
                        Number of subtitles needs to be equal to the number of subplots
  -l Condition1 Condition2, --legends Condition1 Condition2
                        Legends for comparison pairs for every subplot, e.g. BMPRIb BMPRII. Number of legends needs to be equal to twice the number of subplots. For every
                        subplot you need to add more legends by invoking a new -l with two conditions
```

If you want to generate your own simulated data, which will generate shifted or spread experimental data like distributions you can use the sim command. This might be used for example to test which differences in your data ComRed can detect.

```
usage: comred.py sim [-h] -o Output directory -f Output pdf filename -i Input directory -t Title -n List of n-values [List of n-values ...] [-w]
                     [-a Filehandle for writing files] -r Reference directory [-c]

ComRed 0.3 - ComReds Data simulation script - The script uses original data to generate simulated distributions via inverse transform sampling.

optional arguments:
  -h, --help            show this help message and exit
  -o Output directory, --output Output directory
                        Set the output directory, where results will be saved, required parameter
  -f Output pdf filename, --filename Output pdf filename
                        Set an output filename for the checkup plots pdf
  -i Input directory, --input Input directory
                        Input directory filepath. All .comred files in the input directory will be read.
  -t Title, --title Title
                        Title of the overview plot. Wrap the title in '' to allow spaces in the title.
  -n List of n-values [List of n-values ...], --n_list List of n-values [List of n-values ...]
                        List of n-values for doing simulations.
  -w, --write           Toggles writing files to single folder before when using data simulation plots.
  -a Filehandle for writing files, --filehandle Filehandle for writing files
                        Toggles writing files to single folder before when using data simulation plots. Only required, when using -w to write files.
  -r Reference directory, --reference Reference directory
                        Directory filepath for your reference files. All .comred files in the input directory will be read. Only required with -w to write files
  -c, --control         Toggles writing a control plot to check the accuracy of the inverse transform sampling.
```

You can visualize your raw data. After parsing your .csv or .tsv files you can use the visualize command to generate plots for all your parsed files, which shows dots of the found receptor accumulations, where the positions correspond to the positions in their x-y, x-z, or y-z projections and their sizes corresponding to their virtual masses (the weight they have for the final center of mass).

```
usage: comred.py visualize [-h] -o Output directory [-f Output filename] -r Receptor directory -n Nucleus/Reference
                           directory -x X-resolution -y Y-resolution -z Z-resolution

optional arguments:
  -h, --help            show this help message and exit
  -o Output directory, --output_directory Output directory
                        Set the output directory, where results will be saved, required parameter
  -f Output filename, --output_filename Output filename
                        Set a filename scheme for output files
  -r Receptor directory, --receptors Receptor directory
                        Directory containing the receptor files, all .comred files in the directory will be read in,
                        so use a different directory for the reference files, required parameter
  -n Nucleus/Reference directory, --reference Nucleus/Reference directory
                        Directory containing the reference (usually nucleus) files, all .comred files in the directory
                        will be read in, so use a different directory for the receptor files, required parameter
  -x X-resolution, --x_resolution X-resolution
                        x Resolution in micrometer per pixel, required parameter
  -y Y-resolution, --y_resolution Y-resolution
                        y Resolution in micrometer per pixel, required parameter
  -z Z-resolution, --z_resolution Z-resolution
                        Z-Resolution in micrometer per pixel, required parameter
```