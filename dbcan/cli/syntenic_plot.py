import argparse
import os

import matplotlib
import matplotlib.pyplot as plt

# import seaborn as sns
from matplotlib import pyplot
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch, Polygon
from matplotlib.path import Path

matplotlib.use("Agg")

plt.style.use("ggplot")


# blocks = ["0-0-30","2-10-50"]


class blastp_hit:
    """
    A class representing a single hit from a BLASTP search.

    This class encapsulates the details of a BLASTP hit, including various metrics like identity percentage, alignment length, mismatches, gaps, start and end positions in both query and subject sequences, e-value, and bitscore.

    Attributes
    ----------
        qseqid (str): Query sequence ID.
        sseqid (str): Subject sequence ID.
        pident (float): Percentage of identical matches.
        length (int): Alignment length.
        mismatch (int): Number of mismatches.
        gapopen (int): Number of gap openings.
        qstart (int): Start of alignment in query.
        qend (int): End of alignment in query.
        sstart (int): Start of alignment in subject.
        send (int): End of alignment in subject.
        evalue (float): E-value.
        bitscore (float): Bit score.
        qlen (int, optional): Length of the query sequence.
        slen (int, optional): Length of the subject sequence.

    Args:
        lines (list): A list containing the BLASTP hit attributes.
    """

    def __init__(self, lines):
        self.qseqid = lines[0]
        self.sseqid = lines[1]
        self.pident = lines[2]
        self.length = int(lines[3])
        self.mismatch = int(lines[4])
        self.gapopen = int(lines[5])
        self.qstart = int(lines[6])
        self.qend = int(lines[7])
        self.sstart = int(lines[8])
        self.send = int(lines[9])
        self.evalue = float(lines[10])
        self.bitscore = float(lines[11])
        if len(lines) >= 13:
            self.qlen = int(lines[12])
        if len(lines) >= 14:
            self.slen = int(lines[13])

    def __repr__(self):
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])

    def __eq__(self, other):
        if self.evalue == self.evalue:
            if self.bitscore == self.bitscore:
                if self.pident == other.pident:
                    return 1

    def __le__(self, other):
        if self.evalue > other.evalue:
            return 1
        elif self.evalue == other.evalue:
            if self.bitscore < other.bitscore:
                return 1
            elif self.bitscore == other.bitscore:
                if self.pident < other.pident:
                    return 1


def identity_map(seqsim):
    """
    Maps a sequence similarity score to a corresponding color.

    This function takes a sequence similarity score (seqsim) and returns a color based on predefined ranges. The color indicates the range in which the similarity score falls.

    Parameters
    ----------
        seqsim (int or float): The sequence similarity score.

    Returns
    -------
        str: A string representing the color corresponding to the sequence similarity score.

    """
    if 80 <= seqsim <= 100:
        color = "red"
    elif 60 <= seqsim <= 80:
        color = "blue"
    elif 40 <= seqsim <= 60:
        color = "green"
    elif 20 <= seqsim <= 40:
        color = "cyan"
    else:
        color = "gray"
    return color


def Get_Position(starts, ends, strands, maxbp, yshift=0, up=1):
    """
    Calculates the positions of genomic elements for visualization.

    This function computes the positions for genomic elements based on their start and end positions, strand orientation, and other parameters. It returns a list of coordinates to be used for plotting these elements.

    Parameters
    ----------
        starts (list of int): List of start positions of genomic elements.
        ends (list of int): List of end positions of genomic elements.
        strands (list of str): List of strand orientations ('+' or '-') for each element.
        maxbp (int): Maximum base pair position for normalization.
        yshift (int, optional): Vertical shift for the plot positions. Defaults to 0.
        up (int, optional): Determines the cluster group for the elements (1 or 2). Defaults to 1.

    Returns
    -------
        tuple: A tuple containing three lists (polygens, blocks, lines) representing the positions for the genomic elements.

    Notes
    -----
        - The function is specifically tailored for genomic data visualization and assumes a certain format and range of input data.
        - This function is part of a larger genomic data visualization tool and may depend on external libraries like matplotlib.
    """
    Width = 1000
    Height = 160
    # px = 1 / plt.rcParams["figure.dpi"]
    poly_heigth = 5
    Triangle_length = 4
    _, plot_start_y = [0, Height / 2 - poly_heigth - yshift]
    polygens = []
    shfit_pos = starts[0]
    for i in range(len(starts)):
        starts[i] = starts[i] - shfit_pos
        ends[i] = ends[i] - shfit_pos
    # maxbp = max(ends) - min(starts)
    pixeachbp = Width / maxbp
    ###  5     4
    ###            3
    ###  1     2

    ###      5     4
    ###  1
    ###      2     3

    blocks = []
    lines = []
    for i in range(len(starts)):
        if strands[i] == "+":
            positions_str = str(starts[i] * pixeachbp) + " " + str(plot_start_y) + " "  ## first point x,y
            positions_str += str(ends[i] * pixeachbp - Triangle_length) + " " + str(plot_start_y) + " "  ## second point
            if up == 1:  ###cluster 1
                blocks.append(positions_str)

            positions_str += str(ends[i] * pixeachbp) + " " + str(plot_start_y + poly_heigth) + " "  ## 3
            positions_str += (
                str(ends[i] * pixeachbp - Triangle_length) + " " + str(plot_start_y + 2 * poly_heigth) + " "
            )  ### 4
            positions_str += str(starts[i] * pixeachbp) + " " + str(plot_start_y + 2 * poly_heigth)

            positions_str1 = str(starts[i] * pixeachbp) + " " + str(plot_start_y + 2 * poly_heigth) + " "
            positions_str1 += (
                str(ends[i] * pixeachbp - Triangle_length) + " " + str(plot_start_y + 2 * poly_heigth) + " "
            )  ### 5

            if up == 2:  ### cluster 2
                blocks.append(positions_str1)

        if strands[i] == "-":
            positions_str = str(starts[i] * pixeachbp) + " " + str(plot_start_y + poly_heigth) + " "

            positions_str += str(starts[i] * pixeachbp + Triangle_length) + " " + str(plot_start_y) + " "
            positions_str += str(ends[i] * pixeachbp) + " " + str(plot_start_y) + " "

            positions_str1 = str(ends[i] * pixeachbp) + " " + str(plot_start_y) + " "
            positions_str1 += str(starts[i] * pixeachbp + Triangle_length) + " " + str(plot_start_y) + " "
            if up == 1:
                blocks.append(positions_str1)

            positions_str += str(ends[i] * pixeachbp) + " " + str(plot_start_y + 2 * poly_heigth) + " "
            positions_str += str(starts[i] * pixeachbp + Triangle_length) + " " + str(plot_start_y + 2 * poly_heigth)

            positions_str1 = str(ends[i] * pixeachbp) + " " + str(plot_start_y + 2 * poly_heigth) + " "
            positions_str1 += str(starts[i] * pixeachbp + Triangle_length) + " " + str(plot_start_y + 2 * poly_heigth)
            if up == 2:
                blocks.append(positions_str1)

        # print (positions_str)
        polygens.append(positions_str)

        ### for genome line
        if i < len(starts) - 1:
            positions_str = str(ends[i] * pixeachbp) + " " + str(plot_start_y + poly_heigth) + " "
            positions_str += str(starts[i + 1] * pixeachbp) + " " + str(plot_start_y + poly_heigth)
            lines.append(positions_str)

    return polygens, blocks, lines


def plot_Polygon_homologous(polygens1, polygens2, types1, types2, ax):
    """
    Plots polygons representing homologous genomic regions on a given axis.

    This function takes two sets of polygon definitions and their corresponding types, and plots these polygons on the provided matplotlib axis. Each type of polygon is colored differently based on a predefined color map.

    Parameters
    ----------
        polygens1 (list of str): List of space-separated strings defining polygons for the first set of genomic regions.
        polygens2 (list of str): List of space-separated strings defining polygons for the second set of genomic regions.
        types1 (list of str): List of types corresponding to each polygon in `polygens1`. Each type is mapped to a specific color.
        types2 (list of str): List of types corresponding to each polygon in `polygens2`. Each type is mapped to a specific color.
        size (int or tuple): Size of the plot or figure.
        ax (matplotlib.axes.Axes): The matplotlib axis on which the polygons will be plotted.

    Side Effects:
        - Modifies the provided matplotlib axis (`ax`) by adding polygons to it.

    Notes
    -----
        - The function is designed specifically for visualizing genomic data and assumes specific formats for the input data.
        - The color map for polygon types is hardcoded within the function.
    """
    colors_map = {
        "CAZyme": "#FF0000",
        "null": "#808080",
        "other": "#808080",
        "TC": "#9400D3",
        "CDS": "#00FFFF",
        "STP": "#0000FF",
        "TF": "#1E90FF",
    }
    for j in range(len(polygens1)):
        polygen = polygens1[j].split()
        points = []
        color = colors_map[types1[j]]
        for i in range(int(len(polygen) / 2)):
            points.append([float(polygen[2 * i]), float(polygen[2 * i + 1])])
        ax.add_patch(Polygon(points, color=color, alpha=0.5))

    for j in range(len(polygens2)):
        polygen = polygens2[j].split()
        points = []
        color = colors_map[types2[j]]
        for i in range(int(len(polygen) / 2)):
            points.append([float(polygen[2 * i]), float(polygen[2 * i + 1])])
        ax.add_patch(Polygon(points, color=color, alpha=0.5))


def decode_block(block):
    """
    Decodes a synteny block string into its constituent elements.

    This function splits a string representation of a synteny block into its order and similarity components.

    Parameters
    ----------
        block (str): A string representing a synteny block, formatted as 'order1-order2-sim'.

    Returns
    -------
        tuple: A tuple containing two integers (order1, order2) and a float (sim), representing the decoded synteny block.
    """
    order1, order2, sim = block.split("-")
    return int(order1), int(order2), float(sim)


def points2(coord):
    """
    Converts a space-separated string of coordinates into a tuple of floats.

    This function is used for parsing coordinate strings into numerical values.

    Parameters
    ----------
        coord (str): A space-separated string of coordinates, formatted as 'x1 y1 x2 y2'.

    Returns
    -------
        tuple: A tuple containing four float values (x1, y1, x2, y2).
    """
    x1, y1, x2, y2 = coord.split()
    return float(x1), float(y1), float(x2), float(y2)


def Shade_curve(x11, x12, y11, y12, x21, x22, y21, y22, ymid, color):
    """
    Creates a shaded curved patch for plotting.

    This function generates a PathPatch object representing a shaded curve between two sets of coordinates, useful for visualizing connections in genomic data.

    Parameters
    ----------
        x11, x12, y11, y12, x21, x22, y21, y22 (float): Coordinates for defining the curved patch.
        xmid, ymid (float): Midpoint coordinates used in the curve construction.
        color (str): Color code for the patch.

    Returns
    -------
        matplotlib.patches.PathPatch: A PathPatch object representing the shaded curve.
    """
    M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
    pathdata = [
        (M, (x11, y11)),
        (C4, (x11, ymid)),
        (C4, (x21, ymid)),
        (C4, (x21, y21)),
        (L, (x22, y22)),
        (C4, (x22, ymid)),
        (C4, (x12, ymid)),
        (C4, (x12, y12)),
        (CP, (x11, y11)),
    ]
    codes, verts = zip(*pathdata)
    path = Path(verts, codes)
    pp = PathPatch(path, color=color, alpha=0.5, lw=0)

    return pp


def Shade_Line(x11, x12, y11, y12, x21, x22, y21, y22):
    """
    Creates a shaded line patch for plotting.

    This function generates a PathPatch object representing a shaded line between two sets of coordinates.

    Parameters
    ----------
        x11, x12, y11, y12, x21, x22, y21, y22 (float): Coordinates for defining the line patch.
        xmid, ymid (float): Midpoint coordinates, not used in this function but kept for consistency with Shade_curve.

    Returns
    -------
        matplotlib.patches.PathPatch: A PathPatch object representing the shaded line.
    """
    M, _, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
    pathdata = [(M, (x11, y11)), (L, (x21, y21)), (L, (x22, y22)), (L, (x12, y12)), (CP, (x11, y11))]
    codes, verts = zip(*pathdata)
    path = Path(verts, codes)
    pp = PathPatch(path)
    return pp


def plot_syntenic_block(blocks, blocks1_coor, blocks2_coor, ax):
    """
    Plots synteny blocks between two genomic regions on a given axis.

    This function visualizes connections between genomic regions (synteny blocks) on a matplotlib axis.

    Parameters
    ----------
        blocks (list of str): List of synteny block strings to be decoded and plotted.
        blocks1_coor (list of str): List of coordinate strings for the first set of genomic regions.
        blocks2_coor (list of str): List of coordinate strings for the second set of genomic regions.
        ax (matplotlib.axes.Axes): The matplotlib axis on which the synteny blocks will be plotted.

    Side Effects:
        - Modifies the provided matplotlib axis (`ax`) by adding patches to it.
    """
    for block in blocks:
        order1, order2, sim = decode_block(block)
        coord1 = blocks1_coor[order1]
        coord2 = blocks2_coor[order2]
        x11, y11, x12, y12 = points2(coord1)
        x21, y21, x22, y22 = points2(coord2)

        color = identity_map(sim)

        xmid = (x11 + x22 + x21 + x22) / 4
        ymid = (y11 + y22) / 2

        pp = Shade_curve(x11, x12, y11, y12, x21, x22, y21, y22, xmid, ymid, color)
        # pp = Shade_Line(x11,x12,y11,y12,x21,x22,y21,y22,xmid,ymid)
        ax.add_patch(pp)


def plot_genome_line(lines_coor1, lines_coor2, ax):
    """
    Plots genomic lines for two sets of genomic regions on a given axis.

    This function is used for visualizing straight lines representing genomic regions on a matplotlib axis.

    Parameters
    ----------
        lines_coor1 (list of str): List of coordinate strings for the first set of genomic lines.
        lines_coor2 (list of str): List of coordinate strings for the second set of genomic lines.
        ax (matplotlib.axes.Axes): The matplotlib axis on which the genomic lines will be plotted.

    Side Effects:
        - Modifies the provided matplotlib axis (`ax`) by adding lines to it.
    """
    for line in lines_coor1:
        x1, y1, x2, y2 = points2(line)
        ax.add_patch(Polygon([(x1, y1), (x2, y2)], color="black", lw=2))
    for line in lines_coor2:
        x1, y1, x2, y2 = points2(line)
        ax.add_patch(Polygon([(x1, y1), (x2, y2)], color="black", lw=2))


### input: gene cluster1: all starts coordinate, end coodinate, strands,Types,
### input: gene cluster2: all starts coordinate, end coodinate, strands,Types,
### input: blocks information,
### the following is an example of inputs

starts = [63045, 65290, 65974, 69423, 71382, 72920, 73962, 76286, 78315, 80633, 82094, 83100, 85536, 86051]
ends = [65270, 65826, 69375, 71255, 72806, 73918, 76238, 78289, 80606, 81952, 83077, 85391, 86048, 86809]
strands = "+,+,+,+,+,+,+,+,+,+,+,+,-,-".split(",")
Types = "CAZyme,TF,TC,TC,null,STP,CAZyme,CAZyme,CAZyme,null,CAZyme,CAZyme,null,TC".split(",")


starts1 = [
    0,
    2278,
    3299,
    4729,
    7100,
    8466,
    9697,
    10746,
    12338,
    15901,
    17511,
    20825,
    22497,
    23601,
    24981,
    27822,
    29391,
    30525,
    31735,
    32771,
    35152,
    36308,
]
ends1 = [
    2271,
    3247,
    4631,
    6934,
    8450,
    9654,
    10720,
    12330,
    15101,
    17143,
    20805,
    22466,
    23592,
    24792,
    27804,
    29367,
    30519,
    31707,
    32641,
    35093,
    36145,
    36917,
]
strands1 = "-,-,-,-,-,-,-,-,-,-,+,+,+,+,+,+,+,+,+,+,-,+".split(",")
Types1 = [
    "CAZyme",
    "CAZyme",
    "other",
    "CAZyme",
    "other",
    "other",
    "other",
    "TC",
    "TC",
    "other",
    "TC",
    "TC",
    "other",
    "other",
    "TC",
    "TC",
    "other",
    "other",
    "other",
    "CAZyme",
    "STP",
    "TF",
]

### let;s try to give the syntenic block information

### format: cluster order-cluster1 order-sequence indenity

blocks = ["0-0-30", f"1-{len(starts1)-1}-50", f"5-{len(starts1)-2}-84", "13-3-65", "2-11-76"]


def syntenic_plot(starts, starts1, ends, ends1, strands, strands1, Types, Types1, blocks, cgcid, pulid):
    """
    Creates a syntenic plot comparing two genomic regions.

    This function visualizes genomic regions and their synteny, drawing comparisons between two sets of genomic data. It plots polygons representing gene locations and shades areas to represent synteny.

    Parameters
    ----------
        starts (list of int): Start positions for genes in the first genomic region.
        starts1 (list of int): Start positions for genes in the second genomic region.
        ends (list of int): End positions for genes in the first genomic region.
        ends1 (list of int): End positions for genes in the second genomic region.
        strands (list of str): Strand orientations for genes in the first genomic region.
        strands1 (list of str): Strand orientations for genes in the second genomic region.
        Types (list of str): Types of genes in the first genomic region.
        Types1 (list of str): Types of genes in the second genomic region.
        blocks (list of str): Synteny blocks information.
        cgcid (str): Identifier for the first genomic region.
        pulid (str): Identifier for the second genomic region.

    Side Effects:
        - Generates and saves a syntenic plot as an SVG file.
    """
    ### for legends
    custom_lines = [
        Line2D([0], [0], color="red", lw=4, alpha=0.5),
        Line2D([0], [0], color="blue", lw=4, alpha=0.5),
        Line2D([0], [0], color="green", lw=4, alpha=0.5),
        Line2D([0], [0], color="cyan", lw=4, alpha=0.5),
        Line2D([0], [0], color="gray", lw=4, alpha=0.5),
    ]

    # print(starts,starts1,ends,ends1,strands,strands1,Types,Types1,blocks,cgcid,pulid)

    labelcolor = ["red", "blue", "green", "cyan", "gray"]
    labels = ["80-100", "60-80", "40-60", "20-40", "0-20"]

    genecustom_lines = [
        Patch(color="#FF0000", alpha=0.5),
        Patch(color="#808080", alpha=0.5),
        Patch(color="#9400D3", alpha=0.5),
        Patch(color="#0000FF", alpha=0.5),
        Patch(color="#1E90FF", alpha=0.5),
    ]

    genelabelcolor = ["#FF0000", "#808080", "#9400D3", "#0000FF", "#1E90FF"]
    geneslabels = ["CAZyme", "Other", "TC", "STP", "TF"]

    ### for legends

    px = 1 / plt.rcParams["figure.dpi"]  ## px
    Width = 1600
    Height = 320 * 2

    fig = plt.figure(figsize=(Width * px, Height * px * 2 / 4))
    ax = fig.add_subplot(111)

    ### decide which
    maxbp = max([max(ends) - min(starts), max(ends1) - min(starts1)])

    polygens, blocks_coor, lines_coor = Get_Position(starts, ends, strands, maxbp, yshift=0, up=1)
    polygens1, blocks1_coor, lines_coor1 = Get_Position(starts1, ends1, strands1, maxbp, yshift=60, up=2)

    ###

    plot_Polygon_homologous(polygens, polygens1, Types, Types1, 2, ax)

    ###
    plot_syntenic_block(blocks, blocks_coor, blocks1_coor, ax)

    plot_genome_line(lines_coor, lines_coor1, ax)

    ### legend1
    legend1 = pyplot.legend(
        custom_lines,
        labels,
        frameon=False,
        labelcolor=labelcolor,
        loc="upper left",
        title="Identity",
        title_fontsize="x-large",
    )
    ax.add_artist(legend1)

    ### legend2
    legend2 = pyplot.legend(
        genecustom_lines,
        geneslabels,
        frameon=False,
        labelcolor=genelabelcolor,
        loc="lower left",
        title="Gene",
        title_fontsize="x-large",
    )
    ax.add_artist(legend2)

    plt.text(500, 90, cgcid, fontsize=30, horizontalalignment="center")
    plt.text(500, 0, pulid, fontsize=30, horizontalalignment="center")
    plt.ylim(0, 100)
    plt.xlim(-100, 1000)
    plt.axis("off")
    ax.plot()
    plt.tight_layout(pad=0.01)
    cgcid = cgcid.replace("|", "_")  ### need to replace "|" to "_", because | is a special chara for system
    plt.savefig(f"syntenic.svg/{cgcid}.svg")
    plt.close()


def read_blast_result_cgc(filename):
    """
    Reads BLAST results from a file and organizes them into a dictionary.

    This function processes a file containing BLAST results, extracting relevant information and storing it in a structured format.

    Parameters
    ----------
        filename (str): Path to the file containing BLAST results.

    Returns
    -------
        dict: A dictionary where each key is a query identifier and the value is a list of `blastp_hit` objects representing the hits for that query.
    """
    querydict = {}
    for line in open(filename):
        lines = line.split()
        queryids = lines[0].split("|")
        queryid = queryids[0] + "|" + queryids[1]
        hitpulid = lines[1].split(":")[1]
        # lines[0] = lines[2]+"|"+lines[0] ### cgc id
        querydict.setdefault(queryid + ":" + hitpulid, []).append(blastp_hit(lines))
        # print (lines[0])
    return querydict


def syntenic_plot_allpairs(args):
    """
    Generates syntenic plots for all pairs of genomic regions.

    This function reads input data, processes it, and generates syntenic plots for each pair of genomic regions. It handles data reading, processing, and visualization.

    Parameters
    ----------
        args (argparse.Namespace): Command line arguments containing input file paths and parameters for the plot.

    Side Effects:
        - Creates a directory for storing the generated syntenic plots.
        - Generates and saves syntenic plots for each pair of genomic regions.
    """
    # error = 0
    cgcpul_blastp = read_blast_result_cgc(args.blastp)

    _, cgcid2gene, cgcid2geneid = read_UHGG_CGC_stanrdard_out(args.cgc)
    _, PULid2gene, PULid2geneid = read_PUL_cgcgff(args)

    os.makedirs("syntenic.svg", exist_ok=True)

    for line in open(args.input).readlines()[1:]:  ### for each pairs
        lines = line.rstrip().split("\t")
        cgc = lines[0]
        pul = lines[1]
        if not pul:  ### no homologous cgc found from dbCAN-PUL
            continue

        cgcpul = cgc + ":" + pul
        bed_cgc = cgcid2gene[cgc]
        bed_pul = PULid2gene[pul]

        starts1, ends1, strands1, types1 = Get_parameters_for_plot(bed_cgc)
        starts2, ends2, strands2, types2 = Get_parameters_for_plot(bed_pul)
        genes1 = cgcid2geneid[cgc]
        genes2 = PULid2geneid[pul]
        # print (cgc,pul)
        # print (genes1)
        # print (genes2)
        blocks = []
        for record in cgcpul_blastp[cgcpul]:  ### generate block information
            query = record.qseqid
            hit = record.sseqid
            cgc_proteinid = query.split("|")[2]
            pul_proteinid = hit.split(":")[3]
            if not pul_proteinid:
                pul_proteinid = hit.split(":")[2]
            # print (cgcpul,query,hit,cgc_proteinid,pul_proteinid,genes1,genes2)
            try:
                index1 = genes1.index(cgc_proteinid)
                index2 = genes2.index(pul_proteinid)
                blocks.append(f"{index1}-{index2}-{record.pident}")
            except ValueError as e:
                print(cgcpul, query, hit, cgc_proteinid, pul_proteinid, genes1, genes2)
                print(f"Error processing block data: {e}")
                continue
            # print (cgc_proteinid2gene[cgc_proteinid],pul_proteinid2gene[pul_proteinid])
        syntenic_plot(starts1, starts2, ends1, ends2, strands1, strands2, types1, types2, blocks, cgc, pul)


def Get_parameters_for_plot(CGC_stanrdard_list):
    """
    Extracts parameters required for plotting from a list of genomic data.

    This function processes a list of genomic data objects, extracting start and end positions, strand orientations, and gene types for plotting purposes.

    Parameters
    ----------
        CGC_stanrdard_list (list): A list of genomic data objects.

    Returns
    -------
        tuple: A tuple containing lists of starts, ends, strands, and types for the genomic data.
    """
    starts = []
    ends = []
    strands = []
    types = []
    for gene in CGC_stanrdard_list:
        starts.append(gene.Gene_Start)
        ends.append(gene.Gene_END)
        strands.append(gene.Strand)
        types.append(gene.Gene_Type)
    return starts, ends, strands, types


def attribution(desc, mark="DB"):
    """
    Extracts a specific attribute from a semicolon-separated string.

    Parameters
    ----------
        desc (str): The string containing various attributes separated by semicolons.
        mark (str, optional): The marker indicating the start of the desired attribute. Defaults to "DB".

    Returns
    -------
        str: The value of the specified attribute or "Other" if not found.
    """
    descs = desc.split(";")
    for des in descs:
        if des.startswith(mark):
            return des.split("=")[-1]
    return "Other"


def trim_proteinid(proteinid):
    """
    Trims a protein identifier to remove certain parts.

    Parameters
    ----------
        proteinid (str): The protein identifier to be trimmed.

    Returns
    -------
        str: The trimmed protein identifier.
    """
    clear_id = proteinid.split(".")
    return clear_id[0]
    # return ".".join(clear_id[0:-1])


def read_PUL_cgcout(filename="PUL.out"):
    """
    Reads genomic data from a file and organizes it into a dictionary.

    Parameters
    ----------
        filename (str, optional): Path to the file containing genomic data. Defaults to "PUL.out".

    Returns
    -------
        dict: A dictionary mapping protein identifiers to genomic data objects.
    """
    geneid2gene = {}
    for line in open(filename):
        if line.startswith("++++"):
            continue
        lines = line.rstrip("\n").split("\t")
        proteinid = trim_proteinid(lines[8])
        newline = [
            lines[4] + lines[0],
            lines[1],
            lines[5],
            proteinid,
            lines[5],
            lines[6],
            lines[9],
            attribution(lines[-1]),
        ]
        geneid2gene[proteinid] = CGC_stanrdard(newline)
    return geneid2gene


def read_cgcgff(filename, geneid2gene):
    """
    Reads genomic data from a GFF file and updates a dictionary with this information.

    Parameters
    ----------
        filename (str): Path to the GFF file.
        geneid2gene (dict): A dictionary to be updated with data from the GFF file.
    """
    for line in open(filename):
        lines = line.rstrip("\n").split("\t")
        proteinid = attribution(lines[-1], "ID")
        proteinid = trim_proteinid(proteinid)
        PULid = filename.split(".")[0].split("/")[-1]
        if lines[2] == "CDS":
            lines[2] = "other"
        # if PULid == "PUL0380":
        #    print (proteinid)
        newline = [PULid, lines[2], lines[0], proteinid, lines[3], lines[4], lines[6], attribution(lines[-1])]
        geneid2gene[PULid + ":" + proteinid] = CGC_stanrdard(newline)


def read_PUL_cgcgff(args):
    """
    Reads genomic data from multiple GFF files and organizes it into dictionaries.

    Parameters
    ----------
        args (argparse.Namespace): Command line arguments containing paths to GFF files and other parameters.

    Returns
    -------
        tuple: A tuple of three dictionaries containing genomic data organized in different ways.
    """
    PULidgeneid2gene = {}  ### "PUL00035:BT_3559"
    ### list all files

    dbCANPUL_data_folder = args.db_dir if args.db_dir.endswith("/") else args.db_dir + "/"
    if not dbCANPUL_data_folder.startswith("/"):  ### absolute path
        dbCANPUL_data_folder = os.path.abspath(dbCANPUL_data_folder) + "/dbCAN-PUL/"
    else:
        dbCANPUL_data_folder += "/dbCAN-PUL/"
    folder = os.listdir(dbCANPUL_data_folder)
    for filename in folder:
        if filename.startswith("PUL") and filename.endswith(".out"):  ###
            read_cgcgff(dbCANPUL_data_folder + filename + "/cgc.gff", PULidgeneid2gene)

    cgcid2gene = {}
    cgcid2geneid = {}
    for PULidgeneid in PULidgeneid2gene:
        gene = PULidgeneid2gene[PULidgeneid]
        cgcid2gene.setdefault(gene.CGCID, []).append(gene)
        cgcid2geneid.setdefault(gene.CGCID, []).append(gene.Protein_ID)

    return PULidgeneid2gene, cgcid2gene, cgcid2geneid


class CGC_stanrdard:
    """
    Represents a standard format for genomic data related to CGC (Conserved Gene Cluster).

    This class encapsulates the details of a genomic data entry, including information such as CGC ID, gene type, contig ID, protein ID, gene start and end positions, strand, and protein family.

    Attributes
    ----------
        CGCID (str): Identifier for the Conserved Gene Cluster.
        Gene_Type (str): Type of the gene.
        Contig_ID (str): Identifier for the contig.
        Protein_ID (str): Identifier for the protein.
        Gene_Start (int): Starting position of the gene.
        Gene_END (int): Ending position of the gene.
        Strand (str): Strand orientation of the gene.
        Protein_Family (str): Protein family classification.

    Args:
        lines (list): A list of strings representing the genomic data fields.
    """

    def __init__(self, lines):
        """
        Initializes a CGC_stanrdard instance with provided genomic data.

        Args:
            lines (list): A list of strings representing the genomic data fields.
        """
        self.CGCID = lines[0]
        self.Gene_Type = lines[1]
        self.Contig_ID = lines[2]
        self.Protein_ID = lines[3]
        self.Gene_Start = int(lines[4])
        self.Gene_END = int(lines[5])
        self.Strand = lines[6]
        self.Protein_Family = lines[7]

    def __repr__(self):
        """
        Returns a string representation of the CGC_stanrdard instance.

        Returns
        -------
            str: A tab-separated string of all attributes of the CGC_stanrdard instance.
        """
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])


def read_UHGG_CGC_stanrdard_out(filename):
    """
    Reads genomic data from a file and organizes it into structured formats.

    This function processes a file containing genomic data in a standard format, creating CGC_stanrdard objects for each line of data and organizing them into dictionaries by different keys.

    Parameters
    ----------
        filename (str): Path to the file containing genomic data.

    Returns
    -------
        tuple: A tuple containing three dictionaries:
            - geneid2gene: Maps protein IDs to CGC_stanrdard objects.
            - cgcid2gene: Maps CGC IDs to lists of CGC_stanrdard objects.
            - cgcid2geneid: Maps CGC IDs to lists of protein IDs.
    """
    geneid2gene = {}
    cgcid2gene = {}
    cgcid2geneid = {}

    for line in open(filename):
        if line[0:4] == "CGC#":
            continue
        lines = line.rstrip("\n").split("\t")
        lines[0] = lines[2] + "|" + lines[0]
        gene = CGC_stanrdard(lines)
        geneid2gene[gene.Protein_ID] = gene
        cgcid2gene.setdefault(gene.CGCID, []).append(gene)
        cgcid2geneid.setdefault(gene.CGCID, []).append(gene.Protein_ID)
    return geneid2gene, cgcid2gene, cgcid2geneid


def parse_argv():
    """
    Parses command line arguments for the syntenic plot script.

    This function uses argparse to define and parse command line arguments for generating syntenic plots. It specifies the necessary and optional arguments required to run the script.

    Returns
    -------
        argparse.Namespace: An object containing the parsed command line arguments. This includes:
            - function (str): The function to be used for analysis.
            - input (str): Path to the cgc_finder output file.
            - blastp (str): Path to the BLASTP result file for CGC.
            - cgc (str): Additional CGC-related argument.
            - pul (str): Additional PUL-related argument.
            - db_dir (str): Path to the database directory, defaults to 'db'.
    """
    parser = argparse.ArgumentParser(description="syntenic plot for two homologous CGC")
    parser.add_argument("function", help="what function will be used to analyze.")
    parser.add_argument("-i", "--input", help="cgc_finder output")
    parser.add_argument("-b", "--blastp", help="blastp result for cgc")
    parser.add_argument("--cgc")
    parser.add_argument("--pul")
    parser.add_argument("--db_dir", default="db", help="Database directory")
    args = parser.parse_args()
    return args


def main():
    """
    Main function to execute the syntenic plot script.

    This function calls `parse_argv` to parse command line arguments and then uses these arguments to generate syntenic plots for pairs of homologous CGC (Conserved Gene Clusters).

    """
    args = parse_argv()
    syntenic_plot_allpairs(args)


if __name__ == "__main__":
    args = parse_argv()
    # syntenic_plot(starts,starts1,ends,ends1,strands,strands1,Types,Types1,blocks)
    if args.function == "syntenic_plot":
        ## python3 /home/jinfang/libsvm_practise/syntenic.plot.py syntenic_plot -b PUL.pairs.blastp -i cgc.pul.hits.strict --cgc cgc_standard.out
        ## python3 syntenic.plot.py syntenic_plot -b PUL.pairs.blastp --cgc cgc_standard.out -i cgc.pul.hits.strict --pul PUL.out
        ## python3 syntenic.plot.py syntenic_plot -b PUL.pairs.blastp --cgc cgc_standard.out -i cgc.pul.hits.strict
        ## python3 /array1/www/dbCAN3/ty/syntenic.plot.py syntenic_plot -b blastp.out --cgc cgc_standard.out -i sub.prediction.out
        syntenic_plot_allpairs(args)
