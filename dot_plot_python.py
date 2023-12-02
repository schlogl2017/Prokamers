import argparse
import numpy as np
import matplotlib.pyplot as plt
# usage: dot_plot_python.py [-h] w s seqA seqB title output
# dot_plot_python.py: error: the following arguments are required: w, s, seqA, # seqB, title, output

def dotplot(seqA, seqB, w, s):
    # Initialize the dotplot matrix with zeros
    dp = np.zeros((len(seqA), len(seqB)), dtype=int)
    # Iterate over all positions in seqA and seqB
    for i in range(len(seqA)):
        for j in range(len(seqB)):
            # Compute the window around position i in seqA and position j in seqB
            windowA = seqA[max(i-w, 0):min(i+w+1, len(seqA))]
            windowB = seqB[max(j-w, 0):min(j+w+1, len(seqB))]
            # Count the number of matching symbols in the window
            matches = sum([1 for x, y in zip(windowA, windowB) if x == y])
            # Set the matrix element to 1 if the number of matches is at least s
            if matches >= s:
                dp[i,j] = 1
    return dp


def dotplot2Ascii(dp, seqA, seqB, heading, filename):
    # Open the output file
    with open(filename, 'w') as f:
        # Write the heading
        f.write(heading + "\n")
        # Write the x-axis labels
        f.write(" " * 5 + seqB + "\n")
        # Write the dotplot matrix
        for i in range(len(seqA)):
            f.write(seqA[i] + " | ")
            for j in range(len(seqB)):
                if dp[i][j]:
                    f.write("*")
                else:
                    f.write(".")
            f.write("\n")


def dotplot2Graphics(dp, labelA, labelB, heading, filename):
    # create a new figure
    fig, ax = plt.subplots()

    # plot the dots using a scatter plot
    rows, cols = np.where(dp)
    ax.scatter(cols, rows, marker='.', color='black')

    # set the labels and title
    ax.set_xlabel(labelB)
    ax.set_ylabel(labelA)
    ax.set_title(heading)

    # set the tick positions and labels
    xticks = np.arange(0.5, dp.shape[1], 1)
    yticks = np.arange(0.5, dp.shape[0], 1)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(np.arange(1, dp.shape[1] + 1))
    ax.set_yticklabels(np.arange(1, dp.shape[0] + 1)[::-1])

    # save the figure to a file and display it on screen
    plt.savefig(filename)
    plt.show()


def main():
    # Define the command line arguments
    parser = argparse.ArgumentParser(description="Generate a dotplot from two sequences in FASTA format.")
    parser.add_argument("w", type=int, help="the window size")
    parser.add_argument("s", type=int, help="the stringency")
    parser.add_argument("seqA", help="the filename of the first sequence in FASTA format")
    parser.add_argument("seqB", help="the filename of the second sequence in FASTA format")
    parser.add_argument("title", help="the title of the dotplot")
    parser.add_argument("output", help="the output filename for the dotplot")

    # Parse the command line arguments
    args = parser.parse_args()

    # Read the sequences from the input files
    with open(args.seqA) as fileA, open(args.seqB) as fileB:
        seqA = "".join([line.strip() for line in fileA if not line.startswith(">")])
        seqB = "".join([line.strip() for line in fileB if not line.startswith(">")])

    # Generate the dotplot matrix
    dp = dotplot(seqA, seqB, args.w, args.s)

    # Generate the ASCII dotplot
    dotplot2Ascii(dp, seqA, seqB, args.title, args.output + ".txt")

    # Generate the graphical dotplot
    dotplot2Graphics(dp, args.seqA, args.seqB, args.title, args.output + ".png")


if __name__ == "__main__":
    main()







