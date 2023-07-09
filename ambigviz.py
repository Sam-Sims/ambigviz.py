import pysam
import argparse
import pandas as pd
import matplotlib.pyplot as plt

__author__ = "Sam Sims"
__version__ = "0.1.2"


class BamVisualiser:
    def __init__(self, bam):
        self.bam_path = bam
        self.bam_file = pysam.AlignmentFile(bam, "rb")
        self.ref_name = self.bam_file.get_reference_name(0)
        self.ref_length = self.bam_file.get_reference_length(self.ref_name)
        self.check_index()

    def check_index(self):
        if not self.bam_file.has_index():
            print("BAM file is not indexed, attempting to index...")
            # try index, else raise exception
            try:
                pysam.index(self.bam_path)
                print("Done!")
                # Reopen the bam file once indexed
                self.bam_file = pysam.AlignmentFile(self.bam_path, "rb")
            except Exception as e:
                print(f"Could not index BAM file: {e}")
                exit(1)

    def pileup(self, positions, min_depth, quality_threshold):
        rows_list = [
            self.get_base_counts(int(position), min_depth, quality_threshold)
            for position in positions
        ]
        pileup_df = pd.DataFrame(rows_list)
        return pileup_df

    def get_base_counts(self, position, min_depth, quality_threshold):
        base_counts_dict = {"position": position, "A": 0, "T": 0, "C": 0, "G": 0}
        pileup_columns = self.bam_file.count_coverage(
            contig=self.ref_name,
            start=position - 1,
            stop=position,
            quality_threshold=quality_threshold,
        )
        for base, counts in zip("ATCG", pileup_columns):
            base_counts_dict[base] = sum(counts)

        for key in base_counts_dict:
            if base_counts_dict[key] < min_depth:
                base_counts_dict[key] = 0

        return base_counts_dict

    def plot_pileup(self, df, title, fig_width, individual):
        fig, ax = plt.subplots(figsize=(fig_width, 5))
        colors = ["#60935D", "#E63946", "#1B5299", "#F5BB00", "#F7D1CD", "#FBBA72"]

        df.plot(x="position", kind="bar", stacked=True, color=colors, ax=ax)

        # Formatting
        ax.set(xlabel="Position", ylabel="Count")
        ax.legend(bbox_to_anchor=(1.05, 1))
        ax.set_title(f"Ambiguous bases in BAM file at positions {title}")

        # Annotation
        bar_order = df.columns[1:][::-1]
        for idx, row in df.iterrows():
            if individual:
                for container in ax.containers:
                    labels = [
                        f"{container}: {bar.get_height()}"
                        if bar.get_height() > 0
                        else ""
                        for container, bar in zip(bar_order, container)
                    ]
                    ax.bar_label(container, labels=labels, label_type="center", size=8)
            else:
                total_count = row[bar_order].sum()
                annotation_text = "\n".join(
                    f"{container}: {count}"
                    for container, count in zip(bar_order, row[bar_order])
                    if count > 0
                )
                x = idx
                y = total_count / 2
                # make text white
                ax.annotate(annotation_text, xy=(x, y), ha="center", va="center")

        return fig

    def pileup_percentages(self, pileup_df):
        total_counts = pileup_df[["A", "T", "C", "G"]].sum(axis=1)
        percentages_df = pileup_df[["A", "T", "C", "G"]].div(total_counts, axis=0) * 100
        pileup_df[["A", "T", "C", "G"]] = percentages_df.round(2)
        return pileup_df

    def visualise(self, args):
        if args.positions:
            positions = args.positions.split(",")
            title = ",".join(positions)
        elif args.start_pos and args.end_pos:
            positions = range(args.start_pos, args.end_pos + 1)
            title = f"{args.start_pos}-{args.end_pos}"

        pileup_df = self.pileup(positions, args.min_depth, args.quality_threshold)
        if args.indels:
            # call count_indels which returns dict of Insertions: counts and Deletions: counts for each position - add to pileup_df as 2 new columns
            indel_counts = self.count_indels()
            pileup_df["Insertions"] = pileup_df["position"].map(
                lambda x: indel_counts[x]["Insertion"]
            )
            pileup_df["Deletions"] = pileup_df["position"].map(
                lambda x: indel_counts[x]["Deletion"]
            )
        print(pileup_df)

        if args.percentages:
            pileup_df = self.pileup_percentages(pileup_df)

        figure_to_plot = self.plot_pileup(
            pileup_df, title, args.fig_width, args.individual_annotations
        )

        # Outputs
        figure_to_plot.savefig(args.output, bbox_inches="tight")
        if args.save_counts:
            pileup_df.to_csv(args.save_counts, index=False)

    def calculate_all(self, args):
        positions = range(1, self.ref_length + 1)

        # for each position in positions, create a dict with position and base counts set to 0
        base_counts_dict = {
            position: {"A": 0, "T": 0, "C": 0, "G": 0, "Insertion": 0, "Deletion": 0}
            for position in positions
        }

        # count_coverage for all positions
        pileup_columns = self.bam_file.count_coverage(
            contig=self.ref_name,
            quality_threshold=8,
        )

        # pileup_columns contains 4 tuples, each tuple contains a list of counts for each base in the order A, C, G, T
        # for each position, add the counts of each base to the base_counts_dict
        for pos in positions:
            base_counts_dict[pos]["A"] += pileup_columns[0][pos - 1]
            base_counts_dict[pos]["C"] += pileup_columns[1][pos - 1]
            base_counts_dict[pos]["G"] += pileup_columns[2][pos - 1]
            base_counts_dict[pos]["T"] += pileup_columns[3][pos - 1]

        # convert base_counts_dict to a dataframe, with position column
        base_counts_df = pd.DataFrame.from_dict(base_counts_dict, orient="index")
        base_counts_df.index.name = "position"
        base_counts_df.reset_index(inplace=True)

        if args.indels:
            # call count_indels which returns dict of Insertions: counts and Deletions: counts for each position - add to pileup_df as 2 new columns
            indel_counts = self.count_indels()

            # Add "Insertion" and "Deletion" columns to the dataframe
            base_counts_df["Insertion"] = base_counts_df["position"].map(
                lambda x: indel_counts[x]["Insertion"]
            )
            base_counts_df["Deletion"] = base_counts_df["position"].map(
                lambda x: indel_counts[x]["Deletion"]
            )
        # remove positions with 0 coverage
        base_counts_df = base_counts_df[
            base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            > 0
        ]

        # Calculate insertion and deletion percentages
        base_counts_df["Insertion_percent"] = (
            base_counts_df["Insertion"]
            / base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            * 100
        )
        base_counts_df["Deletion_percent"] = (
            base_counts_df["Deletion"]
            / base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            * 100
        )

        # Calculate A, T, C, G percentages
        base_counts_df["A_percent"] = (
            base_counts_df["A"]
            / base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            * 100
        )
        base_counts_df["T_percent"] = (
            base_counts_df["T"]
            / base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            * 100
        )
        base_counts_df["C_percent"] = (
            base_counts_df["C"]
            / base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            * 100
        )
        base_counts_df["G_percent"] = (
            base_counts_df["G"]
            / base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            * 100
        )

        # set any values less than read_fraction to 0
        base_counts_df[base_counts_df < args.read_fraction] = 0

        # remove any rows where only 1 base has a value greater than 0
        base_counts_df = base_counts_df[
            (
                base_counts_df[
                    [
                        "A_percent",
                        "T_percent",
                        "C_percent",
                        "G_percent",
                        "Insertion_percent",
                        "Deletion_percent",
                    ]
                ]
                != 0
            ).sum(axis=1)
            > 1
        ]

        # Remove any rows where all values are 0
        base_counts_df = base_counts_df[
            base_counts_df[["A", "T", "C", "G", "Insertion", "Deletion"]].sum(axis=1)
            > 0
        ]

        base_counts_df.to_csv("mixed_bases.csv", index=False)

        # Remove the percent columns
        base_counts_df.drop(
            columns=[
                "A_percent",
                "T_percent",
                "C_percent",
                "G_percent",
                "Insertion_percent",
                "Deletion_percent",
            ],
            inplace=True,
        )

        # Plot the mixed bases using plot_pileup
        figure_to_plot = self.plot_pileup(
            base_counts_df,
            f"{self.ref_name} - Mixed Bases",
            args.fig_width,
            args.individual_annotations,
        )
        figure_to_plot.savefig("pileup.png", bbox_inches="tight")

    def count_indels(self):
        # Open BAM file for reading
        positions = range(1, self.ref_length + 1)

        # for each position in positions, create a dict with position and base counts set to 0
        indel_count_dict = {
            position: {
                "Deletion": 0,
                "Insertion": 0,
            }
            for position in positions
        }

        # Iterate through each position in the BAM file
        for pileupcolumn in self.bam_file.pileup(
            contig=self.ref_name, min_base_quality=0
        ):
            position = pileupcolumn.pos

            pileupreads = pileupcolumn.pileups

            # Iterate through each pileup read
            for pileupread in pileupreads:
                if pileupread.is_del:
                    indel_count_dict[position + 1]["Deletion"] += 1
                elif pileupread.indel > 0 and not pileupread.is_del:
                    indel_count_dict[position + 2]["Insertion"] += 1
        return indel_count_dict

    def test(self):
        data = {
            "position": [140, 141, 142, 143, 145, 145],
            "A": [140, 500, 100, 347, 0, 300],
            "T": [300, 0, 321, 45, 55, 0],
            "C": [0, 0, 15, 20, 25, 0],
            "G": [0, 0, 40, 50, 60, 0],
        }
        test_df = pd.DataFrame(data)
        self.plot_pileup(test_df, "140-145", 20, False)


# fmt: off
def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualise ambiguous posistions in a BAM file"
    )
    parser.add_argument("-b", "--bam", 
                        help="BAM file to visualise", 
                        required=True)
    
    parser.add_argument("-o", "--output",
                        help="Output file name",
                        default="pileup.png")
    
    parser.add_argument("--positions", 
                        help="Positions to visualise")
    
    parser.add_argument("--start_pos", 
                        help="Start position to visualise", 
                        type=int)
    
    parser.add_argument("--end_pos", 
                        help="End position to visualise", 
                        type=int)
    
    parser.add_argument("--percentages",
                        help="Show percentages instead of counts", 
                        action="store_true"
    )
    parser.add_argument("--min_depth", 
                        help="Minimum depth to visualise", 
                        default=0, 
                        type=int
    )
    parser.add_argument("--save_counts", 
                        help="Save counts to csv file")
    
    parser.add_argument("--fig_width", 
                        help="Adjust width of figure", 
                        default=20, 
                        type=int
    )
    parser.add_argument("--individual_annotations",
                        help="Show individual annotations",
                        action="store_true",
    )
    parser.add_argument("--quality_threshold",
                        help="Quality threshold for pileup",
                        default=0,
                        type=int,
    )
    parser.add_argument("--calculate_all",
                        help="Save list of all positions with mixed bases to csv file.",
                        action="store_true",
    )
    parser.add_argument("--read_fraction",
                        help="Minimum fraction of reads that must contain the base for it to be included in the mixed base csv file.",
                        default=0.0,
                        type=float,
    )
    parser.add_argument("--indels",
                        help="Include indels in the pileup.",
                        action="store_true",
    )
    return parser.parse_args()
# fmt: on


def main():
    args = parse_args()
    bam_vis = BamVisualiser(args.bam)
    if args.calculate_all:
        bam_vis.calculate_all(args)
    else:
        bam_vis.visualise(args)


if __name__ == "__main__":
    main()
