import pysam
import argparse
import pandas as pd
import matplotlib.pyplot as plt

__author__ = "Sam Sims"
__version__ = "0.1.1"


class BamVisualiser:
    def __init__(self, bam):
        self.bam_path = bam
        self.bam_file = pysam.AlignmentFile(bam, "rb")
        self.ref_name = self.bam_file.get_reference_name(0)
        self.check_index()

    def check_index(self):
        if not self.bam_file.has_index():
            print("BAM file is not indexed, attempting to index...")
            pysam.index(self.bam_path)
            print("Done!")
            # Reopen the bam file once indexed
            self.bam_file = pysam.AlignmentFile(self.bam_path, "rb")

    def pileup(self, positions, min_depth):
        rows_list = [
            self.get_base_counts(int(position), min_depth) for position in positions
        ]
        pileup_df = pd.DataFrame(rows_list)
        print(pileup_df)
        return pileup_df

    def get_base_counts(self, position, min_depth):
        rows_dict = {"position": position, "A": 0, "T": 0, "C": 0, "G": 0}
        for col in self.bam_file.pileup(
            contig=self.ref_name,
            start=position,
            stop=position + 1,
            min_base_quality=0,
        ):
            bases = []
            for read in col.pileups:
                if col.pos == (position - 1):
                    if read.query_position is not None:
                        bases.append(read.alignment.query_sequence[read.query_position])
                    else:
                        bases.append(".")
                    rows_dict = {
                        "position": position,
                        "A": bases.count("A"),
                        "T": bases.count("T"),
                        "C": bases.count("C"),
                        "G": bases.count("G"),
                    }
        # check each key in dict and if the value is less than the min_depth, set it to 0
        for key in rows_dict:
            if rows_dict[key] < min_depth:
                rows_dict[key] = 0

        return rows_dict

    def plot_pileup(self, df, title, fig_width, individual):
        fig, ax = plt.subplots(figsize=(fig_width, 5))

        # Define the bar colors
        colors = ["#60935D", "#E63946", "#1B5299", "#F5BB00"]

        df.plot(x="position", kind="bar", stacked=True, color=colors, ax=ax)

        # Formatting
        ax.set_xlabel("Position")
        ax.set_ylabel("Count")
        # Move the legend to outside the plot
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
        ax.set_title(f"Ambiguous bases in BAM file at positions {title}")

        # Annotation to match bar order
        bar_order = df.columns[1:][::-1]
        for j, row in df.iterrows():
            if individual:
                for c in ax.containers:
                    column = c.get_label()
                    labels = [
                        f"{column}: {str(v.get_height())}" if v.get_height() > 0 else ""
                        for v in c
                    ]
                    print(labels)
                    ax.bar_label(
                        c,
                        labels=labels,
                        label_type="center",
                        size=8,
                    )
            else:
                total_count = row[bar_order].sum()
                annotation_text = "\n".join(
                    f"{column}: {count}"
                    for column, count in zip(bar_order, row[bar_order])
                    if count > 0
                )

                x = j
                y = total_count / 2

                ax.annotate(
                    annotation_text,
                    xy=(x, y),
                    xytext=(x, y),
                    ha="center",
                    va="center",
                    color="black",
                )

        plt.savefig("pileup.png")

    def pileup_percentages(self, pileup_df):
        total_counts = pileup_df[["A", "T", "C", "G"]].sum(axis=1)
        percentages_df = pileup_df[["A", "T", "C", "G"]].div(total_counts, axis=0) * 100
        pileup_df[["A", "T", "C", "G"]] = percentages_df.round(2)
        return pileup_df

    def visualise(self, args):
        if args.positions:
            positions = args.positions.split(",")
            pileup_df = self.pileup(positions, args.min_depth)
            title = ",".join(positions)
        elif args.start_pos and args.end_pos:
            positions = range(args.start_pos, args.end_pos + 1)
            pileup_df = self.pileup(positions, args.min_depth)
            title = f"{args.start_pos}-{args.end_pos}"

        if args.percentages:
            pileup_df = self.pileup_percentages(pileup_df)

        self.plot_pileup(pileup_df, title, args.fig_width, args.individual_annotations)

        # save the counts to a csv file
        if args.save_counts:
            pileup_df.to_csv(args.save_counts, index=False)

    def test(self):
        data = {
            "position": [140, 141, 142, 143, 145, 145],
            "A": [140, 500, 100, 347, 0, 300],
            "T": [300, 0, 321, 45, 55, 0],
            "C": [0, 0, 15, 20, 25, 0],
            "G": [0, 0, 40, 50, 60, 0],
        }
        test_df = pd.DataFrame(data)
        self.plot_pileup(test_df, "140-145", 20, True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualise ambiguous posistions in a BAM file"
    )
    parser.add_argument("-b", "--bam", help="BAM file to visualise", required=True)
    parser.add_argument("--positions", help="Positions to visualise")
    parser.add_argument("--start_pos", help="Start position to visualise", type=int)
    parser.add_argument("--end_pos", help="End position to visualise", type=int)
    parser.add_argument(
        "--percentages", help="Show percentages instead of counts", action="store_true"
    )
    parser.add_argument(
        "--min_depth", help="Minimum depth to visualise", default=0, type=int
    )
    parser.add_argument("--save_counts", help="Save counts to csv file")
    parser.add_argument(
        "--fig_width", help="Adjust width of figure", default=20, type=int
    )
    parser.add_argument(
        "--individual_annotations",
        help="Show individual annotations",
        action="store_true",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    bam_vis = BamVisualiser(args.bam)
    bam_vis.visualise(args)
    #bam_vis.test()


if __name__ == "__main__":
    main()
