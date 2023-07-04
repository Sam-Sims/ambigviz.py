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

    def pileup(self, positions, min_depth):
        rows_list = [
            self.get_base_counts(int(position), min_depth) for position in positions
        ]
        pileup_df = pd.DataFrame(rows_list)
        return pileup_df

    def get_base_counts(self, position, min_depth):
        base_counts_dict = {"position": position, "A": 0, "T": 0, "C": 0, "G": 0}
        pileup_columns = self.bam_file.count_coverage(
            contig=self.ref_name,
            start=position - 1,
            stop=position,
            quality_threshold=0,
        )
        for base, counts in zip("ATCG", pileup_columns):
            base_counts_dict[base] = sum(counts)

        for key in base_counts_dict:
            if base_counts_dict[key] < min_depth:
                base_counts_dict[key] = 0

        return base_counts_dict

    def plot_pileup(self, df, title, fig_width, individual):
        fig, ax = plt.subplots(figsize=(fig_width, 5))
        colors = ["#60935D", "#E63946", "#1B5299", "#F5BB00"]

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

        pileup_df = self.pileup(positions, args.min_depth)

        if args.percentages:
            pileup_df = self.pileup_percentages(pileup_df)

        figure_to_plot = self.plot_pileup(
            pileup_df, title, args.fig_width, args.individual_annotations
        )

        # Outputs
        figure_to_plot.savefig(args.output, bbox_inches="tight")
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
    return parser.parse_args()
# fmt: on


def main():
    args = parse_args()
    bam_vis = BamVisualiser(args.bam)
    bam_vis.visualise(args)


if __name__ == "__main__":
    main()
