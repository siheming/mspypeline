from MQInitializer import MQInitializer
from MQPlots import MQPlots
import argparse
import tkinter as tk
from tkinter import filedialog


def main():
    # create arg parser to handle different options
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        '--dir',
        dest='dir',
        action='store',
        help="Path to directory of analysis which should be analyzed."
             "If not set the program will open a prompt and ask to select one."
    )
    parser.add_argument(
        '--yml-file',
        dest='yml_file',
        action='store',
        help="Path to yml file which should be used for analysis, or 'default'."
             "If not set the program will open a prompt and ask to select one,"
             "but if canceled the default yml file will be used."
    )
    args = parser.parse_args()
    args_dict = vars(args)
    # args.dir = "/media/b200-linux/Elements/max_quant_txt_results/xaiomeng_combined_buffers_no_isoform/"
    # print(args)
    # disable most ui things
    root = tk.Tk()
    root.withdraw()
    # if no dir was specified ask for one
    if args.dir is None:
        start_dir = filedialog.askdirectory()
        if not start_dir:
            raise ValueError("Please select a directory")
    else:
        start_dir = args.dir
    mqinit = MQInitializer(start_dir, args.yml_file)
    mqplots = MQPlots.from_MQInitializer(mqinit)
    mqplots.create_results()


if __name__ == "__main__":
    main()