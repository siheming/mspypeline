from MQInitializer import MQInitializer
from MQPlots import MQPlots
import argparse
import tkinter as tk
from tkinter import filedialog
import logging


class MQGUI(tk.Tk):
    def __init__(self, no_gui, loglevel=logging.DEBUG):
        super().__init__()

        self._start_dir = ""
        self.yaml_options = ["default"]
        self.plot_names = ["plot_detection_counts", "plot_number_of_detected_proteins", "plot_intensity_histograms",
            "plot_relative_std", "plot_rank", "plot_pathway_analysis", "plot_pathway_proportions",
            "plot_scatter_replicates", "plot_experiment_comparison", "plot_go_analysis", "plot_venn_results"]
        self.called_update = False

        self.mqinit = MQInitializer("", True, "default", loglevel=loglevel)

        if no_gui:
            self.withdraw()
        else:
            self.make_layout()

    @property
    def start_dir(self):
        return self._start_dir

    @start_dir.setter
    def start_dir(self, start_dir):
        self._start_dir = start_dir
        self.mqinit.start_dir = self.start_dir
        if self.mqinit.has_yml_file():
            self.mqinit.file_path_yaml = "file"
            self.yaml_text.set("file")
            self.yaml_options = ["default", "file"]
            #print(self.mqinit.configs)
            # next steps should also be done hy what the update button does?
        else:
            self.mqinit.file_path_yaml = "default"
            self.yaml_text.set("default")
            self.yaml_options = ["default"]
        self.yaml_button["menu"].delete(0, "end")
        for op in self.yaml_options:
            self.yaml_button['menu'].add_command(label=op, command=tk._setit(self.yaml_text, op))
        self.called_update = False

    def dir_setter(self, *args):
        self.start_dir = self.dir_text.get()

    def yaml_path_setter(self, *args):
        if self.yaml_text.get() == "file":
            self.mqinit.file_path_yaml = "file"
        elif self.yaml_text.get() == "default":
            self.mqinit.file_path_yaml = "default"
        else:
            raise ValueError(f"Unkown setting for yaml text {self.yaml_text.get()}")
        # TODO this should be dynamic
        for plot_name in self.plot_names:
            plot_intensity = plot_name + "_intensity"
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("_var", "_int")
            getattr(self, var_name).set(self.mqinit.configs.get(plot_intensity, "raw"))
            getattr(self, int_name).set(self.mqinit.configs.get(plot_name, False))
        if self.mqinit.configs.get("experiments"):
            self.experiments_list.delete(0, "end")
            for op in self.mqinit.configs.get("experiments"):
                self.experiments_list.insert("end", op)
        else:
            self.experiments_list.delete(0, "end")

    def update_button(self):
        for plot_name in self.plot_names:
            plot_intensity = plot_name + "_intensity"
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("_var", "_int")
            self.mqinit.configs.update({plot_intensity: getattr(self, var_name).get()})
            self.mqinit.configs.update({plot_name: bool(getattr(self, int_name).get())})
        self.mqinit.init_config()
        self.mqinit.prepare_stuff()
        if self.mqinit.configs.get("experiments"):
            self.experiments_list.delete(0, "end")
            for op in self.mqinit.configs.get("experiments"):
                self.experiments_list.insert("end", op)
        else:
            self.experiments_list.delete(0, "end")
        self.called_update = True

    def start_button(self):
        if not self.called_update:
            self.update_button()
        mqplots = MQPlots.from_MQInitializer(self.mqinit)
        mqplots.create_results()

    def make_layout(self):
        self.title("MaxQuant Analyzer Pipeline")

        path_label = tk.Label(self, text="Dir to analyze").grid(row=0, column=0)

        yaml_label = tk.Label(self, text="Yaml file").grid(row=0, column=1)

        self.dir_text = tk.StringVar(value=self.start_dir)
        self.dir_text.trace("w", self.dir_setter)
        dir_button = tk.Button(self, textvariable=self.dir_text,
                               command=lambda: browsefunc(filedialog.askdirectory, self.dir_text, fn_params={
                                   "title": "Please select a directory with MaxQuant result files"}))
        dir_button.grid(row=1, column=0)

        self.yaml_text = tk.StringVar()
        self.yaml_text.trace("w", self.yaml_path_setter)
        self.yaml_button = tk.OptionMenu(self, self.yaml_text, *self.yaml_options)
        self.yaml_button.grid(row=1, column=1)

        replicate_var = tk.IntVar(value=1)
        replicate_button = tk.Checkbutton(self, text="Does the file have replicates?", variable=replicate_var).grid(
            row=2, column=0)

        go_proteins_label = tk.Label(self, text="Go analysis proteins").grid(row=3, column=0)

        experiments_label = tk.Label(self, text="Experiments").grid(row=3, column=1)

        go_proteins_list = tk.Listbox(self, selectmode="multiple", height=3)
        for x in self.mqinit.possible_compartiments:
            go_proteins_list.insert("end", x)

        go_proteins_list.grid(row=4, column=0)

        self.experiments_list = tk.Listbox(self, height=3)
        self.experiments_list.grid(row=4, column=1)

        plot_label = tk.Label(self, text="Which plots should be created").grid(row=5, column=0)

        intensity_label = tk.Label(self, text="Intensity").grid(row=5, column=1)

        self.detection_counts_int, self.detection_counts_var = plot_row(self, "Detection counts", 6 , "raw")
        self.number_of_detected_proteins_int, self.number_of_detected_proteins_var = plot_row(self, "Number of detected proteins", 7, "raw")
        self.intensity_histograms_int, self.intensity_histograms_var = plot_row(self, "Intensity histogram", 8, "raw")
        self.relative_std_int, self.relative_std_var = plot_row(self, "Relative std", 9, "raw")
        self.rank_int, self.rank_var = plot_row(self, "Rank", 10, "raw")
        self.pathway_analysis_int, self.pathway_analysis_var = plot_row(self, "Pathway Analysis", 11, "raw")
        self.pathway_proportions_int, self.pathway_proportions_var = plot_row(self, "Pathway proportions", 12, "raw")
        self.scatter_replicates_int, self.scatter_replicates_var = plot_row(self, "Scatter replicates", 13, "raw")
        self.experiment_comparison_int, self.experiment_comparison_var = plot_row(self, "Experiment comparison", 15, "raw")
        self.go_analysis_int, self.go_analysis_var = plot_row(self, "Go analysis", 16, "raw")
        self.venn_results_int, self.venn_results_var = plot_row(self, "Venn diagrams", 17, "raw")

        update_button = tk.Button(self, text="Update", command=lambda: self.update_button())
        update_button.grid(row=18, column=0)

        start_button = tk.Button(self, text="Start",
                                 command=lambda: self.start_button())
        start_button.grid(row=18, column=1)


def main():
    # create arg parser to handle different options
    parser = argparse.ArgumentParser(description="A pipeline to analyze result files from a MaxQuant report. "
                                                 "The required result files are in the txt directory.")
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
        help="Path to yml file which should be used for analysis, or 'default'. "
             "If not set the program will open a prompt and ask to select one, "
             "but if canceled the default yml file will be used."
    )
    parser.add_argument(
        "--loglevel",
        dest="loglevel",
        action="store",
        default=logging.WARNING,
        help="Logging level of analysis. Should be from options (lowest to highest): DEBUG < INFO < WARNING < ERROR. "
             "The higher the logging level the fewer messages are shown. Default: WARNING"
    )
    parser.add_argument(
        "--has-replicates",
        dest="has_replicates",
        action="store",
        default=None,
        help="If you have replicates of each experiment enter y else enter n"
             "Replicates need to be enumerated in a way that numbers are added at the end of the name"
             "If no replicates are in data set no venn diagrams will be generated"
    )
    parser.add_argument(
        "--no-gui",
        default=False,
        const=True,
        nargs="?",
        help="specify this if no gui should be opened"
    )
    args = parser.parse_args()
    args_dict = vars(args)

    # determine logging level
    try:
        loglevel = getattr(logging, args.loglevel.upper())
    except AttributeError:
        try:
            loglevel = int(args.loglevel)
        except ValueError:
            loglevel = logging.DEBUG

    # disable most ui things
    gui = MQGUI(args.no_gui, loglevel=loglevel)
    if args.no_gui:
        # if no dir was specified ask for one
        if args.dir is None:
            start_dir = filedialog.askdirectory(title="Please select a directory with MaxQuant result files")
            if not start_dir:
                raise ValueError("Please select a directory")
        else:
            start_dir = args.dir

        bool_dict = {"yes": True, "y": True, "true": True, "no": False, "n": False, "false": False}
        # ask if the file has replicates with yes as default
        if args.has_replicates is None:
            has_replicates = input("Please specify if you have replicates in your data [Y/n]: ")
            if not has_replicates:
                has_replicates = True
            else:
                has_replicates = has_replicates.lower()
                if has_replicates not in bool_dict:
                    raise ValueError(f"Please use a valid answer({', '.join(bool_dict.keys())})")
                has_replicates = bool_dict[has_replicates]
        else:
            has_replicates = bool_dict[args.has_replicates.lower()]

        # create initializer which reads all required files
        mqinit = MQInitializer(start_dir, has_replicates, args.yml_file, loglevel=loglevel)
        mqinit.init_config()
        mqinit.prepare_stuff()
        # create plotter from initializer
        mqplots = MQPlots.from_MQInitializer(mqinit, loglevel=loglevel)
        # create all plots and other results
        mqplots.create_results()
    else:
        if args.dir is not None:
            gui.dir_text.set(args.dir)

        gui.mainloop()

def browsefunc(fn, var, fn_params: dict = None):
    if fn_params is None:
        fn_params = {}
    filename = fn(**fn_params)
    var.set(filename)

def plot_row(frame, text: str, row, intensity_default: str):
    int_var = tk.IntVar(value=1)
    plot_button = tk.Checkbutton(frame, text=text, variable=int_var).grid(row=row, column=0)

    intensity_var = tk.StringVar(value=intensity_default)

    intensity_menu = tk.OptionMenu(frame, intensity_var, "lfq", "raw").grid(row=row, column=1)
    return int_var, intensity_var


def plot_tk_buttons():
    pass

if __name__ == "__main__":
    main()