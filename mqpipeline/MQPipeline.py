from mqpipeline.MQInitializer import MQInitializer
from mqpipeline.MQPlots import MQPlots
import argparse
import tkinter as tk
from tkinter import filedialog
import logging


class MQUI(tk.Tk):
    def __init__(self, file_dir, yml_file="default", gui=False, loglevel=logging.DEBUG, configs: dict = None):
        super().__init__()

        # determine logging level
        try:
            loglevel = getattr(logging, loglevel.upper())
        except AttributeError:
            try:
                loglevel = int(loglevel)
            except ValueError:
                loglevel = logging.DEBUG

        base_config = {
            "has_replicates": False,
            "has_groups": False,
        }
        if configs is None:
            configs = {}
        base_config.update(**configs)
        configs = base_config
        self.yaml_options = ["default"]

        self.number_of_plots = 0

        if not gui:
            # get all necessary data, start the analysis and quit
            self.withdraw()
            # create initializer which reads all required files
            mqinit = MQInitializer(file_dir, yml_file, loglevel=loglevel)
            mqinit.init_config()
            mqinit.configs.update(configs)
            mqinit.prepare_stuff()
            # create plotter from initializer
            mqplots = MQPlots.from_MQInitializer(mqinit)
            # create all plots and other results
            mqplots.create_results()
        else:
            self.mqinit = MQInitializer("", yml_file, loglevel=loglevel)
            self.make_layout()
            if file_dir:
                self.dir_text.set(file_dir)
            self.mqinit.configs.update(configs)
            self.mainloop()

    def dir_setter(self, *args):
        self.mqinit.start_dir = self.dir_text.get()
        # write back possible modifications
        self.dir_text.set(self.mqinit.start_dir)
        if self.mqinit.has_yml_file():
            self.mqinit.file_path_yaml = "file"
            self.yaml_text.set("file")
            self.yaml_options = ["default", "file"]
            # next steps should also be done hy what the update button does?
        else:
            self.mqinit.file_path_yaml = "default"
            self.yaml_text.set("default")
            self.yaml_options = ["default"]
        self.yaml_button["menu"].delete(0, "end")
        for op in self.yaml_options:
            self.yaml_button["menu"].add_command(label=op, command=tk._setit(self.yaml_text, op))

    def update_listboxes(self):
        # delete all experiments then add from file
        self.experiments_list.delete(0, "end")
        if self.mqinit.configs.get("experiments"):
            for op in self.mqinit.configs.get("experiments"):
                self.experiments_list.insert("end", op)
        # clear selection then select from configs
        for i, pathway in enumerate(MQInitializer.possible_pathways):
            self.pathway_list.select_clear(i)
        if self.mqinit.configs.get("pathways"):
            for pathway in self.mqinit.configs.get("pathways"):
                self.pathway_list.select_set(MQInitializer.possible_pathways.index(pathway))
        # clear selection then select from configs
        for i, go in enumerate(MQInitializer.possible_gos):
            self.go_proteins_list.select_clear(i)
        if self.mqinit.configs.get("go_terms"):
            for go in self.mqinit.configs.get("go_terms"):
                self.go_proteins_list.select_set(MQInitializer.possible_gos.index(go))

    def yaml_path_setter(self, *args):
        self.mqinit.file_path_yaml = self.yaml_text.get()
        level_names = self.mqinit.configs.get("level_names", [])
        level_names = {i: name for i, name in enumerate(level_names)}
        levels = self.mqinit.configs.get("levels", 3)
        for plot_name in MQPlots.possible_plots:
            plot_settings_name = plot_name + "_settings"
            plot_settings = self.mqinit.configs.get(plot_settings_name, {})
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("var", "int")
            levels_name = var_name.replace("var", "levels")
            getattr(self, var_name).set(plot_settings.get("intensity", "raw"))
            getattr(self, int_name).set(plot_settings.get("create_plot", False))
            selected_levels = getattr(self, levels_name)
            selected_levels.delete(0, "end")
            for level in range(levels):
                selected_levels.insert("end", level_names.get(level, level))
            for level in plot_settings.get("levels", []):
                selected_levels.select_set(level)
        self.replicate_var.set(self.mqinit.configs.get("has_replicates", True))
        self.update_listboxes()

    def update_button(self):
        for plot_name in MQPlots.possible_plots:
            plot_settings = plot_name + "_settings"
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("var", "int")
            levels_name = var_name.replace("var", "levels")
            selected_settings = {
                "create_plot": bool(getattr(self, int_name).get()),
                "df_to_use": getattr(self, var_name).get(),
                "levels": [int(x) for x in getattr(self, levels_name).curselection()]
            }
            self.mqinit.configs.update({plot_settings: selected_settings})
        gos = self.go_proteins_list.curselection()
        gos = [MQInitializer.possible_gos[int(go)] for go in gos]
        pathways = self.pathway_list.curselection()
        pathways = [MQInitializer.possible_pathways[int(pathway)] for pathway in pathways]
        self.mqinit.configs["go_terms"] = gos
        self.mqinit.configs["pathways"] = pathways
        self.mqinit.configs["has_replicates"] = bool(self.replicate_var.get())
        self.mqinit.init_config()
        self.mqinit.prepare_stuff()
        self.update_listboxes()

    def start_button(self):
        self.running_text.set("Creating Plots")
        self.update()
        self.update_button()
        mqplots = MQPlots.from_MQInitializer(self.mqinit)
        mqplots.create_results()
        self.running_text.set("Please press Start")

    def make_layout(self):
        self.title("MaxQuant Analyzer Pipeline")

        path_label = tk.Label(self, text="Dir to analyze").grid(row=0, column=0)

        yaml_label = tk.Label(self, text="Yaml file").grid(row=0, column=1)

        self.dir_text = tk.StringVar(value=self.mqinit.start_dir)
        dir_button = tk.Button(self, textvariable=self.dir_text,
                               command=lambda: browsefunc(filedialog.askdirectory, self.dir_text, fn_params={
                                   "title": "Please select a directory with MaxQuant result files"}))
        dir_button.grid(row=1, column=0)

        self.yaml_text = tk.StringVar()
        self.yaml_button = tk.OptionMenu(self, self.yaml_text, *self.yaml_options)
        self.yaml_button.grid(row=1, column=1)

        self.replicate_var = tk.IntVar(value=1)
        replicate_button = tk.Checkbutton(self, text="Does the file have technical replicates?", variable=self.replicate_var).grid(
            row=2, column=0)

        self.experiments_list = tk.Listbox(self, height=3)
        self.experiments_list.grid(row=2, column=1)

        go_proteins_label = tk.Label(self, text="Go analysis proteins").grid(row=3, column=0)

        experiments_label = tk.Label(self, text="Pathway analysis").grid(row=3, column=1)

        self.go_proteins_list = tk.Listbox(self, selectmode="multiple", height=5, width=len(max(self.mqinit.possible_gos, key=len)))
        self.go_proteins_list.configure(exportselection=False)
        for x in MQInitializer.possible_gos:
            self.go_proteins_list.insert("end", x)

        self.go_proteins_list.grid(row=4, column=0)

        self.pathway_list = tk.Listbox(self, selectmode="multiple", height=5, width=len(max(MQInitializer.possible_pathways, key=len)))
        self.pathway_list.configure(exportselection=False)
        for x in MQInitializer.possible_pathways:
            self.pathway_list.insert("end", x)

        self.pathway_list.grid(row=4, column=1)

        plot_label = tk.Label(self, text="Which plots should be created").grid(row=5, column=0)

        intensity_label = tk.Label(self, text="Intensity").grid(row=5, column=1)

        levels_label = tk.Label(self, text="Levels").grid(row=5, column=2)

        self.heading_length = 6

        self.detection_counts_int, self.detection_counts_var, self.detection_counts_levels = self.plot_row("Detection counts", "raw_log2")
        self.number_of_detected_proteins_int, self.number_of_detected_proteins_var, self.number_of_detected_proteins_levels = self.plot_row("Number of detected proteins", "raw_log2")
        self.intensity_histograms_int, self.intensity_histograms_var, self.intensity_histograms_levels = self.plot_row("Intensity histogram", "raw")
        self.relative_std_int, self.relative_std_var, self.relative_std_levels = self.plot_row("Relative std", "raw_log2")
        self.rank_int, self.rank_var, self.rank_levels = self.plot_row("Rank", "raw_log2")
        self.pathway_analysis_int, self.pathway_analysis_var, self.pathway_analysis_levels = self.plot_row("Pathway Analysis", "raw_log2")
        self.pathway_timeline_int, self.pathway_timeline_var, self.pathway_timeline_levels = self.plot_row("Pathway Timeline", "raw_log2")
        # self.pathway_proportions_int, self.pathway_proportions_var = self.plot_row("Pathway proportions", "raw_log2")
        self.scatter_replicates_int, self.scatter_replicates_var, self.scatter_replicates_levels = self.plot_row("Scatter replicates", "raw_log2")
        self.experiment_comparison_int, self.experiment_comparison_var, self.experiment_comparison_levels = self.plot_row("Experiment comparison", "raw_log2")
        self.go_analysis_int, self.go_analysis_var, self.go_analysis_levels = self.plot_row("Go analysis", "raw_log2")
        self.venn_results_int, self.venn_results_var, self.venn_results_levels = self.plot_row("Venn diagrams", "raw_log2")
        self.venn_groups_int, self.venn_groups_var, self.venn_groups_levels = self.plot_row("Group diagrams", "raw_log2")
        self.r_volcano_int, self.r_volcano_var, self.r_volcano_levels = self.plot_row("Volcano plot (R)", "raw_log2")

        total_length = self.heading_length + self.number_of_plots
        update_button = tk.Button(self, text="Update", command=lambda: self.update_button())
        update_button.grid(row=total_length + 1, column=0)

        start_button = tk.Button(self, text="Start",
                                 command=lambda: self.start_button())
        start_button.grid(row=total_length + 1, column=1)

        self.running_text = tk.StringVar(value="Please press Start")
        self.running_label = tk.Label(self, textvariable=self.running_text).grid(row=total_length + 2, column=1)

        # add all tracing to the variables
        self.dir_text.trace("w", self.dir_setter)
        self.yaml_text.trace("w", self.yaml_path_setter)

    def plot_row(self, text: str, intensity_default: str):
        row = self.heading_length + self.number_of_plots
        int_var = tk.IntVar(value=1)
        plot_button = tk.Checkbutton(self, text=text, variable=int_var).grid(row=row, column=0)

        intensity_var = tk.StringVar(value=intensity_default)

        intensity_menu = tk.OptionMenu(self, intensity_var, "lfq", "raw", "ibaq", "lfq_log2", "raw_log2",
                                       "ibaq_log2").grid(row=row, column=1)
        level_list = tk.Listbox(self, selectmode="multiple", height=2)
        level_list.grid(row=row, column=2)
        level_list.configure(exportselection=False)
        for op in [0, 1, 2]:
            level_list.insert("end", op)
        self.number_of_plots += 1
        return int_var, intensity_var, level_list


class MQParser(argparse.ArgumentParser):
    def __init__(self):
        super().__init__(description="A pipeline to analyze result files from a MaxQuant report. "
                                     "The required result files are in the txt directory.")
        self.add_argument(
            '--dir',
            dest='file_dir',
            action='store',
            default="",
            help="Path to directory of analysis which should be analyzed."
                 "If not set the program will open a prompt and ask to select one."
        )
        self.add_argument(
            '--yml-file',
            dest='yml_file',
            action='store',
            default="file",
            help="Path to yml file which should be used for analysis, or 'default' / 'file'."
        )
        self.add_argument(
            "--loglevel",
            dest="loglevel",
            action="store",
            default=logging.WARNING,
            help="Logging level of analysis. Should be from options (lowest to highest): DEBUG < INFO < WARNING < ERROR. "
                 "The higher the logging level the fewer messages are shown. Default: WARNING"
        )
        self.add_argument(
            "--has-replicates",
            dest="has_replicates",
            default=False,
            const=True,
            nargs="?",
            help="If you have replicates of each experiment specify this"
                 "Replicates need to be enumerated in a way that numbers are added at the end of the name"
                 "If no replicates are in data set no venn diagrams will be generated"
        )
        self.add_argument(
            "--has-groups",
            dest="has_groups",
            default=False,
            const=True,
            nargs="?",
            help="If you have groups of experiments specify this"
                 "Replicates need to be enumerated in a way that numbers are added at the end of the name"
                 "If no replicates are in data set no venn diagrams will be generated"
        )
        self.add_argument(
            "--gui",
            default=False,
            const=True,
            nargs="?",
            help="specify this if a gui should be opened"
        )
        self.args = self.parse_args()
        self.args_dict = vars(self.args)
        move_to_config = ["has_replicates", "has_groups"]
        self.args_dict["configs"] = {k: self.args_dict.pop(k) for k in move_to_config}


def browsefunc(fn, var, fn_params: dict = None):
    if fn_params is None:
        fn_params = {}
    filename = fn(**fn_params)
    var.set(filename)
