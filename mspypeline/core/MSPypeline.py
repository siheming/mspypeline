import argparse
import os
import tkinter as tk
from tkinter import filedialog
import logging
from typing import Optional, Iterable
import webbrowser

from mspypeline import create_app
from mspypeline.core import MSPInitializer
from mspypeline.helpers import get_logger
from mspypeline.modules import default_normalizers
from mspypeline.file_reader import BaseReader, MQReader


class UIHandler:
    """
    | Used to take the mapping of arguments provided by the :class:`~MSPParser` through the command line to organize
      through which point of entry the data analysis should be performed.
    | The UIHandler can initialize the creation of the GUI, perform the whole analysis instantly by running the
      pipeline and creating all plots according to the configs or (not yet available) host the mspypeline on a flask
      server.
    """
    def __init__(self, file_dir, yml_file=None, gui=False, host_flask=False, selected_reader=MQReader.MQReader,
                 loglevel=logging.DEBUG, configs: dict = None):
        """
        Parameters
        ----------
        file_dir
            location where the directory/txt folder to the data can be found.
        yml_file
            path to the yaml config file
        gui
            should a GUI be compiled
        host_flask
            currently not implemented (should a flsk server be hosted)
        selected_reader
            instance of an :class:`~BaseReader` used to process data to an internal format for the plotter.
        loglevel
            level of the logger
        configs
            mapping of configuration
        """
        base_config = {
            "has_techrep": False,
            "has_groups": False,
        }
        if configs is None:
            configs = {}
        base_config.update(**configs)
        configs = base_config

        if gui and host_flask:
            raise ValueError("Can only specify one of host_flask and gui")
        if gui:
            MSPGUI(file_dir=file_dir, yml_file=yml_file, loglevel=loglevel, configs=configs)
        elif host_flask:
            # TODO pass arguments to create app
            app = create_app()
            app.run(debug=True)
        else:
            mspinit = MSPInitializer(file_dir, yml_file, loglevel=loglevel)
            mspinit.init_config()
            mspinit.configs.update(configs)
            mspinit.read_data()
            # create plotter from initializer
            mspplots = selected_reader.plotter.from_MSPInitializer(mspinit)
            # create all plots and other results
            mspplots.create_results()


class MSPGUI(tk.Tk):
    """
    | This class is called from the :class:`~UIHandler` to create the GUI.
    """
    def __init__(self, file_dir, yml_file=None, loglevel=logging.DEBUG, configs: dict = None):
        """
        Parameters
        ----------
        file_dir
            location where the directory/txt folder to the data can be found.
        yml_file
            path to the yaml config file
        loglevel
            level of the logger
        configs
            mapping containing the configurations
        """
        super().__init__()
        self.yaml_options = ["default"]
        self.reader_options = {reader.name: reader for reader in BaseReader.__subclasses__()}
        self.selected_reader = MQReader.MQReader
        self.normalize_options = ["None"] + list(default_normalizers.keys())
        self.mspinit = MSPInitializer(file_dir, yml_file, loglevel=loglevel)
        self.logger = get_logger(self.__class__.__name__, loglevel=loglevel)

        self.number_of_plots = 0

        self.plot_settings = {}
        self.intensity_options = ["lfq_log2", "raw_log2", "ibaq_log2"]
        #,"lfq_normalized_log2", "raw_normalized_log2", "ibaq_normalized_log2]

        self.title("mspypeline")

        path_label = tk.Label(self, text="Dir to analyze", font="Helvetica 10 bold").grid(row=0, column=0)

        yaml_label = tk.Label(self, text="Yaml file", font="Helvetica 10 bold").grid(row=0, column=1)

        reader_label = tk.Label(self, text="File reader", font="Helvetica 10 bold").grid(row=0, column=2)

        self.dir_text = tk.StringVar(value=file_dir)
        dir_button = tk.Button(self, textvariable=self.dir_text,
                               command=lambda: browsefunc(filedialog.askdirectory, self.dir_text, fn_params={
                                   "title": "Please select a directory with MaxQuant result files"}))
        create_tool_tip(dir_button, "Select a directory to analyze")
        dir_button.grid(row=1, column=0)

        self.yaml_text = tk.StringVar()
        self.yaml_button = tk.OptionMenu(self, self.yaml_text, *self.yaml_options)
        self.yaml_button.grid(row=1, column=1)

        self.reader_text = tk.StringVar(value="mqreader")
        self.reader_button = tk.OptionMenu(self, self.reader_text, *self.reader_options.keys())
        self.reader_button.grid(row=1, column=2)

        self.replicate_var = tk.IntVar(value=1)
        replicate_button = tk.Checkbutton(self, text="Does the file have technical replicates?",
                                          variable=self.replicate_var).grid(row=2, column=0)

        go_proteins_label = tk.Label(self, text="Go analysis lists")
        create_tool_tip(go_proteins_label, "For a full list of Proteins see the Documentation."
                                           "\nCustom Gene Lists -> GO Terms")
        go_proteins_label.grid(row=3, column=0)

        experiments_label = tk.Label(self, text="Pathway analysis lists")
        create_tool_tip(experiments_label, "For a full list of Proteins see the Documentation."
                                           "\nCustom Gene Lists -> Pathways")
        experiments_label.grid(row=3, column=1)

        design_label = tk.Label(self, text="Sample names")
        create_tool_tip(design_label, "Inferred sample names of the experiment")
        design_label.grid(row=3, column=2)

        self.go_term_list = tk.Listbox(self, selectmode="multiple", height=5,
                                       width=len(max(self.mspinit.list_full_gos, key=len)))
        self.go_term_list.configure(exportselection=False)
        for x in self.mspinit.list_full_gos:
            self.go_term_list.insert("end", x)

        self.go_term_list.grid(row=4, column=0)

        self.pathway_list = tk.Listbox(self, selectmode="multiple", height=5,
                                       width=len(max(self.mspinit.list_full_pathways, key=len)))
        self.pathway_list.configure(exportselection=False)
        for x in self.mspinit.list_full_pathways:
            self.pathway_list.insert("end", x)

        self.pathway_list.grid(row=4, column=1)

        self.experiments_list = tk.Listbox(self, height=5)
        self.experiments_list.grid(row=4, column=2)

        report_button = tk.Button(self, text="Create Report",
                                  command=lambda: self.report_button())
        report_button.grid(row=5, column=0)
        create_tool_tip(report_button, "MaxQuant report for quality control")

        plot_label = tk.Label(self, text="Which plots should be created").grid(row=6, column=0)

        intensity_label = tk.Label(self, text="Intensities").grid(row=6, column=1)

        levels_label = tk.Label(self, text="Levels").grid(row=6, column=2)

        self.heading_length = 7

        tk.Label(self, text="Normalization plots", font="Helvetica 10 bold").grid(
            row=self.heading_length + self.number_of_plots, column=0)
        self.number_of_plots += 1
        self.plot_row("Normalization overview", "normalization_overview_all_normalizers",
                      "Displays original data and data after different normalization methods")
        self.plot_row("Heatmap overview", "heatmap_overview_all_normalizers",
                      "Displays intensity heatmap of original data and after different normalization methods")

        norm_method_label = tk.Label(self, text="Choose a Normalization Method:", font="Helvetica 10 bold")
        norm_method_label.grid(row=self.heading_length + self.number_of_plots, column=1)
        create_tool_tip(norm_method_label, "For more information about normalization visit the documentation.\n"
                        "Can be found in Workflow -> Data Preprocessing -> Normalization options.")
        self.number_of_plots += 1
        self.plot_intermediate_row("Choose a Normalization Method")

        tk.Label(self, text="Outlier detection / Comparisons", font="Helvetica 10 bold").grid(
            row=self.heading_length + self.number_of_plots, column=0)
        self.number_of_plots += 1
        self.plot_row("Detection counts", "detection_counts",
                      "How many proteins were detected how frequently in the samples of a group?")
        self.plot_row("Number of detected proteins", "detected_proteins_per_replicate",
                      "How many proteins were detected in each of my samples and in total for each group?")
        self.plot_row("Venn diagrams", "venn_results",
                      "How large is the intersection of detected proteins of my samples in each group? How many proteins are uniquely detected in a sample?")
        self.plot_row("Group diagrams", "venn_groups",
                      "How large is the intersection of detected proteins between different groups? How many proteins are uniquely detected in a group?")
        self.plot_row("PCA overview", "pca_overview",
                      "How similar are my samples? Do samples cluster together?")
        self.plot_row("Intensity histogram", "intensity_histograms",
                      "How does the intensity profile of my samples look? How similar are the intensity profiles?")
        self.plot_row("Relative std", "relative_std",
                      "What is the relative standard deviation of the samples of a group?")
        self.plot_row("Scatter replicates", "scatter_replicates",
                      "How well do the overall protein intensities of the samples of each group correlate?")
        self.plot_row("Experiment comparison", "experiment_comparison",
                      "How well do the overall protein intensities of different groups correlate?")
        self.plot_row("Rank", "rank",
                      "Where do my proteins of interest rank in intensity compared to all other proteins?")

        tk.Label(self, text="Statistical inference", font="Helvetica 10 bold").grid(
            row=self.heading_length + self.number_of_plots, column=0)
        self.number_of_plots += 1
        self.plot_row("Pathway Analysis", "pathway_analysis",
                      "What is the intensity of my proteins of interest, and is it significantly different in one group versus the other?")
        #self.plot_row("Pathway Timecourse", "pathway_timecourse")
        self.plot_row("Go analysis", "go_analysis",
                      "Are the proteins of a group enriched for the selected GO terms?")
        self.plot_row("Volcano plot (R)", "r_volcano",
                      "Which proteins are significantly higher or lower in intensity comparing two groups? Which proteins are detected only in one group and not in the other?")
        self.p_val_var = tk.IntVar(value=1)
        pval_button = tk.Checkbutton(self, text="Use adjusted p value", variable=self.p_val_var).grid(
            row=self.heading_length + self.number_of_plots, column=1)

        total_length = self.heading_length + self.number_of_plots

        update_button = tk.Button(self, text="Update", command=lambda: self.update_button())
        update_button.grid(row=total_length + 1, column=1)

        start_button = tk.Button(self, text="Start",
                                 command=lambda: self.start_button())
        start_button.grid(row=total_length + 1, column=2)

        documentation_link = tk.Label(self, text="Documentation link", font="Helvetica 10 underline",
                                      fg="blue", cursor="hand2")
        documentation_link.bind("<ButtonRelease-1>",
                                lambda _: webbrowser.open_new("https://mspypeline.readthedocs.io/en/latest/"))
        documentation_link.grid(row=total_length + 2, column=0)

        self.running_text = tk.StringVar(value="Please press Start")
        self.running_label = tk.Label(self, textvariable=self.running_text).grid(row=total_length + 2, column=2)

        # add all tracing to the variables
        self.dir_text.trace("w", self.dir_setter)
        self.yaml_text.trace("w", self.yaml_path_setter)
        self.reader_text.trace("w", self.reader_setter)
        # make the GUI resizable
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.mspinit.configs.update(configs)
        self.update_yaml_options()
        self.mainloop()

    def dir_setter(self, *args):
        start_dir = self.dir_text.get()
        if os.path.split(start_dir)[1] == "txt":
            start_dir = os.path.split(start_dir)[0]
        self.mspinit.start_dir = start_dir
        self.update_yaml_options()

    def update_yaml_options(self):
        if self.mspinit.has_yml_file():
            self.yaml_text.set("file")
            self.yaml_options = ["default", "file"]
            # next steps should also be done hy what the update button does?
        else:
            self.yaml_text.set("default")
            self.yaml_options = ["default"]
        self.yaml_button["menu"].delete(0, "end")
        for op in self.yaml_options:
            self.yaml_button["menu"].add_command(label=op, command=tk._setit(self.yaml_text, op))

    def update_listboxes(self):
        # delete all experiments then add from file
        self.experiments_list.delete(0, "end")
        for op in self.mspinit.configs.get(self.selected_reader.name, {}).get("all_replicates", []):
            self.experiments_list.insert("end", op)
        # clear selection then select from configs
        for i, pathway in enumerate(self.mspinit.list_full_pathways):
            self.pathway_list.select_clear(i)
        if self.mspinit.configs.get("pathways"):
            for pathway in self.mspinit.configs.get("pathways"):
                try:
                    self.pathway_list.select_set(self.mspinit.list_full_pathways.index(pathway))
                except ValueError:
                    self.logger.warning("Selected pathway file %s not found", pathway)
        # clear selection then select from configs
        for i, go in enumerate(self.mspinit.list_full_gos):
            self.go_term_list.select_clear(i)
        if self.mspinit.configs.get("go_terms"):
            for go in self.mspinit.configs.get("go_terms"):
                try:
                    self.go_term_list.select_set(self.mspinit.list_full_gos.index(go))
                except ValueError:
                    self.logger.warning("Selected go term file %s not found", go)

    def yaml_path_setter(self, *args):
        self.mspinit.file_path_yaml = self.yaml_text.get()
        # get the reader class by the saved name
        self.selected_reader = self.reader_options.get(self.mspinit.configs.get("selected_reader", "mqreader"))
        reader_settings = self.mspinit.configs.get(self.selected_reader.name, {})
        level_names = reader_settings.get("level_names", [])
        level_names = {i: name for i, name in enumerate(level_names)}
        levels = reader_settings.get("levels", 0)
        for plot_name in self.selected_reader.plotter.possible_plots:
            plot_settings_name = plot_name + "_settings"
            plot_settings = self.mspinit.configs.get(plot_settings_name, {})
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("var", "int")
            levels_name = var_name.replace("var", "levels")
            # update all settings in the GUI
            self.plot_settings[int_name].set(plot_settings.get("create_plot", False))

            plot_intensities = plot_settings.get("dfs_to_use", [])
            self.plot_settings[var_name].update_selection(plot_intensities)

            plot_levels = plot_settings.get("levels", [])
            selected_levels = self.plot_settings[levels_name]
            selected_levels.update_options([level_names.get(l, l) for l in range(levels)])
            selected_levels.update_selection([level_names.get(pl, pl) for pl in plot_levels])
        self.replicate_var.set(self.mspinit.configs.get("has_techrep", True))
        self.p_val_var.set(self.mspinit.configs.get("plot_r_volcano_settings", {}).get("adj_pval", False))
        self.normalizer_text.set(self.mspinit.configs.get("selected_normalizer", "None"))
        self.update_listboxes()

    def reader_setter(self, *args):
        self.selected_reader = self.reader_options[self.reader_text.get()]

    def update_button(self):
        self.mspinit.configs["has_techrep"] = bool(self.replicate_var.get())
        self.mspinit.configs["selected_reader"] = str(self.reader_text.get())
        self.mspinit.configs["selected_normalizer"] = str(self.normalizer_text.get())
        reader_settings = self.mspinit.configs.get(self.selected_reader.name, {})
        level_names = reader_settings.get("level_names", [])
        level_names = {name: i for i, name in enumerate(level_names)}
        for plot_name in self.selected_reader.plotter.possible_plots:
            plot_settings = plot_name + "_settings"
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("var", "int")
            levels_name = var_name.replace("var", "levels")
            selected_settings = {
                "create_plot": bool(self.plot_settings[int_name].get()),
                "dfs_to_use": [k for k, v in self.plot_settings[var_name].get_selection().items() if v],
                "levels": [level_names.get(k, k) for k, v in self.plot_settings[levels_name].get_selection().items() if v]
            }
            additional_settings = self.mspinit.configs.get(plot_settings, {})
            additional_settings = {k: v for k, v in additional_settings.items()
                                   if k != "create_plot" and k != "dfs_to_use" and k != "levels"}
            selected_settings.update(additional_settings)
            for k, v in selected_settings.items():
                self.mspinit.configs[plot_settings][k] = v
        gos = self.go_term_list.curselection()
        gos = [self.mspinit.list_full_gos[int(go)] for go in gos]
        pathways = self.pathway_list.curselection()
        pathways = [self.mspinit.list_full_pathways[int(pathway)] for pathway in pathways]
        self.mspinit.configs["plot_r_volcano_settings"]["adj_pval"] = bool(self.p_val_var.get())
        self.mspinit.configs["go_terms"] = gos
        self.mspinit.configs["pathways"] = pathways
        self.mspinit.init_config()
        self.mspinit.read_data()
        self.update_listboxes()
        self.update_yaml_options()

    def start_button(self):
        self.running_text.set("Creating Plots")
        self.update()
        self.update_button()
        mspplots = self.selected_reader.plotter.from_MSPInitializer(self.mspinit)
        mspplots.create_results()
        self.running_text.set("Please press Start")

    def report_button(self):
        self.running_text.set("Creating Report")
        self.update()
        self.update_button()
        mspplots = self.selected_reader.plotter.from_MSPInitializer(self.mspinit)
        mspplots.create_report()
        self.running_text.set("Please press Start")

    def plot_intermediate_row(self, text: str):
        row = self.heading_length + self.number_of_plots
        int_var = tk.IntVar(value=1)
        self.normalizer_text = tk.StringVar(value="None")
        self.normalizer_button = tk.OptionMenu(self, self.normalizer_text, *self.normalize_options)
        self.normalizer_button.grid(row=row, column=1)

        self.number_of_plots += 1

    def plot_row(self, text: str, plot_name: str, plot_tool_tip: str = None):
        row = self.heading_length + self.number_of_plots
        int_var = tk.IntVar(value=1)
        checkbutton = tk.Checkbutton(self, text=text, variable=int_var)
        checkbutton.grid(row=row, column=0)
        if plot_tool_tip:
            create_tool_tip(checkbutton, plot_tool_tip)

        intensity_list = MultiSelectOptionMenu(self, self.intensity_options, "Select Intensities")
        intensity_list.grid(row=row, column=1)

        level_list = MultiSelectOptionMenu(self, button_text="Select Levels")
        level_list.grid(row=row, column=2)

        self.number_of_plots += 1
        self.plot_settings.update({
            f"{plot_name}_int": int_var,
            f"{plot_name}_var": intensity_list,
            f"{plot_name}_levels": level_list
        })


class MultiSelectOptionMenu(tk.Frame):
    def __init__(self, parent, choices: Optional[Iterable] = None, button_text: str = "Default text"):
        super().__init__(parent)
        menubutton = tk.Menubutton(self, text=button_text, indicatoron=True, borderwidth=1, relief="raised")
        self.menu = tk.Menu(menubutton, tearoff=False)
        menubutton.configure(menu=self.menu)
        menubutton.pack(pady=3, padx=3)
        self.choices_dict = {}
        self.choices = choices if choices is not None else ()
        self.update_options()

    def update_options(self, choices: Optional[Iterable] = None):
        if choices is not None:
            self.choices = choices
            self.choices_dict.clear()
            self.menu.delete(0, "end")
        for choice in self.choices:
            self.choices_dict[choice] = tk.BooleanVar(value=False)
            self.menu.add_checkbutton(label=choice, variable=self.choices_dict[choice], onvalue=True, offvalue=False)

    def update_selection(self, choices: Iterable):
        for choice in self.choices:
            if choice in choices:
                self.choices_dict[choice].set(True)
            else:
                self.choices_dict[choice].set(False)

    def get_selection(self):
        return {k: v.get() for k, v in self.choices_dict.items()}


class ToolTip:
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None

    def showtip(self, text: str):
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() + 27
        self.tipwindow = tk.Toplevel(self.widget)
        self.tipwindow.wm_overrideredirect(1)
        self.tipwindow.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tipwindow, text=text, justify=tk.LEFT,
                         background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        if self.tipwindow:
            self.tipwindow.destroy()


def create_tool_tip(widget, text):
    tool_tip = ToolTip(widget)
    widget.bind('<Enter>', lambda x: tool_tip.showtip(text))
    widget.bind('<Leave>', lambda x: tool_tip.hidetip())


class MSPParser(argparse.ArgumentParser):
    """
    | Uses the ``argparse`` module to provide a parser for command line options, arguments and sub-commands through
      which the analysis can be started or the GUI can be called (see :ref:`get-started`).
    """
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
            default=None,
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
            "--has-techrep",
            dest="has_techrep",
            default=False,
            const=True,
            nargs="?",
            help="If you have replicates of each experiment specify this"
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
        self.add_argument(
            "--host-flask",
            dest="host_flask",
            default=False,
            const=True,
            nargs="?",
            help="specify this if a flask server should be hosted"
        )
        self.args = self.parse_args()
        self.args_dict = vars(self.args)
        move_to_config = ["has_techrep"]
        self.args_dict["configs"] = {k: self.args_dict.pop(k) for k in move_to_config}


def browsefunc(fn, var, fn_params: dict = None):
    if fn_params is None:
        fn_params = {}
    filename = fn(**fn_params)
    var.set(filename)
