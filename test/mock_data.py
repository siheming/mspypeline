import pandas as pd
import numpy as np
import os
import shutil


class MockData:
    test_dir = os.path.dirname(os.path.realpath(__file__))
    script_loc = os.path.split(test_dir)[0]
    config_dir = os.path.join(script_loc, "mspypeline", "config")
    go_dir = os.path.join(config_dir, "go_terms")
    pathway_dir = os.path.join(config_dir, "pathways")
    mock_data_dir = os.path.join(test_dir, "mock_data")

    @staticmethod
    def create_mock_data(
            test_cases = (0, 1, 10, 20, 50, 100, 250, 500),
            number_of_non_pathway_genes = 1000,
            average = 28,
            gene_sigma = 2,
            group_sigma = 1,
            experiment_sigma = 1,
            tech_rep_sigma = 0.3,
            noise_sigma = 0.5,
            design_combinations=((True, True), (False, True), (True, False), (False, False)),  #(has group, has technical replicates)
            number_of_technical_replicates = (3, 8, 1, 1),  # bad fix to make it easy
            number_of_experiments = (3, 6, 15, 30),
            number_of_groups = (4, 1, 8, 1),  # bad fix to make it easy
            seed = 100,
            save_to_disk=True
    ):
        # TODO this seed seems insufficient?
        np.seed = seed
        N = sum(test_cases) + number_of_non_pathway_genes

        os.makedirs(MockData.mock_data_dir, exist_ok=True)

        # accumulate all genes
        pathway_genes = {}
        for file_name in os.listdir(MockData.pathway_dir):
            file = os.path.join(MockData.pathway_dir, file_name)
            with open(file) as f:
                pathway = f.readline().strip()
                pathway_genes[pathway] = []
                f.readline()
                for line in f:
                    pathway_genes[pathway].append(line.strip())

        # sort the pathways into the different possible cases
        # TODO this samples a lot of duplicate gene names
        free_pathway_genes = set(pathway_genes)
        test_case_dict = {t: [] for t in test_cases}
        for test_case in reversed(sorted(test_cases)):
            to_rm = set()
            for pathway in free_pathway_genes:
                if len(pathway_genes[pathway]) >= test_case:
                    test_case_dict[test_case].append(pathway)
                    to_rm.add(pathway)
            free_pathway_genes = free_pathway_genes - to_rm

        # randomly sample the genes from the pathways
        genes = []
        for test_case in test_cases:
            if test_case == 0:
                continue
            pathway = np.random.choice(test_case_dict[test_case], 1)[0]
            genes += list(np.random.choice(pathway_genes[pathway], test_case, replace=False))

        # add the non pathway genes
        for i in range(number_of_non_pathway_genes):
            genes.append(f"GENENP{i}")

        assert len(genes) == N
        genes = pd.Series(genes, name="Gene names")

        for i, (has_group, has_tech_rep) in enumerate(design_combinations):
            n_experiments = number_of_experiments[i]
            n_tech_reps = number_of_technical_replicates[i]
            n_groups = number_of_groups[i]
            print(f"n group: {n_groups}, n exp: {n_experiments}, n tech: {n_tech_reps}")
            if has_group:
                assert n_groups > 1
            else:
                assert n_groups == 1
            if has_tech_rep:
                assert n_tech_reps > 1
            else:
                assert n_tech_reps == 1
            assert n_experiments > 1

            group_name = np.array([f"Group{x}_" for x in range(n_groups)]).reshape((n_groups, 1, 1))
            experiment_name = np.array([f"Experiment{x}_" for x in range(n_experiments)]).reshape((1, n_experiments, 1))
            technical_replicate_name = np.array([f"Rep{x}" for x in range(n_tech_reps)]).reshape((1, 1, n_tech_reps))

            names = np.core.defchararray.add(
                np.core.defchararray.add(group_name, experiment_name), technical_replicate_name
            ).reshape((n_groups * n_experiments * n_tech_reps))

            average_effect = np.ones((N, n_groups, n_experiments, n_tech_reps)) * average
            gene_effect = np.random.normal(0, gene_sigma, (N, 1, 1, 1))
            group_effect = np.random.uniform(0, group_sigma, (N, n_groups, 1, 1))
            experiment_effect = np.random.normal(0, experiment_sigma, (1, 1, n_experiments, 1))
            technical_replicate_effect = np.random.normal(0, tech_rep_sigma, (1, 1, 1, n_tech_reps))
            noise = np.random.normal(0, noise_sigma, (N, n_groups, n_experiments, n_tech_reps))
            experiment_data = average_effect + gene_effect + group_effect + experiment_effect + technical_replicate_effect + noise

            assert experiment_data.shape == (N, n_groups, n_experiments, n_tech_reps)

            ex = experiment_data.reshape((N, n_groups * n_experiments * n_tech_reps))
            df = pd.DataFrame(ex, index=genes, columns=names)
            df = df.rename(lambda x: f"Intensity {x}", axis=1)
            df = np.exp2(df)
            df_lfq = df.rename({col: col.replace("Intensity", "LFQ intensity") for col in df}, axis=1)
            df_ibaq = df.rename({col: col.replace("Intensity", "iBAQ") for col in df}, axis=1)
            df = pd.concat([df, df_lfq, df_ibaq], axis=1)
            # TODO drop some unimportant genes randomly

            # required csv columns: "Fasta headers", "Only identified by site", "Reverse", "Potential contaminant", all empty
            # required col "Gene names", which is to be sampled from the pathways
            df["Fasta headers"] = ""
            df["Only identified by site"] = ""
            df["Reverse"] = ""
            df["Potential contaminant"] = ""
            df = df.reset_index()
            df["Protein names"] = df["Gene names"] + "P"
            if save_to_disk:
                dir_name = f"{'has_group' if has_group else 'no_group'}_{'has_tech' if has_tech_rep else 'no_tech'}"
                dir_name = os.path.join(MockData.mock_data_dir, dir_name, "txt")
                os.makedirs(dir_name, exist_ok=True)
                df.to_csv(os.path.join(dir_name, "proteinGroups.txt"), index=False, header=True, sep="\t")
            else:
                return df

    @staticmethod
    def delete_mock_data():
        if os.path.exists(MockData.mock_data_dir):
            shutil.rmtree(MockData.mock_data_dir)


if __name__ == "__main__":
    MockData.delete_mock_data()
    MockData.create_mock_data()
