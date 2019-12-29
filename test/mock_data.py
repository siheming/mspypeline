import pandas as pd
import numpy as np
import os
import shutil


class MockData:
    test_dir = os.path.dirname(os.path.realpath(__file__))
    script_loc = os.path.split(test_dir)[0]
    config_dir = os.path.join(script_loc, "mqpipeline", "config")
    go_dir = os.path.join(config_dir, "go_terms")
    pathway_dir = os.path.join(config_dir, "pathways")
    mock_data_dir = os.path.join(test_dir, "mock_data")

    @staticmethod
    def create_mock_data(
            test_cases = (0, 1, 10, 20, 50, 100, 250, 500),
            number_of_non_pathway_genes = 1000,
            mu = 28,
            sigma = 2,
            tech_rep_sigma = 0.1,
            experiment_sigma = 1,
            group_sigma = 2,
            design_combinations=((True, True), (False, True), (True, False), (False, False)),  #(has group, has technical replicates)
            number_of_technical_replicates = (3, 8, 1, 1),  # bad fix to make it easy
            number_of_experiments = (3, 6, 15, 30),
            number_of_groups = (4, 1, 8, 1),  # bad fix to make it easy
            seed = 100
    ):
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

        # this will be exp2 transformed later to create log normal distributed data
        gene_base_abundance = np.random.normal(mu, sigma, (N,))

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

            group_noise = np.random.normal(0, group_sigma, (N, n_groups))
            experiment_noise = np.random.normal(0, experiment_sigma, (N, n_groups, n_experiments))
            tech_rep_noise = np.random.normal(0, tech_rep_sigma, (N, n_groups, n_experiments, n_tech_reps))

            """group_noise = np.zeros((N, n_groups))
            # then add some constant per dimension
            for j in range(1, n_groups + 1):
                b = 10000
                group_noise[:, j - 1] = b * j
            experiment_noise = np.zeros((N, n_groups, n_experiments))
            # then add some constant per dimension
            for j in range(1, n_experiments + 1):
                b = 1000
                experiment_noise[:, :, j - 1] = b * j


            tech_rep_noise = np.zeros((N, n_groups, n_experiments, n_tech_reps))
            # then add some constant per dimension
            for j in range(1, n_tech_reps + 1):
                b = 100
                tech_rep_noise[:, :, :, j - 1] = b * j"""
            # add the noises togehter, each time add a new axis
            experiment_data = np.broadcast_to(gene_base_abundance[..., np.newaxis],
                                              (N, n_groups)) + group_noise
            experiment_data = np.broadcast_to(experiment_data[..., np.newaxis],
                                              (N, n_groups, n_experiments)) + experiment_noise
            experiment_data = np.broadcast_to(experiment_data[..., np.newaxis],
                                              (N, n_groups, n_experiments, n_tech_reps)) + tech_rep_noise

            assert experiment_data.shape == (N, n_groups, n_experiments, n_tech_reps)

            ex = experiment_data.reshape((N, n_groups * n_experiments * n_tech_reps), order="C")
            df = pd.DataFrame(ex, index=genes, columns=[f"Intensity group{g + 1}_experiment{e + 1}_rep{t + 1}"
                                                        for g in range(n_groups) for e in range(n_experiments) for t in
                                                        range(n_tech_reps)])
            df = np.exp2(df)
            df_lfq = df.rename({col: col.replace("Intensity", "LFQ") for col in df}, axis=1)
            df_ibaq = df.rename({col: col.replace("Intensity", "iBAQ") for col in df}, axis=1)
            df = pd.concat([df, df_lfq, df_ibaq], axis=1)
            # TODO drop some unimportant genes randomly

            # df["Base"] = gene_base_abundance
            # required csv columns: "Fasta headers", "Only identified by site", "Reverse", "Potential contaminant", all empty
            # required col "Gene names", which is to be sampled from the pathways
            df["Fasta headers"] = ""
            df["Only identified by site"] = ""
            df["Reverse"] = ""
            df["Potential contaminant"] = ""
            df = df.reset_index()
            dir_name = f"{'has_group' if has_group else 'no_group'}_{'has_tech' if has_tech_rep else 'no_tech'}"
            dir_name = os.path.join(MockData.mock_data_dir, dir_name, "txt")
            os.makedirs(dir_name, exist_ok=True)
            df.to_csv(os.path.join(dir_name, "proteinGroups.txt"), index=False, header=True)

    @staticmethod
    def delete_mock_data():
        shutil.rmtree(MockData.mock_data_dir)


if __name__ == "__main__":
    MockData.delete_mock_data()
    MockData.create_mock_data()
