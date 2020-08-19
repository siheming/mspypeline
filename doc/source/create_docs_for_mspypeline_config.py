import os


def main():
    go_dir = "../../mspypeline/config/go_terms"
    pathway_dir = "../../mspypeline/config/pathways"

    go_rst_dir = os.path.abspath("./go_terms")
    pathway_rst_dir = os.path.abspath("./pathways")

    def add_dir(source_dir, goal_dir):
        os.makedirs(goal_dir, exist_ok=True)
        for root, dirs, files in os.walk(os.path.abspath(source_dir)):
            for file in files:
                new_name = file.replace(".txt", "").replace("_", " ").lower().title()
                with open(os.path.join(goal_dir, new_name + ".rst"), "w") as out:
                    out.write(new_name + "\n")
                    out.write("=" * len(new_name) + "\n")
                    out.write("\n")
                    out.write(f".. include:: ../{source_dir}/{file}" + "\n")
                    out.write("   :literal:" + "\n")
                    out.write("\n")

    add_dir(go_dir, go_rst_dir)
    add_dir(pathway_dir, pathway_rst_dir)


if __name__ == "__main__":
    main()
