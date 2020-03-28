if __name__ == "__main__":
    from mspypeline import MSPParser, UIHandler
    mspparser = MSPParser()
    gui = UIHandler(**mspparser.args_dict)
