from mspypeline import MSPParser, MSPUI

if __name__ == "__main__":
    mspparser = MSPParser()
    gui = MSPUI(**mspparser.args_dict)
