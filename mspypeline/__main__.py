if __name__ == "__main__":
    from . import MSPParser, UIHandler
    mspparser = MSPParser()
    ui = UIHandler(**mspparser.args_dict)
