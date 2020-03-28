if __name__ == "__main__":
    from mqpipeline import MQParser, UIHandler
    mqparser = MQParser()
    UIHandler(**mqparser.args_dict)
