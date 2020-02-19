from mqpipeline.core import MQParser, UIHandler

if __name__ == "__main__":
    mqparser = MQParser()
    UIHandler(**mqparser.args_dict)
