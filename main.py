from mqpipeline import MQParser, MQUI

if __name__ == "__main__":
    mqparser = MQParser()
    gui = MQUI(**mqparser.args_dict)
