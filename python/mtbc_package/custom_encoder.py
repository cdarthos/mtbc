from json import JSONEncoder


class customEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__
