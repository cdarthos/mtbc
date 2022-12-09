import flask
from mtbc import MtbcRandom

app = flask.Flask(__name__)


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/mtbc")
def mtbc():
    tree = mtbc.MtbcRandom(debug=True, list_length=30)
    return tree.nj_tree

if __name__ == "__main__":
    test = MtbcRandom(debug=True, list_length=30)
    # print(mtbc.construct_search_request())
    # print()
    # print(mtbc.id_list[:10])
    # print()
    # print(mtbc.acc_list)
    # print()
    # print(json.dumps(mtbc.alignement, indent=4))

