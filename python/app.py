import flask
from mtbc import MtbcRandom

app = flask.Flask(__name__)


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/mtbc")
def flask_mtbc():
    tree = MtbcRandom()
    return tree.nj_tree

if __name__ == "__main__":
    mtbc = MtbcRandom(debug=True, list_length=30)