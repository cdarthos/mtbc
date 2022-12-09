from flask import Flask
from mtbc import MtbcRandom

app = Flask(__name__)


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/mtbc")
def mtbc():
    mtbc = MtbcRandom(debug=True, list_length=30)
    return mtbc.nj_tree


