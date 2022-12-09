import os

from flask import Flask, render_template, request
from mtbc import MtbcRandom

app = Flask(__name__)


@app.route("/")
def hello_world():
    return """<form action="/mtbc" method="POST">
  <div>
    <label for="debug">Debug?</label>
    <input name="debug" id="debug" value="True" type="checkbox" />
  </div>
  <div>
    <label for="size">list size</label>
    <input name="size" id="size" value="30" />
  </div>
  <div>
    <button>Send my greetings</button>
  </div>
</form>
"""


@app.route("/mtbc", methods=['GET', 'POST'])
def flask_mtbc():
    debug = False
    if bool(request.form['debug']):
        debug = True
    tree = MtbcRandom(debug=debug, list_length=int(request.form['size']))

    return os.listdir('tree')


if __name__ == "__main__":
    mtbc = MtbcRandom(debug=True, list_length=30)
