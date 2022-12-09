import os

from flask import Flask, render_template, request
from mtbc_package import MtbcRandom

app = Flask(__name__)


@app.route("/mtbc", methods=['GET'])
def hello_world():
    return """<form action="/mtbc_package" method="POST">
  <div>
    <label for="debug">Debug?</label>
    <input name="debug" id="debug" value="True" type="checkbox" />
  </div>
  <div>
    <label for="size">list size</label>
    <input name="size" id="size" value="30" />
  </div>
  <div>
    <button>go</button>
  </div>
</form>
"""


@app.route("/mtbc", methods=['POST'])
def flask_mtbc():
    debug = False
    if request.form.get('debug'):
        debug = True

    print("debug")
    print(debug)
    print(type(debug))
    tree = MtbcRandom(debug=debug, list_length=int(request.form['size']))

    return os.listdir('tree')


if __name__ == "__main__":
    mtbc = MtbcRandom(debug=True, list_length=30)
