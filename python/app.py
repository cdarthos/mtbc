import os

from flask import Flask, render_template, request
from mtbc import MtbcRandom

app = Flask(__name__)


@app.route("/")
def hello_world():
    return """<form action="/mtbc" method="POST">
  <div>
    <label for="debug">Degug?</label>
    <input name="debug" id="debug" value="True" />
  </div>
  <div>
    <label for="to">Who do you want to say it to?</label>
    <input name="to" id="to" value="Mom" />
  </div>
  <div>
    <button>Send my greetings</button>
  </div>
</form>
"""


@app.route("/mtbc", methods=['GET', 'POST'])
def flask_mtbc():
    tree = MtbcRandom(debug=request.form['debug'], list_length=int(request.form['to']))

    return os.listdir('tree')


if __name__ == "__main__":
    mtbc = MtbcRandom(debug=True, list_length=30)
