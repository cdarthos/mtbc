import os

import flask
from mtbc import MtbcRandom

app = flask.Flask(__name__)


@app.route("/")
def hello_world():
    return """<form method="POST">
  <div>
    <label for="say">What greeting do you want to say?</label>
    <input name="say" id="say" value="Hi" />
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


@app.route("/mtbc")
def flask_mtbc():
    tree = MtbcRandom(debug=True, list_length=30)

    return os.listdir('tree')

if __name__ == "__main__":
    mtbc = MtbcRandom(debug=True, list_length=30)