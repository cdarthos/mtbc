FROM python:3.9

RUN apt-get update && apt install -y raxml


WORKDIR /code

COPY requirements.txt /code/requirements.txt

RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

RUN pip install --no-cache-dir pymongo

#COPY . /code/


WORKDIR /code/python

CMD ["uvicorn", "fast_api:test", "--host", "0.0.0.0", "--port", "80", "--workers", "4", "--log-level", "debug"]

# If running behind a proxy like Nginx or Traefik add --proxy-headers
# CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--proxy-headers"]