import os

from pydantic import BaseSettings


class mongoSettings(BaseSettings):
    host: str = os.environ['MONGO_HOST']
    port: int = 27017