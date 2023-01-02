from pydantic import BaseSettings


class mongoSettings(BaseSettings):
    host: str = "localhost"
    port: int = 27017