from fastapi import FastAPI
import uvicorn
from python.mtbc_package import mtbc_ncbi

test = FastAPI()


@test.get("/")
async def root():
    return {"message": "Hello test"}


@test.get("/mtbc")
async def set_param(debug: bool = True):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(debug=debug)
    return mtbc_inst.to_json()


if __name__ == "__main__":
    uvicorn.run(test, host="0.0.0.0", port=8000)
