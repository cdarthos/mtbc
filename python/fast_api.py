from fastapi import FastAPI
import uvicorn
from python.mtbc import mtbc_ncbi as MT

app = FastAPI()


@app.get("/")
async def root():
    return {"message": "Hello test"}

@app.get("/mtbc")
async def set_param():
    mtbc_inst = mtbc.MtbcGetRandomSRA()
    return mtbc_inst.to_json()

if __name__ == "__main__":
    # mtbc = MtbcGetRandomSRA(debug=True, list_length=1000,select_mycobacterium_canettii=True,
    # select_mycobacterium_mungi=True, select_mycobacterium_orygis=True)
    mtbc = MT.MtbcGetRandomSRA(debug=True, list_length=10000, id1=2)
