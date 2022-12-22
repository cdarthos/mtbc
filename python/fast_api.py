from fastapi import FastAPI
import uvicorn
import json
from fastapi.responses import FileResponse

from mtbc_package import mtbc_ncbi, mtbc_tools
import uuid



test = FastAPI()


@test.get("/")
async def root():
    return {"message": "Hello test"}


@test.get("/mtbc_sra_list")
async def set_param(debug: bool = False,
                    select_mycobacterium_canettii: bool = False,
                    select_mycobacterium_mungi: bool = False,
                    select_mycobacterium_orygis: bool = False,
                    select_mycobacterium_tuberculosis: bool = False,
                    outgroup: str = '',
                    retmax: int = 10000000,
                    list_length: int = 10,
                    email: str = 'A.N.Other@example.com'):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(debug=debug,
                                           select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           retmax=retmax,
                                           list_length=list_length,
                                           email=email)
    print(json.dumps(mtbc_inst.to_json()))
    return mtbc_inst.to_json()


@test.get("/mtbc_fasta_align")
async def fasta_align(debug: bool = False,
                      select_mycobacterium_canettii: bool = False,
                      select_mycobacterium_mungi: bool = False,
                      select_mycobacterium_orygis: bool = False,
                      select_mycobacterium_tuberculosis: bool = False,
                      outgroup: str = '',
                      retmax: int = 10000000,
                      list_length: int = 10,
                      email: str = 'A.N.Other@example.com',
                      id1: int = 1):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(debug=debug,
                                           select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           retmax=retmax,
                                           list_length=list_length,
                                           email=email)
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(mtbc_inst)
    mtbc_fasta.align_reconstruct()
    return FileResponse(path='alignement/{0}.fasta'.format(mtbc_fasta.id), media_type='text/plain',
                        filename="{0}.fasta".format(mtbc_fasta.id))


@test.get("/mtbc_nj_tree")
async def nj_tree(debug: bool = False,
                  select_mycobacterium_canettii: bool = False,
                  select_mycobacterium_mungi: bool = False,
                  select_mycobacterium_orygis: bool = False,
                  select_mycobacterium_tuberculosis: bool = False,
                  outgroup: str = '',
                  retmax: int = 10000000,
                  list_length: int = 10,
                  email: str = 'A.N.Other@example.com'):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(debug=debug,
                                           select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           retmax=retmax,
                                           list_length=list_length,
                                           email=email)
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(mtbc_inst)
    mtbc_fasta.align_reconstruct()
    mtbc_fasta.create_nj_tree()
    return FileResponse(path='tree/{0}.nwk'.format(mtbc_fasta.id), media_type='text/plain', filename="{0}".nwk)


if __name__ == "__main__":
    uvicorn.run(test, host="0.0.0.0", port=8000)