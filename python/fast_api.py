from fastapi import FastAPI, Query, UploadFile, File, Request
import uvicorn
import json
from fastapi.responses import FileResponse
from fastapi.templating import Jinja2Templates
from mtbc_package import mtbc_ncbi, mtbc_tools
import uuid
import os

templates = Jinja2Templates(directory="templates")

test = FastAPI()


@test.get("/")
async def root(request: Request):
    fasta = os.listdir("alignement/")
    nj_tree = os.listdir("nj_tree/")
    return templates.TemplateResponse("index.j2", {"request": request, "fasta": fasta, "nj_tree": nj_tree})


@test.get("/download_fasta")
async def download_fasta(fasta: str = ''):
    return FileResponse(path='alignement/{0}'.format(fasta), media_type='text/plain',
                        filename="{0}.fasta".format(fasta))

@test.get("/download_nj_tree")
async def download_fasta(nj_tree: str = ''):
    return FileResponse(path='nj_tree/{0}'.format(nj_tree), media_type='text/plain',
                        filename="{0}.nwk".format(nj_tree))



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
                      email: str = 'A.N.Other@example.com'
                      ):
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
    return FileResponse(path='alignement/{0}'.format(mtbc_fasta.id), media_type='text/plain',
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
    return FileResponse(path='nj_tree/{0}'.format(mtbc_fasta.id), media_type='text/plain', filename="{0}.nwk".format(mtbc_fasta.id))


if __name__ == "__main__":
    uvicorn.run(test, host="0.0.0.0", port=80)
