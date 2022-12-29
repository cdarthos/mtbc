from fastapi import FastAPI, Query, UploadFile, File, Request
import uvicorn
import json
from fastapi.responses import FileResponse
from fastapi.templating import Jinja2Templates
from uvicorn import logging

from mtbc_package import mtbc_ncbi, mtbc_tools, custom_encoder
import uuid
import os
from types import SimpleNamespace
import logging

templates = Jinja2Templates(directory="templates")

test = FastAPI()

logger = logging.getLogger()
logger.setLevel(logging.INFO)

test.json_sra = os.listdir("request/")


@test.get("/")
async def root(request: Request):
    sra_list = os.listdir("request/")
    fasta = os.listdir("alignement/")
    nj_tree = os.listdir("nj_tree/")
    ml_tree = os.listdir("ml_tree/")
    return templates.TemplateResponse("index.j2", {"request": request, "fasta": fasta, "nj_tree": nj_tree, "sra_list": sra_list, "ml_tree": ml_tree})

@test.get("/download_sra")
async def download_sra(id: str = ''):
    with open("request/{0}".format(id), 'r') as json_request:
        return json.load(json_request)


@test.get("/download_fasta")
async def download_fasta(id: str = ''):
    return FileResponse(path='alignement/{0}'.format(id), media_type='text/plain',
                        filename="{0}.fasta".format(id))

@test.get("/download_nj_tree")
async def download_nj_tree(id: str = ''):
    return FileResponse(path='nj_tree/{0}'.format(id), media_type='text/plain',
                        filename="{0}.nwk".format(id))

@test.get("/download_ml_tree")
async def download_ml_tree(id: str = ""):
    return FileResponse(path='ml_tree/{0}'.format(id), media_type='text/plain',
                        filename="{0}.nwk".format(id))

@test.get("/mtbc_sra_list")
async def set_param(select_mycobacterium_canettii: bool = False,
                    select_mycobacterium_mungi: bool = False,
                    select_mycobacterium_orygis: bool = False,
                    select_mycobacterium_tuberculosis: bool = False,
                    outgroup: str = '',
                    list_length: int = 10,
                    email: str = 'A.N.Other@example.com'):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           list_length=list_length,
                                           email=email)

    logging.info("LOGGING")
    test.json_sra = os.listdir("request/")
    logger.info(test.json_sra)
    return mtbc_inst.to_json()


@test.get("/mtbc_fasta_align")
async def fasta_align(select_mycobacterium_canettii: bool = False,
                      select_mycobacterium_mungi: bool = False,
                      select_mycobacterium_orygis: bool = False,
                      select_mycobacterium_tuberculosis: bool = False,
                      outgroup: str = '',
                      list_length: int = 10,
                      email: str = 'A.N.Other@example.com'
                      ):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           list_length=list_length,
                                           email=email)
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(mtbc_inst)
    return FileResponse(path='alignement/{0}'.format(mtbc_fasta.id), media_type='text/plain',
                        filename="{0}.fasta".format(mtbc_fasta.id))


@test.get("/mtbc_fasta_align_from_json")
async def fasta_align_from_json(json_file: str = Query("id",enum=os.listdir("request/"))):
    with open("request/{0}".format(json_file), 'r') as json_request:
        mtbc_json = json.load(json_request, object_hook=lambda d: SimpleNamespace(**d))
    logger.info(mtbc_json)
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(mtbc_json)
    return FileResponse(path='alignement/{0}'.format(mtbc_fasta.id), media_type='text/plain',filename="{0}.fasta".format(mtbc_fasta.id))

@test.get("/mtbc_nj_tree_from_fasta")
async def nj_tree_from_fasta(json_file: str = Query("id",enum=os.listdir("request/"))):
    with open("request/{0}".format(json_file), 'r') as json_request:
        mtbc_json = json.load(json_request, object_hook=lambda d: SimpleNamespace(**d))
    logger.info(mtbc_json)
    nj_tree = mtbc_tools.MtbcTree(mtbc_json)
    nj_tree.create_nj_tree()
    return FileResponse(path='nj_tree/{0}'.format(nj_tree.id), media_type='text/plain',
                        filename="{0}-nj.nwk".format(nj_tree.id))


@test.get("/mtbc_ml_tree_from_fasta")
async def ml_tree_from_fasta(json_file: str = Query("id",enum=os.listdir("request/"))):
    with open("request/{0}".format(json_file), 'r') as json_request:
        mtbc_json = json.load(json_request, object_hook=lambda d: SimpleNamespace(**d))
    logger.info(mtbc_json)
    ml_tree = mtbc_tools.MtbcTree(mtbc_json)
    ml_tree.create_ml_tree()
    return FileResponse(path='ml_tree/{0}'.format(ml_tree.id), media_type='text/plain',
                        filename="{0}-ml.nwk".format(ml_tree.id))









@test.get("/mtbc_nj_tree")
async def nj_tree(select_mycobacterium_canettii: bool = False,
                  select_mycobacterium_mungi: bool = False,
                  select_mycobacterium_orygis: bool = False,
                  select_mycobacterium_tuberculosis: bool = False,
                  outgroup: str = '',
                  list_length: int = 10,
                  email: str = 'A.N.Other@example.com'):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           list_length=list_length,
                                           email=email)
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(mtbc_inst)
    mtbc_fasta.align_reconstruct()
    mtbc_fasta.create_nj_tree()
    return FileResponse(path='nj_tree/{0}-nj'.format(mtbc_fasta.id), media_type='text/plain', filename="{0}.nwk".format(mtbc_fasta.id))

@test.get("/mtbc_ml_tree")
async def ml_tree(select_mycobacterium_canettii: bool = False,
                  select_mycobacterium_mungi: bool = False,
                  select_mycobacterium_orygis: bool = False,
                  select_mycobacterium_tuberculosis: bool = False,
                  outgroup: str = '',
                  list_length: int = 10,
                  email: str = 'A.N.Other@example.com'):
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           list_length=list_length,
                                           email=email)
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(mtbc_inst)
    mtbc_fasta.align_reconstruct()
    mtbc_fasta.create_ml_tree()
    return FileResponse(path='ml_tree/{0}'.format(mtbc_fasta.id), media_type='text/plain',
                        filename="{0}-ml.nwk".format(mtbc_fasta.id))


if __name__ == "__main__":
    uvicorn.run(test, host="0.0.0.0", port=80)
