import json
import logging
import os
from types import SimpleNamespace

import uvicorn
from fastapi import FastAPI, Request
from fastapi.templating import Jinja2Templates
from pymongo import MongoClient
from starlette.responses import Response, RedirectResponse
from mtbc_package import mtbc_ncbi, mtbc_tools
from settings import mongoSettings


FORMAT = "%(levelname)s:%(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

mongosettings = mongoSettings()

templates = Jinja2Templates(directory="templates")
#logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)
logging.info("Start fastapi")
test = FastAPI()


@test.get("/")
async def root(request: Request):
    client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
    db_mtbc = client.db_mtbc
    request_data = db_mtbc.request_data
    db_sra_list = request_data.find().distinct("_id")

    db_nj_tree = request_data.find({"nj_tree": {"$ne": None}}).distinct("_id")

    db_ml_tree = request_data.find({"ml_tree": {"$ne": None}}).distinct("_id")

    db_fasta = request_data.find({"fasta": {"$ne": None}}).distinct("_id")

    client.close()
    return templates.TemplateResponse("index.j2",
                                      {"request": request, "fasta": db_fasta, "nj_tree": db_nj_tree,
                                       "sra_list": db_sra_list,
                                       "ml_tree": db_ml_tree})


@test.get("/get_sra_list_form")
async def get_sra_list_form(request: Request):
    return templates.TemplateResponse("sra_list_form.j2",
                                      {"request": request}
                                      )


@test.get("/download_sra/{id}")
async def download_sra(id: str = ''):
    client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
    db_mtbc = client.db_mtbc
    request_data = db_mtbc.request_data

    resultat = request_data.find_one({"_id": id})
    client.close()

    return resultat


@test.get("/download_fasta/{id}")
async def download_fasta(id: str = ''):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
    except:
        logging.error("error to connect mongo db")

    resultat = request_data.find_one({"_id": id})["fasta"]
    client.close()

    response = Response(resultat, media_type='text/plain')
    response.headers["Content-Disposition"] = 'attachment; filename={0}.fasta'.format(id)
    return response


@test.get("/download_nj_tree/{id}")
async def download_nj_tree(id: str = ''):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
    except:
        logging.error("error to connect mongo db")

    resultat = request_data.find_one({"_id": id})["nj_tree"]
    client.close()

    response = Response(resultat, media_type='text/plain')
    response.headers["Content-Disposition"] = 'attachment; filename={0}-nj.nwk'.format(id)
    return response


@test.get("/download_ml_tree/{id}")
async def download_ml_tree(id: str = ""):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
    except:
        logging.error("error to connect mongo db")

    resultat = request_data.find_one({"_id": id})["ml_tree"]
    client.close()

    response = Response(resultat, media_type='text/plain')
    response.headers["Content-Disposition"] = 'attachment; filename={0}-ml.nwk'.format(id)
    return response


@test.post("/mtbc_sra_list")
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
    client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
    db_mtbc = client.db_mtbc
    request_data = db_mtbc.request_data
    get_id = request_data.insert_one(mtbc_inst.to_json()).inserted_id
    logging.info(get_id)
    client.close()
    return mtbc_inst.to_json()


@test.get("/mtbc_fasta_align_from_json")
async def fasta_align_from_json(id: str = ""):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
    except:
        logging.error("error to connect mongo db")
    resultat = request_data.find_one({"_id": id})
    client.close()

    resultat_json = json.dumps(resultat)
    resultat_obj = json.loads(resultat_json, object_hook=lambda d: SimpleNamespace(**d))
    logging.info(resultat_obj._id)

    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(resultat_obj)

    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
    except:
        logging.error("error to connect mongo db")
    # request_data.update_one(
    #     {"_id": id},
    #     {"$set": {"fasta": mtbc_fasta.fasta}})
    # logging.info("""request_data.update_one(
    #     {"_id": id},
    #     {"$set": {"fasta": mtbc_fasta.fasta}})""")
    request_data.update_one(
       {"_id": id},
       {"$set": {"sequence_dict": mtbc_fasta.sequence_dict}})
    logging.info("""request_data.update_one(
       {"_id": id},
       {"$set": {"sequence_dict": mtbc_fasta.sequence_dict}})""")


    client.close()
    return mtbc_fasta.fasta


@test.get("/mtbc_nj_tree_from_db/{id}")
async def nj_tree_from_db(id: str = ""):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
    except:
        logging.error("error to connect mongo db")
    resultat = request_data.find_one({"_id": id})["nj_tree"]
    fasta = request_data.find_one({"_id": id})["fasta"]
    client.close()

    if resultat is not None:
        return resultat
    else:
        if fasta is None:

            fasta_coroutine = fasta_align_from_json(id)
            await fasta_coroutine
            try:
                client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
                db_mtbc = client.db_mtbc
                request_data = db_mtbc.request_data
            except:
                logging.error("error to connect mongo db")
            fasta = request_data.find_one({"_id": id})["fasta"]
            client.close()

        nj_tree = mtbc_tools.MtbcTree.create_nj_tree_static(id, fasta)
        try:
            client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
        except:
            logging.error("error to connect mongo db")
        request_data.update_one(
            {"_id": id},
            {"$set": {"nj_tree": nj_tree}})
        client.close()
        return RedirectResponse("/download_nj_tree/{0}".format(id))


@test.get("/mtbc_ml_tree_from_db/{id}")
async def ml_tree_from_db(id: str = ""):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
    except:
        logging.error("error to connect mongo db")
    resultat = request_data.find_one({"_id": id})["ml_tree"]
    fasta = request_data.find_one({"_id": id})["fasta"]
    client.close()
    if resultat is not None:
        return resultat
    else:
        if fasta is None:
            fasta_coroutine = fasta_align_from_json(id)
            await fasta_coroutine
            try:
                client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
                db_mtbc = client.db_mtbc
                request_data = db_mtbc.request_data
            except:
                logging.error("error to connect mongo db")
            fasta = request_data.find_one({"_id": id})["fasta"]
            client.close()

        ml_tree = mtbc_tools.MtbcTree.create_ml_tree_static(id, fasta)
        try:
            client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
        except:
            logging.error("error to connect mongo db")
        request_data.update_one(
            {"_id": id},
            {"$set": {"ml_tree": ml_tree}})
        client.close()
        return RedirectResponse("/download_ml_tree/{0}".format(id))


if __name__ == "__main__":
    uvicorn.run(test, host="0.0.0.0", port=8080)
