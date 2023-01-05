import json
import logging
import os
from types import SimpleNamespace

import uvicorn
from fastapi import FastAPI, Request
from fastapi.templating import Jinja2Templates
from pymongo import MongoClient
from starlette.responses import Response, RedirectResponse
from uvicorn.loops import asyncio

from mtbc_package import mtbc_ncbi, mtbc_tools, mtbc_tree
from settings import mongoSettings
import time

FORMAT = "%(levelname)s:%(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

mongosettings = mongoSettings()

templates = Jinja2Templates(directory="templates")
# logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
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

    tests = request_data.find({}, {"_id": 1, "fasta": 1, "nj_tree": 1, "ml_tree": 1, "final_acc_list_length": 1})
    # for test in tests:
    #    logging.info(test)
    # logging.info("\n")

    test2 = list(tests)
    logging.debug(test2)

    client.close()
    return templates.TemplateResponse("index.j2",
                                      {"request": request, "fasta": db_fasta, "nj_tree": db_nj_tree,
                                       "sra_list": db_sra_list,
                                       "ml_tree": db_ml_tree,
                                       "test": test2})


@test.get("/get_sra_list_form")
async def get_sra_list_form(request: Request):
    return templates.TemplateResponse("sra_list_form.j2",
                                      {"request": request}
                                      )


@test.get("/download_sra/{id}")
async def download_sra(id: str = ''):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        resultat = request_data.find_one({"_id": id})
    except:
        logging.error("error to connect mongo db")
    finally:
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

    resultat = request_data.find_one({"_id": id}, {"_id": 0, "fasta": 1})["fasta"]
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
                    # ncbi_list_length: int = 100,

                    email: str = 'A.N.Other@example.com',
                    snp_select: list = [],
                    snp_reject: list = [],
                    target_list_length=100
                    ):
    start_time = time.time()
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           outgroup=outgroup,
                                           # ncbi_list_length=ncbi_list_length,
                                           ncbi_list_length=3 * target_list_length,
                                           email=email,
                                           snp_select=snp_select,
                                           snp_reject=snp_reject,
                                           target_list_length=target_list_length
                                           )
    logging.info("LOGGING")
    mtbc_sra_list_time = time.time() - start_time
    logging.info("mtbc_sra_list_time : " + str(mtbc_sra_list_time))
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        get_id = request_data.insert_one(mtbc_inst.to_json()).inserted_id
        logging.info("get_id : " + str(get_id))
        request_data.update_one({"_id": get_id}, {"$set": {"mtbc_sra_list_time": mtbc_sra_list_time}})
    except:
        logging.error("error to connect mongo db")
    finally:
        client.close()
    return mtbc_inst.to_json()


@test.get("/mtbc_fasta_align_from_json")
def fasta_align_from_json(id: str = ""):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        resultat = request_data.find_one({"_id": id})
    except:
        logging.error("error to connect mongo db")
    finally:
        client.close()

    start_time = time.time()

    resultat_json = json.dumps(resultat)
    resultat_obj = json.loads(resultat_json, object_hook=lambda d: SimpleNamespace(**d))
    logging.info(resultat_obj._id)
    sequence_dict = resultat["sequence_dict"]
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(resultat_obj, sequence_dict,
                                               target_list_length=resultat["target_list_length"],
                                               final_acc_list=resultat["final_acc_list"])
    logging.info("Preparation envoi à MongoDB")

    mtbc_fasta_align_from_json_time = time.time() - start_time
    logging.info("mtbc_fasta_align_from_json : " + str(mtbc_fasta_align_from_json_time))

    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        request_data.update_one(
            {"_id": id},
            {"$set": {"mtbc_fasta_align_from_json_time": mtbc_fasta_align_from_json_time}})
        request_data.update_one(
            {"_id": id},
            {"$set": {"final_acc_list_length": mtbc_fasta.final_acc_list_length}})
        request_data.update_one(
            {"_id": id},
            {"$set": {"final_acc_list": mtbc_fasta.final_acc_list}})
        request_data.update_one(
            {"_id": id},
            {"$set": {"fasta": mtbc_fasta.fasta}})
        logging.info("""request_data.update_one(
            {"_id": id},
            {"$set": {"fasta": mtbc_fasta.fasta}})""")
        request_data.update_one(
            {"_id": id},
            {"$set": {"sequence_dict": mtbc_fasta.sequence_dict}})
        logging.info("""request_data.update_one(
            {"_id": id},
            {"$set": {"sequence_dict": mtbc_fasta.sequence_dict}})""")
    except Exception as e:
        logging.error("error to connect mongo db")
        logging.error(e)
    finally:
        client.close()
        logging.info("Mongodb connection close")

    response = Response(mtbc_fasta.fasta, media_type='text/plain')
    response.headers["Content-Disposition"] = 'attachment; filename={0}.fasta'.format(id)

    return response


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

            fasta_ = fasta_align_from_json(id)
            logging.info(fasta_)
            fasta = fasta_.body.decode("utf-8")
            logging.debug(fasta)
            #await fasta_coroutine
            #try:
            #    client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
            #    db_mtbc = client.db_mtbc
            #    request_data = db_mtbc.request_data
            #except:
            #    logging.error("error to connect mongo db")
            #fasta = request_data.find_one({"_id": id})["fasta"]
            #client.close()

        start_time = time.time()
        nj_tree = mtbc_tree.MtbcTree.create_nj_tree_static(id, fasta)
        nj_tree_time = time.time() - start_time

        try:
            client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
        except:
            logging.error("error to connect mongo db")
        request_data.update_one(
            {"_id": id},
            {"$set": {"nj_tree": nj_tree}})
        request_data.update_one(
            {"_id": id},
            {"$set": {"nj_tree_time": nj_tree_time}})
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

        start_time = time.time()
        ml_tree = mtbc_tree.MtbcTree.create_ml_tree_static(id, fasta)
        ml_tree_time = time.time() - start_time

        try:
            client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
            request_data.update_one(
                {"_id": id},
                {"$set": {"ml_tree": ml_tree}})
            request_data.update_one(
                {"_id": id},
                {"$set": {"ml_tree_time": ml_tree_time}})
        except:
            logging.error("error to connect mongo db")
        finally:
            client.close()
        return RedirectResponse("/download_ml_tree/{0}".format(id))


if __name__ == "__main__":
    uvicorn.run(test, host="0.0.0.0", port=8080)
