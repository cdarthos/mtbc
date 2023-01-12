from __future__ import annotations
import json
import logging
from datetime import datetime
from types import SimpleNamespace
from typing import Union, List
import uvicorn
from fastapi import FastAPI, Request, Query
from fastapi.templating import Jinja2Templates
from pymongo import MongoClient
from starlette import status
from starlette.responses import Response, RedirectResponse
from mtbc_package import mtbc_ncbi, mtbc_tools, mtbc_tree
from settings import mongoSettings
import time

FORMAT = "%(levelname)s:%(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

mongosettings = mongoSettings()

templates = Jinja2Templates(directory="templates")
logging.info("Start fastapi")

tags_metadata = [
    {
        "name": "mtbc_sra_list",
        "description": "Select parameters",
    },
    {
        "name": "id_interface",
        "description": "Go to basic html interface summarize id result",
    },
    {
        "name": "interface",
        "description": "Go to basic html interface summarize store result",
    },
    {
        "name": "download_sra",
        "description": "Get json file with all store information",
    },
    {
        "name": "download_fasta",
        "description": "Get fasta file alignment",
    },
    {
        "name": "download_fasta",
        "description": "Get fasta file alignment",
    }
]

test = FastAPI(openapi_tags=tags_metadata)


@test.get("/old")
async def root(request: Request):
    client = None
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        db_sra_list = request_data.find().distinct("_id")
        db_nj_tree = request_data.find({"nj_tree": {"$ne": None}}).distinct("_id")
        db_ml_tree = request_data.find({"ml_tree": {"$ne": None}}).distinct("_id")
        db_fasta = request_data.find({"fasta": {"$ne": None}}).distinct("_id")
        tests = request_data.find({}, {"_id": 1, "fasta": 1, "nj_tree": 1, "ml_tree": 1, "final_acc_list_length": 1,
                                       "create_date": 1})
        test2 = list(tests)
        logging.debug(test2)
    finally:
        client.close()
    return templates.TemplateResponse("index.j2",
                                      {"request": request, "fasta": db_fasta, "nj_tree": db_nj_tree,
                                       "sra_list": db_sra_list,
                                       "ml_tree": db_ml_tree,
                                       "test": test2})


@test.get("/", tags=["interface"])
async def test_root(request: Request):
    client = None
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data

        tests = request_data.find({}, {"_id": 1,
                                       "create_date": 1,
                                       "final_acc_list_length": 1,
                                       "target_list_length": 1,                                       
                                       "snp_select": 1,
                                       "snp_reject": 1,
                                       "email": 1})
        tests = list(tests)
        logging.info(tests)
    finally:
        client.close()
    return templates.TemplateResponse("test.j2",
                                      {"request": request,
                                       "test": tests})

@test.get("/id/{id}", tags=["id_interface"])
async def id_interface(id: str,request: Request):
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result_cursor = request_data.find_one({"_id": id})
        result = list(result_cursor)
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()
    logging.info(result)
    logging.info(result_cursor["_id"])

    fasta = False
    if result_cursor["fasta"] :
        fasta = True
    ml = False
    if result_cursor["ml_tree"] :
        ml = True
    nj = False
    if result_cursor["nj_tree"] :
        nj = True

    final_length = result_cursor["final_acc_list_length"]
    target_length = result_cursor["target_list_length"]
    snp_select = result_cursor["snp_select"]
    snp_reject = result_cursor["snp_reject"]
    create_date = None
    if create_date in result_cursor:
        create_date = result_cursor["create_date"]


    logging.info(final_length)


    return templates.TemplateResponse("id_interface.j2",
                                      {"request": request,
                                       "id": id,
                                       "fasta": fasta,
                                       "ml": ml,
                                       "nj": nj,
                                       "final_length":final_length,
                                       "target_length":target_length,
                                       "snp_select":snp_select,
                                       "snp_reject":snp_reject,
                                       "create_date":create_date})



@test.get("/download_sra/{id}", tags=["download_sra"])
async def download_sra(id: str = None):
    result = None
    client = None
    if id is None:
        try:
            client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
        finally:
            client.close()
        return request_data.find().distinct("_id")

    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result = request_data.find_one({"_id": id})
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()

    return result


@test.get("/download_fasta/{id}", tags=["download_fasta"])
async def download_fasta(id: Union[str, None] = None):
    result = None
    client = None
    if id is None:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        client.close()
        return request_data.find({"fasta": {"$ne": None}}).distinct("_id")

    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result = request_data.find_one({"_id": id}, {"_id": 0, "fasta": 1})["fasta"]
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()

    response = Response(result, media_type='text/plain')
    response.headers["Content-Disposition"] = 'attachment; filename={0}.fasta'.format(id)
    return response


@test.get("/download_nj_tree/{id}")
async def download_nj_tree(id: Union[str, None] = None):
    client = None
    result = None
    if id is None:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        client.close()
        return request_data.find({"nj_tree": {"$ne": None}}).distinct("_id")

    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result = request_data.find_one({"_id": id})["nj_tree"]
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()

    response = Response(result, media_type='text/plain')
    response.headers["Content-Disposition"] = 'attachment; filename={0}-nj.nwk'.format(id)
    return response


@test.get("/download_ml_tree/{id}")
async def download_ml_tree(id: Union[str, None] = None):
    result = None
    client = None
    if id is None:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        client.close()
        return request_data.find({"ml_tree": {"$ne": None}}).distinct("_id")

    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result = request_data.find_one({"_id": id})["ml_tree"]
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()

    response = Response(result, media_type='text/plain')
    response.headers["Content-Disposition"] = 'attachment; filename={0}-ml.nwk'.format(id)
    return response


@test.post("/mtbc_sra_list", tags=["mtbc_sra_list"])
async def set_param(select_mycobacterium_canettii: bool = False,
                    select_mycobacterium_mungi: bool = False,
                    select_mycobacterium_orygis: bool = False,
                    select_mycobacterium_tuberculosis: bool = False,
                    # outgroup: str = '',
                    all_id_to_acc: bool = False,

                    email: str = 'A.N.Other@example.com',
                    snp_select: Union[List[str], None] = Query(default=None),
                    snp_reject: Union[List[str], None] = Query(default=None),
                    # raxml_parameter: Union[List[str], None] = Query(default=None),
                    target_list_length: int = 100
                    ):
    client = None
    start_time = time.time()
    mtbc_inst = mtbc_ncbi.MtbcGetRandomSRA(select_mycobacterium_canettii=select_mycobacterium_canettii,
                                           select_mycobacterium_mungi=select_mycobacterium_mungi,
                                           select_mycobacterium_orygis=select_mycobacterium_orygis,
                                           select_mycobacterium_tuberculosis=select_mycobacterium_tuberculosis,
                                           # outgroup=outgroup,
                                           email=email,
                                           snp_select=snp_select,
                                           snp_reject=snp_reject,
                                           target_list_length=target_list_length,
                                           all_id_to_acc=all_id_to_acc
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
        # request_data.update_one({"_id": get_id}, {"$set": {"raxml_parameter": raxml_parameter}})
        request_data.update_one({"_id": get_id}, {"$set": {"mtbc_sra_list_time": mtbc_sra_list_time}})
        request_data.update_one({"_id": get_id}, {"$set": {"create_date": datetime.utcnow()}})
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()
    #return mtbc_inst.to_json()
    url = "/id/{0}".format(get_id)
    response = RedirectResponse(url=url, status_code=status.HTTP_302_FOUND)
    return response


@test.get("/mtbc_fasta_align_from_json")
@test.get("/mtbc_fasta_align_from_json/{id}")
def fasta_align_from_json(id: str = ""):
    result = None
    client = None
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result = request_data.find_one({"_id": id})
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()

    start_time = time.time()
    result_json = json.dumps(result)
    result_obj = None
    logging.info(result_obj._id)
    mtbc_fasta = mtbc_tools.MtbcAcclistToFASTA(id = id,
                                               ncbi_random_acc_list = result["ncbi_random_acc_list"],                                               
                                               sequence_dict=result["sequence_dict"],
                                               target_list_length=result["target_list_length"],
                                               final_acc_list=result["final_acc_list"],
                                               snp_reject=result["snp_reject"],
                                               snp_select=result["snp_select"])
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
    client = None
    fasta = None
    result = None
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result = request_data.find_one({"_id": id})["nj_tree"]
        fasta = request_data.find_one({"_id": id})["fasta"]
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()

    if result is not None:
        return result
    else:
        if fasta is None:
            fasta_ = fasta_align_from_json(id)
            logging.info(fasta_)
            fasta = fasta_.body.decode("utf-8")
            logging.debug(fasta)

        start_time = time.time()
        nj_tree = mtbc_tree.MtbcTree.create_nj_tree_static(fasta)
        nj_tree_time = time.time() - start_time

        try:
            client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
            request_data.update_one(
                {"_id": id},
                {"$set": {"nj_tree": nj_tree}})
            request_data.update_one(
                {"_id": id},
                {"$set": {"nj_tree_time": nj_tree_time}})
        except Exception:
            logging.error("error to connect mongo db")
        finally:
            client.close()
        return RedirectResponse("/download_nj_tree/{0}".format(id))


@test.get("/mtbc_ml_tree_from_db/{id}")
async def ml_tree_from_db(id: str = ""):
    client = None
    fasta = None
    result = None
    try:
        client = MongoClient('mongodb://{0}:{1}/'.format(mongosettings.host, mongosettings.port))
        db_mtbc = client.db_mtbc
        request_data = db_mtbc.request_data
        result = request_data.find_one({"_id": id})["ml_tree"]
        fasta = request_data.find_one({"_id": id})["fasta"]
    except Exception:
        logging.error("error to connect mongo db")
    finally:
        client.close()
    if result is not None:
        return result
    else:
        if fasta is None:
            fasta_ = fasta_align_from_json(id)
            logging.info(fasta_)
            fasta = fasta_.body.decode("utf-8")
            logging.debug(fasta)

        start_time = time.time()
        ml_tree = mtbc_tree.MtbcTree.create_ml_tree_static(fasta)
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
        except Exception:
            logging.error("error to connect mongo db")
        finally:
            client.close()
        return RedirectResponse("/download_ml_tree/{0}".format(id))


if __name__ == "__main__":
    uvicorn.run(test, host="0.0.0.0", port=8080)
