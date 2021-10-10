# David R Thompson
import pymongo
from emit_main.config.config import Config
import os, sys
from datetime import datetime
import argparse


def main():

    parser = argparse.ArgumentParser(description='Estimate the dark state at a time.')
    parser.add_argument('timestamp', metavar='T', type=str, help='timestamp of the form %Y%m%d%H%M%S')
    parser.add_argument('config_path', metavar='P', type=str, 
            help='path to confiuration file', 
            default='../emit-main/emit_main/config/test_sds_config.json')
    args = parser.parse_args()

    config_path = args.config_path

    # Get config properties
    config = Config(config_path).get_dictionary()
    client = pymongo.MongoClient(config["mongodb_host"], config["mongodb_port"])
    db = client[config["mongodb_db_name"]]

    print('connected')
    acquisitions = db.acquisitions
    dt = datetime.strptime(args.timestamp,'%Y%m%d%H%M%S')
    isodate = dt.strftime('%Y-%m-%dT%H:%M:%S.000+00:00')
    print(isodate)
    higher = acquisitions.find({"submode": 'science', 'start_time':{'$gte':dt}})
    higher = higher.sort([("start_time",pymongo.DESCENDING)])
    higher = higher.limit(1)
    for h in higher:
       print('---')
       print(h)
    print('==================================')

    lower = acquisitions.find({"submode": 'science', 'start_time':{'$lte':dt}})
    lower = lower.sort([("start_time",pymongo.ASCENDING)])
    lower = lower.limit(1)
    for l in lower:
       print('---')
       print(l)
    print('---')

if __name__ == '__main__':
    main()
